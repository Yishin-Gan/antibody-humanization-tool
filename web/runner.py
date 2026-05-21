"""
Single-submission pipeline runner.

Wraps pipeline.generate_sequence + evaluation.score_sequences so the Flask
web app can submit one mouse VH+VL pair and any combination of processing
modes (pipeline, preferred, sapiens, lab) and get back a single dict
ready for the interactive report.

Heavy work (ANARCI, Sapiens, OASis 9-mer DB, CamSol) lives in the existing
modules — this file just chooses which steps to run and stitches the
per-mode results together.
"""

from __future__ import annotations

import sys
import os
import io
import contextlib
import time
import traceback
import hashlib
import ast
from dataclasses import dataclass, field
from typing import Optional

sys.path.insert(0, "/workspace/antibody-humanization-tool")

from pipeline.step_a_numbering import number_sequence
from pipeline.step_b_germline_scoring import (
    rank_germlines, _GERMLINE_FR_DB, normalize_germline_name,
)
from pipeline.generate_sequence import graft, humanize_sapiens, get_germline_seq_by_region
from evaluation.score_sequences import score_one


# Spec mode → seq_id used inside score_one. score_one branches on these
# IDs for germline-identity logic (seq0 detects post-hoc, seq2/seq6 also
# compute post-Sapiens drift, seq5 uses the canonical pathway).
MODE_TO_SID = {
    "pipeline":  "2",
    "preferred": "6",
    "sapiens":   "0",
    "lab":       "5",
}
ALL_MODES = ["pipeline", "preferred", "sapiens", "lab"]


@dataclass
class RunResult:
    """In-memory bundle for one submission. Indexed by job_id in the Flask app."""
    job_id:           str
    timestamp:        float
    mouse_vh:         str
    mouse_vl:         str
    vl_chain_type:    str
    active_modes:     list                          = field(default_factory=list)
    pipeline_germ_vh: Optional[str]                 = None
    pipeline_germ_vl: Optional[str]                 = None
    preferred_germ_vh: Optional[str]                = None
    preferred_germ_vl: Optional[str]                = None
    lab_hu_vh:        str                           = ""
    lab_hu_vl:        str                           = ""
    lab_final_vh:     str                           = ""
    lab_final_vl:     str                           = ""
    lab_germ_vh:      Optional[str]                 = None
    lab_germ_vl:      Optional[str]                 = None
    structure:        bool                          = False
    # per-mode outputs: { mode -> {vh, vl, vh_germline, vl_germline, scores_row, grafted_vh, grafted_vl} }
    modes:            dict                          = field(default_factory=dict)
    # per-mode germline residues from get_germline_seq_by_region: { mode -> {VH: {region: {pos: aa}}, VL: ...} }
    germ_seq:         dict                          = field(default_factory=dict)
    # mouse-scored row (for OASis FR bar in summary)
    mouse_scores:     dict                          = field(default_factory=dict)
    error:            Optional[str]                 = None
    log:              list                          = field(default_factory=list)


# ── input sanitisation ────────────────────────────────────────────────────────

_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean_aa(raw: str) -> str:
    """Strip FASTA headers/whitespace/numbers, uppercase, keep only valid AAs."""
    if not raw:
        return ""
    out = []
    for line in raw.splitlines():
        if line.startswith(">"):
            continue
        for ch in line:
            if ch.isalpha():
                up = ch.upper()
                if up in _VALID_AA:
                    out.append(up)
    return "".join(out)


def validate_germline_name(name: str, chain: str) -> dict:
    """Look up a germline name in _GERMLINE_FR_DB.

    Returns {found: bool, resolved: <name in DB or None>, note: str}.
    chain is 'H', 'K', or 'L'. If name has *XX suffix, try exact; otherwise
    fall back to first matching gene-level allele.
    """
    if not name:
        return {"found": False, "resolved": None, "note": ""}
    db = _GERMLINE_FR_DB.get(chain, {})
    if name in db:
        return {"found": True, "resolved": name, "note": "exact match"}
    if "*" not in name:
        # gene-level — find first allele
        for entry in db:
            if entry.split("*")[0] == name:
                return {"found": True, "resolved": entry,
                        "note": f"resolved gene-level → {entry}"}
    # try *01 fallback
    fallback = f"{name.split('*')[0]}*01"
    if fallback in db:
        return {"found": True, "resolved": fallback,
                "note": f"allele not in DB → using {fallback}"}
    return {"found": False, "resolved": None,
            "note": "not in database"}


def detect_vl_chain_type(vl_seq: str) -> dict:
    """Run ANARCI on a VL sequence and report the inferred chain type.

    Returns {chain_type: 'K'|'L'|None, ok: bool, error: str|None}.
    """
    seq = clean_aa(vl_seq)
    if len(seq) < 60:
        return {"chain_type": None, "ok": False, "error": "VL too short to number"}
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            numbered = number_sequence(seq, chain_type=None)
        ct = numbered.get("chain_type")
        return {"chain_type": ct, "ok": ct in ("K", "L"), "error": None}
    except Exception as e:
        return {"chain_type": None, "ok": False, "error": str(e)}


# ── main runner ───────────────────────────────────────────────────────────────


def _make_job_id(mouse_vh: str, mouse_vl: str) -> str:
    h = hashlib.sha1(f"{mouse_vh}|{mouse_vl}|{time.time()}".encode()).hexdigest()
    return h[:12]


def _log(result: RunResult, msg: str) -> None:
    print(msg, flush=True)
    result.log.append(msg)


def run(
    mouse_vh:           str,
    mouse_vl:           str,
    modes:              list,
    vl_chain_type:      Optional[str] = None,
    preferred_germ_vh:  Optional[str] = None,
    preferred_germ_vl:  Optional[str] = None,
    lab_hu_vh:          str = "",
    lab_hu_vl:          str = "",
    lab_final_vh:       str = "",
    lab_final_vl:       str = "",
    lab_germ_vh:        Optional[str] = None,
    lab_germ_vl:        Optional[str] = None,
    run_structure:      bool = False,
) -> RunResult:
    """Run the requested subset of the humanization pipeline on a single pair.

    `modes` is a subset of ALL_MODES. At least one must be chosen.
    Lab reference mode requires lab_hu_vh, lab_hu_vl, lab_final_vh, lab_final_vl
    and lab_germ_vh, lab_germ_vl (so the lab seq5 can be scored against its
    grafted baseline).

    Returns a RunResult with raw sequences and per-mode score rows. Reporting
    helpers in report_data.py convert this into the JSON the report consumes.
    """
    mouse_vh = clean_aa(mouse_vh)
    mouse_vl = clean_aa(mouse_vl)
    lab_hu_vh = clean_aa(lab_hu_vh)
    lab_hu_vl = clean_aa(lab_hu_vl)
    lab_final_vh = clean_aa(lab_final_vh)
    lab_final_vl = clean_aa(lab_final_vl)

    result = RunResult(
        job_id=_make_job_id(mouse_vh, mouse_vl),
        timestamp=time.time(),
        mouse_vh=mouse_vh, mouse_vl=mouse_vl,
        vl_chain_type=vl_chain_type or "",
        active_modes=list(modes),
        preferred_germ_vh=preferred_germ_vh,
        preferred_germ_vl=preferred_germ_vl,
        lab_hu_vh=lab_hu_vh, lab_hu_vl=lab_hu_vl,
        lab_final_vh=lab_final_vh, lab_final_vl=lab_final_vl,
        lab_germ_vh=lab_germ_vh, lab_germ_vl=lab_germ_vl,
        structure=run_structure,
    )

    try:
        # ── Step A: number mouse and infer chain type if not provided ─────────
        _log(result, "[stage] Numbering mouse VH/VL")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            vh_num = number_sequence(mouse_vh, chain_type="H")
            vl_num = number_sequence(mouse_vl, chain_type=vl_chain_type or None)
        ct = vl_chain_type or vl_num.get("chain_type") or "K"
        result.vl_chain_type = ct
        _log(result, f"[stage] VL chain type: {ct}")

        # ── Step B: pipeline germline ranking (always — also used to display
        #    the pipeline-detected germline even when only Sapiens mode runs) ──
        _log(result, "[stage] Ranking germlines")
        with contextlib.redirect_stdout(io.StringIO()):
            vh_ranks = rank_germlines(vh_num["fr_residues"], "H", top_n=5)
            vl_ranks = rank_germlines(vl_num["fr_residues"], ct, top_n=5)
        pipe_vh = vh_ranks[0]["germline"] if vh_ranks else None
        pipe_vl = vl_ranks[0]["germline"] if vl_ranks else None
        result.pipeline_germ_vh = pipe_vh
        result.pipeline_germ_vl = pipe_vl
        _log(result, f"[stage] Pipeline germline VH={pipe_vh} VL={pipe_vl}")

        # Resolve preferred germlines if requested
        if "preferred" in modes:
            if preferred_germ_vh:
                v = validate_germline_name(preferred_germ_vh, "H")
                result.preferred_germ_vh = v["resolved"] or preferred_germ_vh
            if preferred_germ_vl:
                v = validate_germline_name(preferred_germ_vl, ct)
                result.preferred_germ_vl = v["resolved"] or preferred_germ_vl

        # Resolve / detect lab germlines for lab mode
        if "lab" in modes:
            # If user gave lab germlines, normalise + validate. Otherwise detect.
            if lab_germ_vh:
                v = validate_germline_name(lab_germ_vh, "H")
                result.lab_germ_vh = v["resolved"] or lab_germ_vh
            else:
                with contextlib.redirect_stdout(io.StringIO()):
                    rk = rank_germlines(
                        number_sequence(lab_hu_vh, chain_type="H")["fr_residues"],
                        "H", top_n=1) if lab_hu_vh else []
                result.lab_germ_vh = rk[0]["germline"] if rk else None
            if lab_germ_vl:
                v = validate_germline_name(lab_germ_vl, ct)
                result.lab_germ_vl = v["resolved"] or lab_germ_vl
            else:
                with contextlib.redirect_stdout(io.StringIO()):
                    rk = rank_germlines(
                        number_sequence(lab_hu_vl, chain_type=ct)["fr_residues"],
                        ct, top_n=1) if lab_hu_vl else []
                result.lab_germ_vl = rk[0]["germline"] if rk else None

        # ── For each active mode, build sequences ─────────────────────────────
        per_mode = {}

        if "pipeline" in modes:
            _log(result, "[stage] Pipeline mode: graft + Sapiens")
            with contextlib.redirect_stdout(io.StringIO()):
                vh1, _ = graft(mouse_vh, pipe_vh, chain_type="H") if pipe_vh else (None, None)
                vl1, _ = graft(mouse_vl, pipe_vl, chain_type=ct) if pipe_vl else (None, None)
                vh2, _ = humanize_sapiens(vh1, "H") if vh1 else (None, None)
                vl2, _ = humanize_sapiens(vl1, ct) if vl1 else (None, None)
            per_mode["pipeline"] = {
                "vh": vh2, "vl": vl2,
                "vh_germline": pipe_vh, "vl_germline": pipe_vl,
                "grafted_vh": vh1, "grafted_vl": vl1,
            }

        if "preferred" in modes:
            _log(result, "[stage] Preferred mode: graft + Sapiens")
            pv = result.preferred_germ_vh
            pl = result.preferred_germ_vl
            with contextlib.redirect_stdout(io.StringIO()):
                vh4, _ = graft(mouse_vh, pv, chain_type="H") if pv else (None, None)
                vl4, _ = graft(mouse_vl, pl, chain_type=ct) if pl else (None, None)
                vh6, _ = humanize_sapiens(vh4, "H") if vh4 else (None, None)
                vl6, _ = humanize_sapiens(vl4, ct) if vl4 else (None, None)
            per_mode["preferred"] = {
                "vh": vh6, "vl": vl6,
                "vh_germline": pv, "vl_germline": pl,
                "grafted_vh": vh4, "grafted_vl": vl4,
            }

        if "sapiens" in modes:
            _log(result, "[stage] Sapiens mode: Sapiens on mouse")
            with contextlib.redirect_stdout(io.StringIO()):
                vh0, _ = humanize_sapiens(mouse_vh, "H")
                vl0, _ = humanize_sapiens(mouse_vl, ct)
            per_mode["sapiens"] = {
                "vh": vh0, "vl": vl0,
                "vh_germline": None, "vl_germline": None,
                "grafted_vh": mouse_vh, "grafted_vl": mouse_vl,  # seq0 baseline = mouse
            }

        if "lab" in modes:
            _log(result, "[stage] Lab mode: using supplied sequences")
            per_mode["lab"] = {
                "vh": lab_final_vh, "vl": lab_final_vl,
                "vh_germline": result.lab_germ_vh, "vl_germline": result.lab_germ_vl,
                # seq5 in the existing code uses seq3 (lab_hu) as the grafted baseline
                "grafted_vh": lab_hu_vh, "grafted_vl": lab_hu_vl,
            }

        # ── Score each mode using score_one ───────────────────────────────────
        for mode, data in per_mode.items():
            sid = MODE_TO_SID[mode]
            _log(result, f"[stage] Scoring mode={mode} (seq_id={sid})")
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    row = score_one(
                        clone_id="query",
                        seq_id=sid,
                        seq_label=mode,
                        vh_seq=data["vh"] or "",
                        vl_seq=data["vl"] or "",
                        vl_chain_type=ct,
                        mouse_vh=mouse_vh,
                        mouse_vl=mouse_vl,
                        lab_hu_vh=lab_hu_vh, lab_hu_vl=lab_hu_vl,
                        lab_final_vh=lab_final_vh, lab_final_vl=lab_final_vl,
                        vh_germline=data["vh_germline"] or "",
                        vl_germline=data["vl_germline"] or "",
                        grafted_vh=data["grafted_vh"] or "",
                        grafted_vl=data["grafted_vl"] or "",
                        run_structure=run_structure,
                    )
                data["scores_row"] = row
            except Exception as e:
                _log(result, f"[error] scoring {mode} failed: {e}")
                data["scores_row"] = {"clone": "query", "seq_id": sid, "_error": str(e)}

        # ── Score the mouse input separately (for OASis FR baseline bar) ──────
        _log(result, "[stage] Scoring mouse baseline")
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mouse_row = score_one(
                    clone_id="query", seq_id="m", seq_label="mouse",
                    vh_seq=mouse_vh, vl_seq=mouse_vl,
                    vl_chain_type=ct,
                    mouse_vh=mouse_vh, mouse_vl=mouse_vl,
                    lab_hu_vh="", lab_hu_vl="",
                    lab_final_vh="", lab_final_vl="",
                    vh_germline="", vl_germline="",
                    grafted_vh="", grafted_vl="",
                    run_structure=False,
                )
            result.mouse_scores = mouse_row
        except Exception as e:
            _log(result, f"[error] mouse baseline scoring failed: {e}")
            result.mouse_scores = {"_error": str(e)}

        result.modes = per_mode

        # ── Retrieve germline residue sets per active mode for alignment panels
        for mode, data in per_mode.items():
            try:
                gvh = get_germline_seq_by_region(data["vh_germline"], "H") if data["vh_germline"] else {}
                gvl = get_germline_seq_by_region(data["vl_germline"], ct) if data["vl_germline"] else {}
                result.germ_seq[mode] = {"VH": gvh, "VL": gvl}
            except Exception as e:
                result.germ_seq[mode] = {"VH": {}, "VL": {}, "_error": str(e)}

        _log(result, "[stage] Done")
    except Exception as e:
        result.error = f"{e}\n{traceback.format_exc()}"
        _log(result, f"[fatal] {e}")
    return result
