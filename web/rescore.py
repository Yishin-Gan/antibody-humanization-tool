"""
Re-score a sequence after the user has applied residue-level edits in
the report's tweak mode.

Re-runs only the cheap-ish components the spec calls out:
    - OASis FR/CDR identity   (compute_oasis_per_position + paired compute_oasis)
    - CamSol intrinsic scores (compute_camsol)
    - Germline FR identity    (compute_germline_identity)
    - Vernier / FR mutable classification for affected positions
    - Liabilities (cheap, useful when edits add/remove NG, M, DG, etc.)

Does NOT re-run Sapiens, ANARCI numbering (length unchanged → reuse cached
positions), or ABodyBuilder2.
"""

from __future__ import annotations

import ast
import contextlib
import io
import sys
from typing import Optional

sys.path.insert(0, "/workspace/antibody-humanization-tool")

from pipeline.step_a_numbering import number_sequence
from evaluation.score_sequences import (
    compute_oasis_per_position, compute_oasis,
    compute_germline_identity, compute_camsol,
    compute_vernier, find_liabilities, compute_physicochemical,
    compute_structure_scores,
    AA_CHARGE_PH7, _estimate_pi,
)


def apply_edits(seq: str, edits: dict) -> str:
    """Apply a {linear_index_1based: new_aa} dict to a sequence string.

    Linear index is the 1-based residue index in the displayed sequence (not
    IMGT) — this matches how the front-end emits the edit positions.
    """
    if not edits:
        return seq
    chars = list(seq)
    for k, aa in edits.items():
        i = int(k) - 1
        if 0 <= i < len(chars) and aa and aa.upper() in "ACDEFGHIKLMNPQRSTVWY":
            chars[i] = aa.upper()
    return "".join(chars)


def rescore_structure(vh_seq: str, vl_seq: str, clone_id: str) -> dict:
    """Re-run ABodyBuilder2 paired-Fv structure prediction. ~30–60s typical.

    Returns the conf_* fields ready to patch into the report payload.
    """
    if not vh_seq or not vl_seq:
        return {}
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            return compute_structure_scores(vh_seq, vl_seq, clone_id)
        except Exception as e:
            return {"_error": f"structure prediction failed: {e}"}


def rescore_chain(
    edited_seq:    str,
    mouse_seq:     str,
    chain_type:    str,
    germline_name: str,
    grafted_seq:   str,
) -> dict:
    """Re-score one chain. Returns the subset of fields the report consumes.

    Pure pipeline functions — no DB/network beyond OASis SQLite and the
    in-process CamSol predictor.
    """
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            num = number_sequence(edited_seq, chain_type=chain_type)
        except Exception as e:
            return {"_error": f"numbering failed: {e}"}

        mouse_num = number_sequence(mouse_seq, chain_type=chain_type) if mouse_seq else None
        grafted_num = number_sequence(grafted_seq, chain_type=chain_type) if grafted_seq else None

        oasis_res = compute_oasis_per_position(edited_seq, num, chain_type)
        # Build query oasis FR map (for Vernier consistency)
        try:
            pp_detail = ast.literal_eval(oasis_res.get("oasis_per_position_detail", "[]"))
            oasis_fr = {d["imgt_pos"]: d["nmer"][0] for d in pp_detail
                        if "FR" in d.get("region", "")}
        except Exception:
            pp_detail = []
            oasis_fr = None

        germline_id, _ = (None, None)
        if germline_name:
            germline_id, _ = compute_germline_identity(num, germline_name, chain_type)

        camsol = compute_camsol(edited_seq, num)
        liab = find_liabilities(edited_seq, num)
        physico = compute_physicochemical(edited_seq)

        vern = {}
        if mouse_num:
            try:
                mouse_oasis = compute_oasis_per_position(mouse_seq, mouse_num, chain_type)
                mouse_pp = ast.literal_eval(mouse_oasis.get("oasis_per_position_detail", "[]"))
                mouse_fr = {d["imgt_pos"]: d["nmer"][0] for d in mouse_pp}
            except Exception:
                mouse_fr = None
            try:
                graft_oasis = compute_oasis_per_position(grafted_seq, mouse_num, chain_type) if grafted_seq else None
                if graft_oasis is not None:
                    g_pp = ast.literal_eval(graft_oasis.get("oasis_per_position_detail", "[]"))
                    graft_fr = {d["imgt_pos"]: d["nmer"][0] for d in g_pp}
                else:
                    graft_fr = None
            except Exception:
                graft_fr = None
            try:
                vern = compute_vernier(num, mouse_num, grafted_num, chain_type,
                                       oasis_fr, mouse_fr, graft_fr)
            except Exception:
                vern = {}

    pp_map = {}
    for d in pp_detail:
        pos = d.get("imgt_pos")
        nmer = d.get("nmer", "")
        if pos is not None and nmer:
            pp_map[str(pos)] = {"aa": nmer[0], "region": d.get("region", "")}

    backmut_detail = []
    raw_bm = vern.get("backmut_detail", "[]")
    if isinstance(raw_bm, str):
        try:
            raw_bm = ast.literal_eval(raw_bm)
        except Exception:
            raw_bm = []
    for d in raw_bm or []:
        backmut_detail.append({
            "pos":     int(d.get("imgt_pos")),
            "vernier": bool(d.get("is_vernier")),
            "mouse":   d.get("mouse_aa"),
            "grafted": d.get("grafted_aa"),
            "query":   d.get("query_aa"),
            "status":  d.get("status"),
        })

    return {
        "seq":          edited_seq,
        "oasis":        round(oasis_res.get("oasis_identity") or 0.0, 4),
        "oasis_fr":     round(oasis_res.get("oasis_fr_identity") or 0.0, 4),
        "oasis_cdr":    round(oasis_res.get("oasis_cdr_identity") or 0.0, 4),
        "germline_id":  round(germline_id, 4) if germline_id is not None else None,
        "camsol":       round(camsol.get("camsol_score") or 0.0, 4),
        "camsol_fr":    round(camsol.get("camsol_fr_score") or 0.0, 4),
        "camsol_cdr":   round(camsol.get("camsol_cdr_score") or 0.0, 4),
        "hs_fr":        camsol.get("camsol_hotspot_fr_count"),
        "hs_cdr":       camsol.get("camsol_hotspot_cdr_count"),
        "vern_mut":     vern.get("vernier_mutable_count"),
        "vern_back":    vern.get("vernier_backmut_count"),
        "vern_hum":     vern.get("vernier_humanized_count"),
        "vern_other":   vern.get("vernier_other_count"),
        "fr_mut":       vern.get("fr_mutable_count"),
        "fr_back":      vern.get("fr_backmut_count"),
        "fr_hum":       vern.get("fr_humanized_count"),
        "fr_other":     vern.get("fr_other_count"),
        "pp":           pp_map,
        "mutable":      backmut_detail,
        "lia": {
            "dc": liab.get("deamidation_cdr_count"),
            "df": liab.get("deamidation_fr_count"),
            "oc": liab.get("oxidation_cdr_count"),
            "of": liab.get("oxidation_fr_count"),
            "if": liab.get("isomerization_fr_count"),
        },
        "pi":           round(_estimate_pi(edited_seq), 2),
        "ch":           round(sum(AA_CHARGE_PH7.get(a, 0) for a in edited_seq), 2),
    }
