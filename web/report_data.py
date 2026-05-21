"""
Convert a RunResult into the JSON the interactive report template consumes.

The report is fully client-rendered; this module produces one large
serialisable dict (DATA) that is dumped as JSON into the page.
"""

from __future__ import annotations

import ast
import contextlib
import io
import json
import sys
from typing import Optional

sys.path.insert(0, "/workspace/antibody-humanization-tool")
from pipeline.step_a_numbering import number_sequence

from web.runner import RunResult


# IMGT canonical region boundaries
IMGT_REGIONS = {
    "FR1":  (1, 26),
    "CDR1": (27, 38),
    "FR2":  (39, 55),
    "CDR2": (56, 65),
    "FR3":  (66, 104),
    "CDR3": (105, 117),
    "FR4":  (118, 128),
}

VERNIER_VH = {2, 27, 29, 30, 47, 48, 67, 69, 71, 78, 80, 93, 94}
VERNIER_VL = {2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71}

# Colour assignments per mode (spec §4.2)
MODE_LABEL = {
    "pipeline":  "Pipeline",
    "preferred": "Preferred",
    "sapiens":   "Sapiens",
    "lab":       "Lab ref",
}


def _f(v, ndigits: Optional[int] = 6):
    try:
        if v is None or v == "" or v == "None":
            return None
        f = float(v)
        return round(f, ndigits) if ndigits is not None else f
    except (ValueError, TypeError):
        return None


def _i(v):
    try:
        if v is None or v == "" or v == "None":
            return None
        return int(float(v))
    except (ValueError, TypeError):
        return None


def _safe_eval(s):
    if not s or s in ("[]", "{}", "None"):
        return []
    try:
        return ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return []


def find_regions(seq: str, cdr1: str, cdr2: str, cdr3: str) -> list:
    """Substring-locate CDRs and return [{t,start,end}] for FR/CDR/FR... segments.

    Falls back to fuzzy match (≤3 mismatches) when CDR string isn't found
    exactly — matches the existing build_report.py behaviour.
    """
    if not seq:
        return []
    regions = []
    pos = 0
    for cs, cn in [(cdr1, "CDR1"), (cdr2, "CDR2"), (cdr3, "CDR3")]:
        if not cs:
            continue
        idx = seq.find(cs, pos)
        if idx == -1:
            best_i, best_d = -1, len(cs)
            for i in range(pos, len(seq) - len(cs) + 1):
                d = sum(1 for a, b in zip(seq[i:i+len(cs)], cs) if a != b)
                if d < best_d:
                    best_i, best_d = i, d
            if best_d <= 3:
                idx = best_i
            else:
                continue
        if idx > pos:
            regions.append({"t": "FR", "start": pos, "end": idx})
        regions.append({"t": cn, "start": idx, "end": idx + len(cs)})
        pos = idx + len(cs)
    if pos < len(seq):
        regions.append({"t": "FR", "start": pos, "end": len(seq)})
    return regions


def _fr_label_for(pos: int) -> str:
    if pos <= 26: return "FR1"
    if pos <= 55: return "FR2"
    if pos <= 104: return "FR3"
    return "FR4"


def _cdr_label_for(pos: int) -> str:
    if pos <= 38: return "CDR1"
    if pos <= 65: return "CDR2"
    return "CDR3"


def _abnumber_chain(seq: str):
    """Return an abnumber Chain numbered with IMGT, or None on failure.

    Important: we must use abnumber (not ANARCI's number_sequence) because the
    per-position OASis data in scores.csv is also abnumber-numbered. ANARCI
    and abnumber can assign DIFFERENT IMGT positions to the same residue when
    CDR1 has non-standard length (e.g. 8C11 VL), causing alignment offsets.
    Using abnumber consistently for the linear→IMGT map and the scaffold
    residue lookup guarantees both agree with the OASis per-position data.
    """
    if not seq:
        return None
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            from abnumber import Chain
            return Chain(seq, scheme="imgt", cdr_definition="imgt")
    except Exception:
        return None


def _scaffold_imgt_residues(seq: str, chain_type: str) -> dict:
    """Return {str(imgt_pos): {aa, region}} for the grafted scaffold sequence,
    numbered by abnumber so it matches the mode sequence's numbering.

    Spec §7.3 says the alignment "Germline" row shows the scaffold sequence
    used during grafting (seq1 for pipeline, seq4 for preferred, seq3 for
    lab) — NOT the canonical DB residues.
    """
    c = _abnumber_chain(seq)
    if c is None:
        return {}
    out = {}
    for pos, aa in c:
        base = pos.number  # insertion-coded positions still expose .number
        # Determine region from position label if available, else infer
        region = None
        try:
            region = pos.get_region()  # 'FR1' .. 'FR4', 'CDR1'..'CDR3'
        except Exception:
            pass
        if not region:
            if base <= 26: region = "FR1"
            elif base <= 38: region = "CDR1"
            elif base <= 55: region = "FR2"
            elif base <= 65: region = "CDR2"
            elif base <= 104: region = "FR3"
            elif base <= 117: region = "CDR3"
            else: region = "FR4"
        out[str(base)] = {"aa": aa, "region": region}
    return out


def _linear_to_imgt(seq: str, chain_type: str) -> list:
    """Return list[int|None] of length len(seq): IMGT position per linear residue.

    Uses abnumber to match the OASis per-position numbering. Insertion-coded
    positions (e.g. CDR3 112A) collapse to the base IMGT integer.
    """
    c = _abnumber_chain(seq)
    if c is None:
        return [None] * len(seq)
    out = []
    for pos, aa in c:
        out.append(int(pos.number))
    if len(out) > len(seq): out = out[:len(seq)]
    while len(out) < len(seq): out.append(None)
    return out


def _per_pos_map(detail_list):
    """Build {imgt_pos: {aa, region}} from a per-position OASis detail list."""
    out = {}
    for d in detail_list:
        pos = d.get("imgt_pos")
        if pos is None:
            continue
        nmer = d.get("nmer", "")
        if not nmer:
            continue
        out[str(pos)] = {"aa": nmer[0], "region": d.get("region", "")}
    return out


def _backmut_list(detail_list):
    """Convert backmut_detail records to JSON-safe form, sorted by IMGT pos."""
    out = []
    for d in detail_list:
        out.append({
            "pos":      int(d.get("imgt_pos")),
            "vernier":  bool(d.get("is_vernier")),
            "mouse":    d.get("mouse_aa"),
            "grafted":  d.get("grafted_aa"),
            "query":    d.get("query_aa"),
            "status":   d.get("status"),
        })
    out.sort(key=lambda x: x["pos"])
    return out


def _germ_db_from_detail(detail_list):
    """Build {imgt_pos: {aa, region}} from vh/vl_germline_fr_detail."""
    out = {}
    for pos, aa in detail_list:
        out[str(pos)] = {"aa": aa, "region": "FR"}
    return out


def _scores_block(sc: dict, mode: str) -> dict:
    """Extract the metric subset the report consumes from a scores row."""
    sid = {"pipeline": "2", "preferred": "6", "sapiens": "0", "lab": "5"}[mode]
    block = {
        # OASis
        "oa_vh":      _f(sc.get("oasis_vh_identity")),
        "oa_vl":      _f(sc.get("oasis_vl_identity")),
        "oa_fr_vh":   _f(sc.get("vh_oasis_fr_identity")),
        "oa_fr_vl":   _f(sc.get("vl_oasis_fr_identity")),
        "oa_cdr_vh":  _f(sc.get("vh_oasis_cdr_identity")),
        "oa_cdr_vl":  _f(sc.get("vl_oasis_cdr_identity")),
        # Germline identity
        "g_vh":       _f(sc.get("vh_germline_identity")),
        "g_vl":       _f(sc.get("vl_germline_identity")),
        # CamSol
        "cs_vh":      _f(sc.get("vh_camsol_score")),
        "cs_vl":      _f(sc.get("vl_camsol_score")),
        "cs_fr_vh":   _f(sc.get("vh_camsol_fr_score")),
        "cs_fr_vl":   _f(sc.get("vl_camsol_fr_score")),
        "cs_cdr_vh":  _f(sc.get("vh_camsol_cdr_score")),
        "cs_cdr_vl":  _f(sc.get("vl_camsol_cdr_score")),
        "hs_vh":      _i(sc.get("vh_camsol_hotspot_count")),
        "hs_vl":      _i(sc.get("vl_camsol_hotspot_count")),
        "hs_fr_vh":   _i(sc.get("vh_camsol_hotspot_fr_count")),
        "hs_fr_vl":   _i(sc.get("vl_camsol_hotspot_fr_count")),
        "hs_cdr_vh":  _i(sc.get("vh_camsol_hotspot_cdr_count")),
        "hs_cdr_vl":  _i(sc.get("vl_camsol_hotspot_cdr_count")),
        # Vernier counts
        "vh_vern_mut":       _i(sc.get("vh_vernier_mutable_count")),
        "vh_vern_back":      _i(sc.get("vh_vernier_backmut_count")),
        "vh_vern_hum":       _i(sc.get("vh_vernier_humanized_count")),
        "vh_vern_other":     _i(sc.get("vh_vernier_other_count")),
        "vl_vern_mut":       _i(sc.get("vl_vernier_mutable_count")),
        "vl_vern_back":      _i(sc.get("vl_vernier_backmut_count")),
        "vl_vern_hum":       _i(sc.get("vl_vernier_humanized_count")),
        "vl_vern_other":     _i(sc.get("vl_vernier_other_count")),
        # All FR counts
        "vh_fr_mut":         _i(sc.get("vh_fr_mutable_count")),
        "vh_fr_back":        _i(sc.get("vh_fr_backmut_count")),
        "vh_fr_hum":         _i(sc.get("vh_fr_humanized_count")),
        "vh_fr_other":       _i(sc.get("vh_fr_other_count")),
        "vl_fr_mut":         _i(sc.get("vl_fr_mutable_count")),
        "vl_fr_back":        _i(sc.get("vl_fr_backmut_count")),
        "vl_fr_hum":         _i(sc.get("vl_fr_humanized_count")),
        "vl_fr_other":       _i(sc.get("vl_fr_other_count")),
        # Physicochemical
        "pi":         _f(sc.get("fv_pi"), 2),
        "ch":         _f(sc.get("fv_net_charge_ph7"), 2),
        # Liabilities
        "lia": {
            "vh_dc": _i(sc.get("vh_deamidation_cdr_count")),
            "vh_df": _i(sc.get("vh_deamidation_fr_count")),
            "vh_oc": _i(sc.get("vh_oxidation_cdr_count")),
            "vh_of": _i(sc.get("vh_oxidation_fr_count")),
            "vh_if": _i(sc.get("vh_isomerization_fr_count")),
            "vl_dc": _i(sc.get("vl_deamidation_cdr_count")),
            "vl_df": _i(sc.get("vl_deamidation_fr_count")),
            "vl_oc": _i(sc.get("vl_oxidation_cdr_count")),
            "vl_of": _i(sc.get("vl_oxidation_fr_count")),
            "vl_if": _i(sc.get("vl_isomerization_fr_count")),
        },
        # Structure (if computed)
        "cf_vh":       _f(sc.get("conf_mean_vh"), 4),
        "cf_vl":       _f(sc.get("conf_mean_vl"), 4),
        "cf_fr_vh":    _f(sc.get("conf_fr_mean_vh"), 4),
        "cf_fr_vl":    _f(sc.get("conf_fr_mean_vl"), 4),
        "cf_cdr_vh":   _f(sc.get("conf_cdr_mean_vh"), 4),
        "cf_cdr_vl":   _f(sc.get("conf_cdr_mean_vl"), 4),
        "cf_cdr1_vh":  _f(sc.get("conf_cdr1_mean_vh"), 4),
        "cf_cdr2_vh":  _f(sc.get("conf_cdr2_mean_vh"), 4),
        "cf_cdr3_vh":  _f(sc.get("conf_cdr3_mean_vh"), 4),
        "cf_cdr1_vl":  _f(sc.get("conf_cdr1_mean_vl"), 4),
        "cf_cdr2_vl":  _f(sc.get("conf_cdr2_mean_vl"), 4),
        "cf_cdr3_vl":  _f(sc.get("conf_cdr3_mean_vl"), 4),
        "cf_min_vh":   _f(sc.get("conf_min_vh"), 4),
        "cf_min_vl":   _f(sc.get("conf_min_vl"), 4),
        # Drift (only present on sid 2, 6)
        "drift_vh_post":     sc.get(f"vh_detected_germline_post_sap_seq{sid}") or None,
        "drift_vh_post_id":  _f(sc.get(f"vh_detected_germline_post_sap_seq{sid}_identity")),
        "drift_vh_flag":     str(sc.get(f"vh_germline_drifted_seq{sid}") or "") == "True",
        "drift_vh_delta":    _f(sc.get(f"vh_germline_identity_delta_seq{sid}")),
        "drift_vl_post":     sc.get(f"vl_detected_germline_post_sap_seq{sid}") or None,
        "drift_vl_post_id":  _f(sc.get(f"vl_detected_germline_post_sap_seq{sid}_identity")),
        "drift_vl_flag":     str(sc.get(f"vl_germline_drifted_seq{sid}") or "") == "True",
        "drift_vl_delta":    _f(sc.get(f"vl_germline_identity_delta_seq{sid}")),
    }
    return block


def _flatten_germline_seq(germ_seq_regions: dict) -> dict:
    """Convert {region: {pos: aa}} → {pos: {aa, region}}."""
    out = {}
    for region, mapping in (germ_seq_regions or {}).items():
        for pos, aa in mapping.items():
            out[str(pos)] = {"aa": aa, "region": region}
    return out


def build(result: RunResult) -> dict:
    """Top-level: returns a fully JSON-serialisable dict for the report."""
    chain_type = result.vl_chain_type or "K"

    modes_out = {}
    mutable_out = {"VH": {"vern": [], "non": []}, "VL": {"vern": [], "non": []}}
    # mutable_out is filled per mode using a list-of-lists by position;
    # the report merges across modes via the per-mode mutable map below.

    # Per-mode mutable maps: {mode: {VH/VL: [{pos, status, query, ...}]}}
    mut_by_mode = {}
    for mode, d in result.modes.items():
        sc = d.get("scores_row", {})
        bm_vh = _safe_eval(sc.get("vh_backmut_detail", "[]"))
        bm_vl = _safe_eval(sc.get("vl_backmut_detail", "[]"))
        mut_by_mode[mode] = {
            "VH": _backmut_list(bm_vh),
            "VL": _backmut_list(bm_vl),
        }

    for mode, d in result.modes.items():
        sc = d.get("scores_row", {})
        vh = d.get("vh") or ""
        vl = d.get("vl") or ""
        regions_vh = find_regions(vh,
                                  sc.get("vh_cdr1_sequence", ""),
                                  sc.get("vh_cdr2_sequence", ""),
                                  sc.get("vh_cdr3_sequence", ""))
        regions_vl = find_regions(vl,
                                  sc.get("vl_cdr1_sequence", ""),
                                  sc.get("vl_cdr2_sequence", ""),
                                  sc.get("vl_cdr3_sequence", ""))
        vh_pp_detail = _safe_eval(sc.get("vh_oasis_per_position_detail", "[]"))
        vl_pp_detail = _safe_eval(sc.get("vl_oasis_per_position_detail", "[]"))
        vh_pp = _per_pos_map(vh_pp_detail)
        vl_pp = _per_pos_map(vl_pp_detail)

        vh_germ_db = _germ_db_from_detail(_safe_eval(sc.get("vh_germline_fr_detail", "[]")))
        vl_germ_db = _germ_db_from_detail(_safe_eval(sc.get("vl_germline_fr_detail", "[]")))

        # Per-spec §7.3: alignment "Germline" row shows the scaffold sequence
        # used during grafting, not the canonical DB residues. For pipeline /
        # preferred / lab modes the scaffold is the grafted seq (seq1 / seq4 /
        # seq3). For sapiens mode there is no grafted scaffold — fall back to
        # the post-hoc detected germline's canonical DB residues.
        grafted_vh = d.get("grafted_vh") or ""
        grafted_vl = d.get("grafted_vl") or ""
        if mode == "sapiens":
            # Use canonical DB for the post-hoc detected germline
            germ_seq = result.germ_seq.get(mode, {})
            scaffold_vh = _flatten_germline_seq(germ_seq.get("VH", {}))
            scaffold_vl = _flatten_germline_seq(germ_seq.get("VL", {}))
        else:
            scaffold_vh = _scaffold_imgt_residues(grafted_vh, "H")
            scaffold_vl = _scaffold_imgt_residues(grafted_vl, chain_type)

        modes_out[mode] = {
            "label":         MODE_LABEL[mode],
            "vh":            vh,
            "vl":            vl,
            "vh_regions":    regions_vh,
            "vl_regions":    regions_vl,
            "vh_imgt":       _linear_to_imgt(vh, "H"),
            "vl_imgt":       _linear_to_imgt(vl, chain_type),
            "vh_pp":         vh_pp,
            "vl_pp":         vl_pp,
            "vh_germ_db":    vh_germ_db,     # canonical germline residues from DB
            "vl_germ_db":    vl_germ_db,
            "germ_seq_vh":   scaffold_vh,    # SCAFFOLD residues — per spec §7.3
            "germ_seq_vl":   scaffold_vl,
            "vh_germline":   d.get("vh_germline"),
            "vl_germline":   d.get("vl_germline"),
            "grafted_vh":    d.get("grafted_vh") or "",
            "grafted_vl":    d.get("grafted_vl") or "",
            "scores":        _scores_block(sc, mode),
            "mutable_vh":    mut_by_mode[mode]["VH"],
            "mutable_vl":    mut_by_mode[mode]["VL"],
            "cdr": {
                "vh1": sc.get("vh_cdr1_sequence", ""),
                "vh2": sc.get("vh_cdr2_sequence", ""),
                "vh3": sc.get("vh_cdr3_sequence", ""),
                "vl1": sc.get("vl_cdr1_sequence", ""),
                "vl2": sc.get("vl_cdr2_sequence", ""),
                "vl3": sc.get("vl_cdr3_sequence", ""),
            },
        }

    # Mouse baseline: minimal block — used for OASis FR bar and as reference
    mouse_sc = result.mouse_scores or {}
    mouse_regions_vh = find_regions(result.mouse_vh,
                                    mouse_sc.get("vh_cdr1_sequence", ""),
                                    mouse_sc.get("vh_cdr2_sequence", ""),
                                    mouse_sc.get("vh_cdr3_sequence", ""))
    mouse_regions_vl = find_regions(result.mouse_vl,
                                    mouse_sc.get("vl_cdr1_sequence", ""),
                                    mouse_sc.get("vl_cdr2_sequence", ""),
                                    mouse_sc.get("vl_cdr3_sequence", ""))
    mouse_pp_vh = _per_pos_map(_safe_eval(mouse_sc.get("vh_oasis_per_position_detail", "[]")))
    mouse_pp_vl = _per_pos_map(_safe_eval(mouse_sc.get("vl_oasis_per_position_detail", "[]")))
    mouse_out = {
        "label":      "Mouse",
        "vh":         result.mouse_vh,
        "vl":         result.mouse_vl,
        "vh_regions": mouse_regions_vh,
        "vl_regions": mouse_regions_vl,
        "vh_imgt":    _linear_to_imgt(result.mouse_vh, "H"),
        "vl_imgt":    _linear_to_imgt(result.mouse_vl, chain_type),
        "vh_pp":      mouse_pp_vh,
        "vl_pp":      mouse_pp_vl,
        "scores": {
            "oa_vh":     _f(mouse_sc.get("oasis_vh_identity")),
            "oa_vl":     _f(mouse_sc.get("oasis_vl_identity")),
            "oa_fr_vh":  _f(mouse_sc.get("vh_oasis_fr_identity")),
            "oa_fr_vl":  _f(mouse_sc.get("vl_oasis_fr_identity")),
            "oa_cdr_vh": _f(mouse_sc.get("vh_oasis_cdr_identity")),
            "oa_cdr_vl": _f(mouse_sc.get("vl_oasis_cdr_identity")),
            "cs_vh":     _f(mouse_sc.get("vh_camsol_score")),
            "cs_vl":     _f(mouse_sc.get("vl_camsol_score")),
            "cs_fr_vh":  _f(mouse_sc.get("vh_camsol_fr_score")),
            "cs_fr_vl":  _f(mouse_sc.get("vl_camsol_fr_score")),
            "pi":        _f(mouse_sc.get("fv_pi"), 2),
            "ch":        _f(mouse_sc.get("fv_net_charge_ph7"), 2),
        },
    }

    payload = {
        "job_id":       result.job_id,
        "timestamp":    result.timestamp,
        "chain_type":   chain_type,
        "active_modes": list(result.modes.keys()),
        "mouse_vh":     result.mouse_vh,
        "mouse_vl":     result.mouse_vl,
        "vernier_vh":   sorted(VERNIER_VH),
        "vernier_vl":   sorted(VERNIER_VL),
        "regions":      IMGT_REGIONS,
        "modes":        modes_out,
        "mouse":        mouse_out,
        "pipeline_germ_vh":  result.pipeline_germ_vh,
        "pipeline_germ_vl":  result.pipeline_germ_vl,
        "preferred_germ_vh": result.preferred_germ_vh,
        "preferred_germ_vl": result.preferred_germ_vl,
        "lab_germ_vh":       result.lab_germ_vh,
        "lab_germ_vl":       result.lab_germ_vl,
        "structure":         result.structure,
        "error":             result.error,
    }
    return payload


def to_json(result: RunResult) -> str:
    return json.dumps(build(result), default=str)
