"""
evaluate.py

Evaluates the pipeline against lab ground truth data.

Two evaluations:
  1. Germline selection — does the pipeline's top-ranked germline match
     the germline embedded in the lab's Hu sequence?
  2. Grafting accuracy — does the pipeline's grafted sequence match
     the lab's Hu sequence, and where do differences fall?

Input CSV columns:
    clone, mouse_vh, mouse_vl, hu_vh, hu_vl, final_vh, final_vl

Usage (from project root):
    python3 tests/evaluate.py --csv data/benchmarks/humanization_benchmark.csv
"""

# isort: skip_file
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E402
from abnumber import Chain
from pipeline.step_b_germline_scoring import rank_germlines, normalize_germline_name
from pipeline.step_a_numbering import number_sequence, IMGT_REGIONS, ALL_CDR_POSITIONS
from typing import Optional
from dataclasses import dataclass, field
from collections import Counter
import csv
import argparse

# ── Lab germline name parsing ─────────────────────────────────────────────────


def parse_lab_germline(name: str) -> str:
    """
    Normalize a lab-provided germline name to gene level.
    Handles:
      - Functional suffix: IGHV1-2*02F → IGHV1-2  (strips *allele and F suffix)
      - Standard allele:   IGHV3-23*01 → IGHV3-23  (strips *allele)
      - No allele:         IGHV3-23    → IGHV3-23   (unchanged)

    Returns gene-level name only (no allele) for consistent comparison.
    """
    import re
    name = name.strip()
    if not name:
        return name
    # Strip allele including any trailing functional/ORF suffix (F, P, ORF)
    name = re.sub(r'[*]\d+[A-Za-z]*$', '', name)
    return name.strip()


# IMGT CDR/FR boundary positions — differences here suggest CDR definition mismatch
CDR_BOUNDARY_POSITIONS = {26, 27, 38, 39, 55, 56, 65, 66, 104, 105, 117, 118}


# ── Data class for one clone's evaluation result ─────────────────────────────

@dataclass
class CloneEvalResult:
    clone_id: str

    # Germline selection (VH)
    vh_true_germline:             Optional[str] = None
    # FR identity of Hu vs detected germline
    vh_true_germline_fr_identity: Optional[float] = None
    # always from detect_germline(hu_vh)
    vh_detected_germline:        Optional[str] = None
    _hu_vh_det_rankings: list = field(
        default_factory=list)  # Hu-based detection rankings
    vh_pipeline_rank1:            Optional[str] = None
    vh_true_germline_rank:        Optional[int] = None
    vh_recall_at_1:               Optional[bool] = None
    vh_recall_at_3:               Optional[bool] = None
    vh_recall_at_5:               Optional[bool] = None

    # Germline selection (VL)
    vl_true_germline:             Optional[str] = None
    # FR identity of Hu vs detected germline
    vl_true_germline_fr_identity: Optional[float] = None
    # always from detect_germline(hu_vl)
    vl_detected_germline:        Optional[str] = None
    _hu_vl_det_rankings: list = field(
        default_factory=list)  # Hu-based detection rankings
    vl_pipeline_rank1:            Optional[str] = None
    vl_true_germline_rank:        Optional[int] = None
    vl_recall_at_1:               Optional[bool] = None
    vl_recall_at_3:               Optional[bool] = None
    vl_recall_at_5:               Optional[bool] = None

    # Grafting accuracy (VH)
    vh_grafted_seq:        Optional[str] = None
    vh_hu_seq:             Optional[str] = None
    vh_seq_match:          Optional[bool] = None
    vh_seq_identity:       Optional[float] = None
    vh_n_differences:      Optional[int] = None
    vh_diffs:              list = field(default_factory=list)

    # Grafting accuracy (VL)
    vl_grafted_seq:        Optional[str] = None
    vl_hu_seq:             Optional[str] = None
    vl_seq_match:          Optional[bool] = None
    vl_seq_identity:       Optional[float] = None
    vl_n_differences:      Optional[int] = None
    vl_diffs:              list = field(default_factory=list)

    error:                 Optional[str] = None

    # Debug data — populated when verbose=True, used for side-by-side output
    # Stored as dicts so CloneEvalResult stays serialisable
    _mouse_vh_numbered:    dict = field(default_factory=dict)
    _mouse_vl_numbered:    dict = field(default_factory=dict)
    _hu_vh_numbered:       dict = field(default_factory=dict)
    _hu_vl_numbered:       dict = field(default_factory=dict)
    _vh_rankings:          list = field(default_factory=list)
    _vl_rankings:          list = field(default_factory=list)
    _vl_chain_type:        str = "K"


# ── Helper: detect germline from Hu sequence using FR-only Step B scoring ────

def detect_germline(hu_sequence: str, chain_type: str) -> tuple[Optional[str], Optional[float], list]:
    """
    Detect the germline used in a CDR-grafted Hu sequence using FR-only identity.

    More reliable than ANARCI whole-sequence assignment because:
    - Hu sequences are chimeric (mouse CDRs + human framework)
    - ANARCI scores the full sequence, so mouse CDRs add noise to scoring
    - This approach extracts FR positions first (Step A), then scores FR-only
      against all germlines (Step B) — CDR noise never enters the comparison

    The detected germline FR identity against the Hu sequence should be near 100%
    if detection is correct, since the Hu framework IS that germline. We return
    this identity as a detection quality check.

    Returns:
        (gene, fr_identity) — gene name without allele, and FR identity score
        Both are None if detection fails.
    """
    try:
        numbered = number_sequence(hu_sequence, chain_type=chain_type)
        rankings = rank_germlines(
            query_fr=numbered["fr_residues"],
            chain_type=chain_type,
            top_n=20,           # return top 20 so we can find lab germline score
            min_identity=0.0,   # always return top hit regardless of score
        )
        if rankings:
            return rankings[0]["gene"], rankings[0]["fr_identity"], rankings
        return None, None, []
    except Exception:
        return None, None, []


# ── Helper: compare two sequences position by position ───────────────────────

def compare_sequences(seq_a: str, seq_b: str, chain_type: str) -> dict:
    """
    Compare two sequences using IMGT-numbered positions.

    Returns:
        match         — bool, True if sequences are identical
        identity      — float, fraction of matching positions
        n_differences — int, number of differing positions
        diffs         — list of per-position diff dicts:
                          position    : IMGT position number
                          pipeline_aa : residue in pipeline's grafted sequence
                          lab_aa      : residue in lab's Hu sequence
                          region      : FR1/CDR1/FR2/... region name
                          is_boundary : True if position is at a CDR/FR boundary
                          is_cdr      : True if position is inside a CDR
    """
    try:
        numbered_a = number_sequence(seq_a, chain_type=chain_type)
        numbered_b = number_sequence(seq_b, chain_type=chain_type)

        residues_a = {**numbered_a["fr_residues"],
                      **{k: v for k, v in numbered_a["cdr_residues"].items()
                         if isinstance(k, int)}}
        residues_b = {**numbered_b["fr_residues"],
                      **{k: v for k, v in numbered_b["cdr_residues"].items()
                         if isinstance(k, int)}}

        all_positions = set(residues_a.keys()) | set(residues_b.keys())
        comparable_positions = sorted(
            p for p in all_positions
            if residues_a.get(p) and residues_b.get(p)
        )

        if not comparable_positions:
            return {"match": False, "identity": 0.0, "n_differences": -1, "diffs": []}

        diffs = []
        matches = 0
        for p in comparable_positions:
            aa_a = residues_a[p]
            aa_b = residues_b[p]
            if aa_a == aa_b:
                matches += 1
            else:
                region = next(
                    (name for name, positions in IMGT_REGIONS.items() if p in positions),
                    "unknown"
                )
                diffs.append({
                    "position":    p,
                    "pipeline_aa": aa_a,
                    "lab_aa":      aa_b,
                    "region":      region,
                    "is_boundary": p in CDR_BOUNDARY_POSITIONS,
                    "is_cdr":      p in ALL_CDR_POSITIONS,
                })

        identity = matches / len(comparable_positions)
        return {
            "match":         len(diffs) == 0,
            "identity":      identity,
            "n_differences": len(diffs),
            "diffs":         diffs,
        }
    except Exception as e:
        return {"match": False, "identity": 0.0, "n_differences": -1,
                "diffs": [], "error": str(e)}


# ── Helper: get FR sequences for a named germline from the database ───────────

def get_germline_fr_by_region(germline_name: str, chain_type: str,
                              prefer_allele: str = None) -> dict:
    """
    Retrieve FR residues for a named germline, split by region.
    Returns {"FR1": {pos: aa}, "FR2": ..., "FR3": ..., "FR4": ...}

    Args:
        germline_name:  gene-level (IGHV3-23) or allele-level (IGHV3-23*01)
        chain_type:     'H', 'K', or 'L'
        prefer_allele:  if provided, use this specific allele instead of first match
                        Use this to ensure consistency with rank_germlines scoring.
    """
    from anarci.germlines import all_germlines
    human_germlines = all_germlines["V"][chain_type]["human"]

    # If specific allele requested, try exact match first
    matched_seq = None
    if prefer_allele and prefer_allele in human_germlines:
        matched_seq = human_germlines[prefer_allele]
    else:
        # Try exact match first, then gene-level prefix match
        # Use *01 allele preferentially for consistency
        gene = germline_name.split("*")[0]
        allele01 = f"{gene}*01"
        if germline_name in human_germlines:
            matched_seq = human_germlines[germline_name]
        elif allele01 in human_germlines:
            matched_seq = human_germlines[allele01]
        else:
            for name, aligned in human_germlines.items():
                if name.split("*")[0] == gene:
                    matched_seq = aligned
                    break

    if matched_seq is None:
        return {}

    # Parse aligned string into {pos: aa}
    all_residues = {
        pos_idx + 1: aa
        for pos_idx, aa in enumerate(matched_seq)
        if aa != "-"
    }

    # Split into FR regions
    fr_by_region = {}
    for region, positions in IMGT_REGIONS.items():
        if region.startswith("FR"):
            fr_by_region[region] = {
                p: all_residues[p] for p in positions if p in all_residues
            }
    return fr_by_region


# ── Helper: format side-by-side region comparison (same as verify.py) ─────────

def format_region_comparison(region: str, seq_a: dict, seq_b: dict,
                             label_a: str = "A", label_b: str = "B") -> list[str]:
    """
    Format a side-by-side comparison of two region sequences with diff markers.
    seq_a and seq_b are {imgt_pos: aa} dicts.
    Returns list of lines.
    """
    all_positions = sorted(set(seq_a.keys()) | set(seq_b.keys()))
    if not all_positions:
        return [f"    {region}: (empty in both)"]

    mouse_str = hu_str = diff_str = ""
    for pos in all_positions:
        aa_a = seq_a.get(pos, "-")
        aa_b = seq_b.get(pos, "-")
        mouse_str += aa_a
        hu_str += aa_b
        diff_str += " " if aa_a == aa_b else "^"

    lines = [f"    {region}:"]
    lines.append(f"      {label_a:<12}: {mouse_str}")
    lines.append(f"      {label_b:<12}: {hu_str}")

    if "^" in diff_str:
        lines.append(f"      {'Diff':<12}  {diff_str}")
        diffs = [
            f"pos{pos}({seq_a.get(pos, '-')}→{seq_b.get(pos, '-')})"
            for pos in all_positions
            if seq_a.get(pos, "-") != seq_b.get(pos, "-")
        ]
        lines.append(f"      {'Changes':<12}: {', '.join(diffs)}")
    else:
        lines.append(f"      {'Diff':<12}  (identical)")

    return lines


# ── Helper: print full debug for one clone ────────────────────────────────────

def print_clone_debug(r: CloneEvalResult) -> None:
    """
    Print three side-by-side comparisons for one clone:
      1. Mouse vs Hu — CDR and FR split verification
      2. Hu FRs vs top-1 detected germline FRs (ground truth quality check)
      3. Pipeline top-1 germline FRs vs lab ground truth germline FRs
    """
    if not r._mouse_vh_numbered and not r._mouse_vl_numbered:
        print("    (no debug data — run with --verbose)")
        return

    for chain_label, mouse_num, hu_num, rankings, chain_type, true_germ, pipe_rank1 in [
        ("VH", r._mouse_vh_numbered, r._hu_vh_numbered,
         r._vh_rankings, "H", r.vh_true_germline, r.vh_pipeline_rank1),
        ("VL", r._mouse_vl_numbered, r._hu_vl_numbered,
         r._vl_rankings, r._vl_chain_type, r.vl_true_germline, r.vl_pipeline_rank1),
    ]:
        if not mouse_num or not hu_num:
            continue

        print(
            f"\n  ── {chain_label} ──────────────────────────────────────────────────")

        # ── Pipeline top-N germline rankings ──────────────────────────────────
        print(f"\n  Pipeline top-{len(rankings)} germline rankings "
              f"(FR identity vs mouse, '← lab used' marks detected lab germline):")
        print(f"    {'Rank':<5} {'Germline':<22} {'FR%':>6} {'Match':>8} "
              f"{'Preferred':>10} {'':>12}")
        print(f"    {'-'*5} {'-'*22} {'-'*6} {'-'*8} {'-'*10} {'-'*12}")
        for rank, entry in enumerate(rankings, 1):
            preferred_mark = "✓" if entry["preferred"] else ""
            true_mark = "← lab used" if entry["gene"] == (
                true_germ or "") else ""
            print(
                f"    {rank:<5} {entry['germline']:<22} "
                f"{entry['fr_identity']:>5.1%} "
                f"{entry['matched']:>3}/{entry['comparable']:<4} "
                f"{preferred_mark:>10} "
                f"{true_mark:>12}"
            )
        if true_germ and not any(r["gene"] == true_germ for r in rankings):
            print(
                f"    ⚠ Lab germline '{true_germ}' not found in top {len(rankings)}")

        # ── Normalization report ───────────────────────────────────────────────
        from pipeline.step_b_germline_scoring import print_normalization_report
        print_normalization_report(rankings, chain_type)

        # ── Comparison 1: Mouse vs Hu — CDR and FR split verification ─────────
        print(
            f"\n  [1] Mouse vs Hu — CDR/FR split (confirms grafting and CDR preservation)")

        mouse_fr = mouse_num.get("fr_by_region",  {})
        mouse_cdr = mouse_num.get("cdr_by_region", {})
        hu_fr = hu_num.get("fr_by_region",  {})
        hu_cdr = hu_num.get("cdr_by_region", {})

        print(f"\n    CDR regions:")
        for region in ["CDR1", "CDR2", "CDR3"]:
            for line in format_region_comparison(
                region,
                mouse_cdr.get(region, {}),
                hu_cdr.get(region, {}),
                label_a="Mouse", label_b="Hu",
            ):
                print(line)

        print(f"\n    FR regions:")
        for region in ["FR1", "FR2", "FR3", "FR4"]:
            for line in format_region_comparison(
                region,
                mouse_fr.get(region, {}),
                hu_fr.get(region, {}),
                label_a="Mouse", label_b="Hu",
            ):
                print(line)

        # Legend for comparisons 2-4:
        #   a = detected ground truth germline FRs
        #   b = lab Hu FRs (extracted from grafted sequence)
        #   c = pipeline top-1 germline FRs

        pipe_germ_fr = get_germline_fr_by_region(
            pipe_rank1, chain_type) if pipe_rank1 else {}
        true_germ_fr = get_germline_fr_by_region(
            true_germ,  chain_type) if true_germ else {}

        # ── Comparison 2: a vs b — detected germline FRs vs lab Hu FRs ──────────
        # If detection is correct, these should be near identical.
        # Large differences here mean germline detection is unreliable.
        print(f"\n  [2] a vs b — Detected germline FRs vs Lab Hu FRs")
        print(f"      a (detected germline) : {true_germ or 'N/A'}")
        print(f"      b (lab Hu sequence)   : extracted FR residues")
        if true_germ and true_germ_fr:
            for region in ["FR1", "FR2", "FR3", "FR4"]:
                for line in format_region_comparison(
                    region,
                    true_germ_fr.get(region, {}),
                    hu_fr.get(region, {}),
                    label_a=f"a:{true_germ[:10]}", label_b="b:Hu FR",
                ):
                    print(line)
        else:
            print(f"      (no detected germline — cannot compare)")

        # ── Comparison 3: b vs c — lab Hu FRs vs pipeline top-1 germline FRs ────
        # Shows how far the pipeline's chosen germline is from what the lab used.
        # Differences here directly explain grafting identity gaps.
        print(f"\n  [3] b vs c — Lab Hu FRs vs Pipeline top-1 germline FRs")
        print(f"      b (lab Hu sequence)   : extracted FR residues")
        print(f"      c (pipeline top-1)    : {pipe_rank1 or 'N/A'}")
        if pipe_rank1 and pipe_germ_fr:
            for region in ["FR1", "FR2", "FR3", "FR4"]:
                for line in format_region_comparison(
                    region,
                    hu_fr.get(region, {}),
                    pipe_germ_fr.get(region, {}),
                    label_a="b:Hu FR", label_b=f"c:{pipe_rank1[:10]}",
                ):
                    print(line)
        else:
            print(f"      (no pipeline germline — cannot compare)")

        # ── Comparison 4: a vs c — detected germline FRs vs pipeline top-1 FRs ──
        # Shows whether the two germlines agree even if they have different names.
        # If a==c, the pipeline chose the right germline but detection assigned
        # a different allele or very similar gene.
        print(
            f"\n  [4] a vs c — Detected germline FRs vs Pipeline top-1 germline FRs")
        print(f"      a (detected germline) : {true_germ or 'N/A'}")
        print(f"      c (pipeline top-1)    : {pipe_rank1 or 'N/A'}")
        if pipe_rank1 and true_germ:
            if pipe_rank1 == true_germ:
                print(f"      → Same germline — FRs are identical by definition ✓")
            elif pipe_germ_fr and true_germ_fr:
                for region in ["FR1", "FR2", "FR3", "FR4"]:
                    for line in format_region_comparison(
                        region,
                        true_germ_fr.get(region, {}),
                        pipe_germ_fr.get(region, {}),
                        label_a=f"a:{true_germ[:10]}", label_b=f"c:{pipe_rank1[:10]}",
                    ):
                        print(line)
            else:
                print(f"      (one or both germlines not found in database)")
        else:
            print(f"      (missing germline info — cannot compare)")

        # ── Comparison 5: mouse FRs vs pipeline top-1 germline FRs ──────────────
        # Shows the structural distance between the mouse sequence and the pipeline's
        # chosen human germline. Large differences = challenging humanization case.
        # Small differences = mouse sequence is already close to human germline.
        print(f"\n  [5] Mouse FRs vs Pipeline top-1 germline FRs")
        print(f"      Mouse sequence  : input")
        print(f"      c (pipeline top-1) : {pipe_rank1 or 'N/A'}")
        if pipe_rank1 and pipe_germ_fr:
            mouse_fr_flat = mouse_num.get("fr_by_region", {})
            for region in ["FR1", "FR2", "FR3", "FR4"]:
                for line in format_region_comparison(
                    region,
                    mouse_fr_flat.get(region, {}),
                    pipe_germ_fr.get(region, {}),
                    label_a="Mouse", label_b=f"c:{pipe_rank1[:10]}",
                ):
                    print(line)
        else:
            print(f"      (no pipeline germline — cannot compare)")

        # ── Comparison 6: mouse FRs vs detected lab germline FRs ─────────────────
        # Shows how close the mouse sequence is to the germline the lab actually used.
        # If [5] and [6] show similar difference counts, the pipeline's germline choice
        # is as good as the lab's in terms of structural distance from the mouse.
        # If [6] shows fewer differences, the lab's germline is a better fit.
        print(f"\n  [6] Mouse FRs vs Detected lab germline FRs")
        print(f"      Mouse sequence       : input")
        print(f"      a (detected germline): {true_germ or 'N/A'}")
        if true_germ and true_germ_fr:
            mouse_fr_flat = mouse_num.get("fr_by_region", {})
            for region in ["FR1", "FR2", "FR3", "FR4"]:
                for line in format_region_comparison(
                    region,
                    mouse_fr_flat.get(region, {}),
                    true_germ_fr.get(region, {}),
                    label_a="Mouse", label_b=f"a:{true_germ[:10]}",
                ):
                    print(line)
        else:
            print(f"      (no detected germline — cannot compare)")


# ── Core evaluation: Phase 1 — FR analysis only (no abnumber needed) ─────────

def evaluate_clone_fr(clone_id: str, row: dict, top_n: int = 10) -> CloneEvalResult:
    """
    Phase 1 evaluation — FR analysis only.
    Uses ANARCI database throughout. No abnumber dependency.

    Covers:
      - Germline detection from lab Hu sequence (FR-only scoring)
      - Pipeline germline ranking from mouse sequence
      - Recall@K — does pipeline rank lab's germline in top N?
      - FR-level comparisons: a vs b, b vs c, a vs c
    """
    result = CloneEvalResult(clone_id=clone_id)

    mouse_vh = row.get("mouse_vh",        "").strip()
    mouse_vl = row.get("mouse_vl",        "").strip()
    hu_vh = row.get("hu_vh",           "").strip()
    hu_vl = row.get("hu_vl",           "").strip()
    lab_vh_germline = row.get("vh_germline", "").strip()
    lab_vl_germline = row.get("vl_germline", "").strip()

    if not mouse_vh or not mouse_vl:
        result.error = "Missing mouse sequences"
        return result
    if not hu_vh or not hu_vl:
        result.error = "Missing Hu sequences"
        return result

    try:
        # Step A: number all sequences
        vh_numbered = number_sequence(mouse_vh, chain_type="H")
        vl_numbered = number_sequence(mouse_vl, chain_type=None)
        vl_chain_type = vl_numbered["chain_type"]
        hu_vh_numbered = number_sequence(hu_vh, chain_type="H")
        hu_vl_numbered = number_sequence(hu_vl, chain_type=vl_chain_type)

        # Step B: rank germlines from mouse FRs
        vh_rankings = rank_germlines(
            vh_numbered["fr_residues"], "H",           top_n=top_n)
        vl_rankings = rank_germlines(
            vl_numbered["fr_residues"], vl_chain_type, top_n=top_n)

        # Store intermediates for verbose debug output
        result._mouse_vh_numbered = vh_numbered
        result._mouse_vl_numbered = vl_numbered
        result._hu_vh_numbered = hu_vh_numbered
        result._hu_vl_numbered = hu_vl_numbered
        result._vh_rankings = vh_rankings
        result._vl_rankings = vl_rankings
        result._vl_chain_type = vl_chain_type

        # Always detect germline from Hu sequence — stored separately from lab-provided
        # Always detect germline from Hu sequence — stored separately from lab-provided
        result.vh_detected_germline, result.vh_true_germline_fr_identity, result._hu_vh_det_rankings = \
            detect_germline(hu_vh, "H")
        result.vl_detected_germline, result.vl_true_germline_fr_identity, result._hu_vl_det_rankings = \
            detect_germline(hu_vl, vl_chain_type)
        # Ground truth germline — prefer lab-provided over detected
        # Lab names may include functional suffixes (e.g. *02F) — strip to gene level
        if lab_vh_germline:
            result.vh_true_germline = parse_lab_germline(lab_vh_germline)
            print(
                f"  VH: lab-provided germline: {lab_vh_germline} parsed to {result.vh_true_germline} | detected: {result.vh_detected_germline}")
        else:
            result.vh_true_germline = result.vh_detected_germline

        if lab_vl_germline:
            result.vl_true_germline = parse_lab_germline(lab_vl_germline)
            print(
                f"  VL: lab-provided germline: {lab_vl_germline} parsed to {result.vl_true_germline} | detected: {result.vl_detected_germline}")
        else:
            result.vl_true_germline = result.vl_detected_germline

        if vh_rankings:
            result.vh_pipeline_rank1 = vh_rankings[0]["gene"]
        if vl_rankings:
            result.vl_pipeline_rank1 = vl_rankings[0]["gene"]

        # Recall@K — does pipeline rank lab's germline?
        vh_genes = [r["gene"] for r in vh_rankings]
        vl_genes = [r["gene"] for r in vl_rankings]

        for attr, true_g, genes in [
            ("vh", result.vh_true_germline, vh_genes),
            ("vl", result.vl_true_germline, vl_genes),
        ]:
            if true_g in genes:
                rank = genes.index(true_g) + 1
                setattr(result, f"{attr}_true_germline_rank", rank)
                setattr(result, f"{attr}_recall_at_1", rank <= 1)
                setattr(result, f"{attr}_recall_at_3", rank <= 3)
                setattr(result, f"{attr}_recall_at_5", rank <= 5)
            else:
                setattr(result, f"{attr}_true_germline_rank", None)
                setattr(result, f"{attr}_recall_at_1", False)
                setattr(result, f"{attr}_recall_at_3", False)
                setattr(result, f"{attr}_recall_at_5", False)

    except Exception as e:
        result.error = str(e)

    return result


# ── Core evaluation: Phase 2 — Grafting accuracy (requires abnumber) ──────────

def evaluate_clone_grafting(result: CloneEvalResult, row: dict) -> CloneEvalResult:
    """
    Phase 2 evaluation — grafting accuracy.
    Requires abnumber. Run after evaluate_clone_fr().

    Grafts mouse CDRs onto pipeline's top-1 germline (using abnumber),
    then compares the resulting full sequence against the lab's Hu sequence.
    """
    if result.error:
        return result  # skip if Phase 1 failed

    mouse_vh = row.get("mouse_vh", "").strip()
    mouse_vl = row.get("mouse_vl", "").strip()
    hu_vh = row.get("hu_vh",    "").strip()
    hu_vl = row.get("hu_vl",    "").strip()

    vh_rankings = getattr(result, "_vh_rankings",   [])
    vl_rankings = getattr(result, "_vl_rankings",   [])
    vl_chain_type = getattr(result, "_vl_chain_type", "K")

    try:
        if vh_rankings:
            vh_germ_normalized = normalize_germline_name(
                vh_rankings[0]["germline"])
            mouse_vh_chain = Chain(
                mouse_vh, scheme="imgt", cdr_definition="imgt")
            grafted_vh = mouse_vh_chain.graft_cdrs_onto_human_germline(
                v_gene=vh_germ_normalized, backmutate_vernier=False)
            result.vh_grafted_seq = grafted_vh.seq
            result.vh_hu_seq = hu_vh
            if vh_germ_normalized != vh_rankings[0]["germline"]:
                print(f"  ⚠ VH normalized: {vh_rankings[0]['germline']} "
                      f"→ {vh_germ_normalized}", end=" ")

            cmp = compare_sequences(result.vh_grafted_seq, hu_vh, "H")
            result.vh_seq_match = cmp["match"]
            result.vh_seq_identity = cmp["identity"]
            result.vh_n_differences = cmp["n_differences"]
            result.vh_diffs = cmp["diffs"]

        if vl_rankings:
            vl_germ_normalized = normalize_germline_name(
                vl_rankings[0]["germline"])
            mouse_vl_chain = Chain(
                mouse_vl, scheme="imgt", cdr_definition="imgt")
            grafted_vl = mouse_vl_chain.graft_cdrs_onto_human_germline(
                v_gene=vl_germ_normalized, backmutate_vernier=False)
            result.vl_grafted_seq = grafted_vl.seq
            result.vl_hu_seq = hu_vl
            if vl_germ_normalized != vl_rankings[0]["germline"]:
                print(f"  ⚠ VL normalized: {vl_rankings[0]['germline']} "
                      f"→ {vl_germ_normalized}", end=" ")

            cmp = compare_sequences(
                result.vl_grafted_seq, hu_vl, vl_chain_type)
            result.vl_seq_match = cmp["match"]
            result.vl_seq_identity = cmp["identity"]
            result.vl_n_differences = cmp["n_differences"]
            result.vl_diffs = cmp["diffs"]

    except Exception as e:
        result.error = f"Grafting failed: {e}"

    return result


# ── Convenience wrapper — runs both phases ────────────────────────────────────

def evaluate_clone(clone_id: str, row: dict, top_n: int = 10,
                   run_grafting: bool = True) -> CloneEvalResult:
    """
    Run FR analysis (Phase 1) and optionally grafting accuracy (Phase 2).

    Args:
        run_grafting: if False, skips Phase 2 (no abnumber needed).
                      Use --mode fr in CLI to set this.
    """
    result = evaluate_clone_fr(clone_id, row, top_n=top_n)
    if run_grafting and not result.error:
        result = evaluate_clone_grafting(result, row)
    return result


# ── Output: per-clone position-level breakdown ────────────────────────────────

def print_position_diffs(result: CloneEvalResult) -> None:
    """Print a detailed breakdown of where differences fall for one clone."""
    for chain, diffs in [("VH", result.vh_diffs), ("VL", result.vl_diffs)]:
        if not diffs:
            continue
        print(f"    {chain} differences ({len(diffs)} positions):")
        print(
            f"      {'Pos':>5}  {'Region':<8}  {'Pipeline':>8}  {'Lab':>8}  {'Note'}")
        print(f"      {'-'*5}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*20}")
        for d in diffs:
            note = ""
            if d["is_boundary"]:
                note = "⚠ CDR/FR boundary"
            elif d["is_cdr"]:
                note = "CDR region"
            print(f"      {d['position']:>5}  {d['region']:<8}  "
                  f"{d['pipeline_aa']:>8}  {d['lab_aa']:>8}  {note}")


# ── Output: aggregate summary ─────────────────────────────────────────────────


# ── Germline match table (diagnose mode) ──────────────────────────────────────

def print_germline_match_table(results: list[CloneEvalResult]) -> None:
    """
    Print a simple germline match table:
      clone | VH detected | VH lab | match | VL detected | VL lab | match
    For mismatches, print side-by-side FR comparison of detected vs lab germline.
    """
    print(f"\nGERMLINE DETECTION — Detected from Hu sequence vs Lab actual")
    print("-" * 75)
    print(f"{'Clone':<10} {'VH Detected':>14} {'VH Lab':>14} {'Match':>6}  "
          f"{'VL Detected':>14} {'VL Lab':>14} {'Match':>6}")
    print("-" * 75)

    failures = []
    for r in results:
        if r.error:
            print(f"{r.clone_id:<10} ERROR: {r.error}")
            continue

        vh_detected = r.vh_detected_germline or "N/A"
        vl_detected = r.vl_detected_germline or "N/A"
        vh_lab = r.vh_true_germline or "N/A"
        vl_lab = r.vl_true_germline or "N/A"
        vh_match = vh_detected == vh_lab
        vl_match = vl_detected == vl_lab

        # Check tie using Hu-based detection rankings
        vh_hu_det = getattr(r, "_hu_vh_det_rankings", [])
        vl_hu_det = getattr(r, "_hu_vl_det_rankings", [])

        def get_tie_info(hu_det_rankings, lab_germ, match):
            """Return (show_tie, n_tied) based on rules:
              - match=True:  show ~x (ties at top regardless)
              - match=False: show ~x only if lab and detected share same score
            """
            if not hu_det_rankings:
                return False, 0
            det_score = hu_det_rankings[0]["fr_identity"]
            lab_entry = next(
                (e for e in hu_det_rankings if e["gene"] == lab_germ), None
            )
            all_at_top = [e for e in hu_det_rankings
                          if e["fr_identity"] == det_score]
            n = len(set(e["gene"] for e in all_at_top))
            if match:
                return n > 1, n
            else:
                # Only show tie if lab germline shares the detected score
                lab_tied = (lab_entry is not None and
                            lab_entry["fr_identity"] == det_score)
                return lab_tied, n if lab_tied else 0

        vh_show_tie, vh_n_tied = get_tie_info(vh_hu_det, vh_lab, vh_match)
        vl_show_tie, vl_n_tied = get_tie_info(vl_hu_det, vl_lab, vl_match)

        vh_mark = ("✓" if vh_match else "✗") + \
            (f" ~{vh_n_tied}" if vh_show_tie else "")
        vl_mark = ("✓" if vl_match else "✗") + \
            (f" ~{vl_n_tied}" if vl_show_tie else "")

        print(f"{r.clone_id:<10} {vh_detected:>14} {vh_lab:>14} {vh_mark:>8}  "
              f"{vl_detected:>14} {vl_lab:>14} {vl_mark:>8}")

        if not vh_match or not vl_match:
            failures.append((r, vh_match, vl_match))

    if not failures:
        print("\n✓ All clones matched correctly")
        return

    # For each failure, print FR comparison of detected vs lab germline
    print(f"\n{'='*65}")
    print(f"MISMATCH DETAIL — FR comparison: detected germline vs lab actual")
    print(f"{'='*65}")

    for r, vh_match, vl_match in failures:
        debug_data = getattr(r, "_vh_rankings", None)
        if debug_data is None:
            print(f"\n  Clone {r.clone_id}: no debug data available")
            continue

        print(f"\n  Clone: {r.clone_id}")
        vl_chain_type = getattr(r, "_vl_chain_type", "K")

        for chain_label, match, detected, lab_germ, chain_type in [
            ("VH", vh_match, r.vh_detected_germline, r.vh_true_germline, "H"),
            ("VL", vl_match, r.vl_detected_germline,
             r.vl_true_germline, vl_chain_type),
        ]:
            if match:
                continue  # only expand failures

            print(
                f"\n  ── {chain_label} MISMATCH ──────────────────────────────")
            print(f"     Detected         : {detected or 'N/A'}")
            print(f"     Lab actual        : {lab_germ or 'N/A'}")

            # Show tie warning ONLY when detected and lab germline share the same score
            # Use Hu-based detection rankings (not mouse-based pipeline rankings)
            hu_det_rankings = (
                r._hu_vh_det_rankings if chain_label == "VH"
                else r._hu_vl_det_rankings
            )
            det_score = hu_det_rankings[0]["fr_identity"] if hu_det_rankings else None
            lab_entry = next(
                (e for e in hu_det_rankings if e["gene"] == lab_germ), None
            )
            lab_score = lab_entry["fr_identity"] if lab_entry else None

            if det_score is not None and lab_score is not None and det_score == lab_score:
                # Detected and lab germline tied — selection was arbitrary
                all_tied = list(dict.fromkeys(
                    e["germline"] for e in hu_det_rankings
                    if e["fr_identity"] == det_score
                ))
                n = len(all_tied)
                print(
                    f"     ⚠ TIE: detected and lab germline share same score ({det_score:.1%})")
                print(f"       Selection was arbitrary among {n} germlines:")
                marked = [
                    f"{g} ← lab used" if g.split("*")[0] == lab_germ else g
                    for g in all_tied
                ]
                print(
                    f"       {', '.join(marked) if marked else '(none)'}")

            # Use Hu-based detection rankings to find correct allele
            hu_det_rankings = (
                r._hu_vh_det_rankings if chain_label == "VH"
                else r._hu_vl_det_rankings
            )
            det_allele = next(
                (e["germline"]
                 for e in hu_det_rankings if e["gene"] == detected), None
            )
            det_fr = get_germline_fr_by_region(detected, chain_type,
                                               prefer_allele=det_allele) if detected else {}
            lab_fr = get_germline_fr_by_region(
                lab_germ,  chain_type) if lab_germ else {}

            # If lab germline not in database, skip Options A and B
            # but still show Option C (Hu FRs vs detected germline)
            lab_found = bool(lab_fr)
            det_found = bool(det_fr)

            if not det_found:
                print(
                    f"     (detected germline {detected} not found in database)")
                continue

            if not lab_found:
                print(
                    f"     (lab germline {lab_germ} not found in database — showing Hu vs detected only)")

            if lab_found:
                print(f"\n     FR comparison (detected germline vs lab actual):")
                for region in ["FR1", "FR2", "FR3", "FR4"]:
                    for line in format_region_comparison(
                        region,
                        det_fr.get(region, {}),
                        lab_fr.get(region, {}),
                        label_a=f"{(detected or '')[:12]}",
                        label_b=f"{(lab_germ or '')[:12]}",
                    ):
                        print(f"  {line}")

            # Option B: Hu FR residues (from input file) vs lab stated germline
            # Shows how the actual grafted sequence compares to the lab's stated germline
            hu_numbered = (
                r._hu_vh_numbered if chain_label == "VH"
                else r._hu_vl_numbered
            )
            hu_fr_regions = hu_numbered.get("fr_by_region", {})

            if hu_fr_regions and lab_found:
                print(
                    f"\n     Hu sequence FRs (from input) vs Lab stated germline ({lab_germ}):")
                print(f"     (FR4 excluded — comes from J gene, not in V gene database)")
                for region in ["FR1", "FR2", "FR3"]:
                    for line in format_region_comparison(
                        region,
                        hu_fr_regions.get(region, {}),
                        lab_fr.get(region, {}),
                        label_a="Hu FR",
                        label_b=f"{(lab_germ or '')[:12]}",
                    ):
                        print(f"  {line}")

            # Option C: Hu FR residues (from input file) vs detected germline
            # Shows how the actual grafted sequence compares to what the pipeline detected
            if hu_fr_regions and det_fr:
                print(
                    f"\n     Hu sequence FRs (from input) vs Detected germline ({detected}):")
                print(f"     (FR4 excluded — comes from J gene, not in V gene database)")
                for region in ["FR1", "FR2", "FR3"]:
                    for line in format_region_comparison(
                        region,
                        hu_fr_regions.get(region, {}),
                        det_fr.get(region, {}),
                        label_a="Hu FR",
                        label_b=f"{(detected or '')[:12]}",
                    ):
                        print(f"  {line}")
            else:
                print(f"\n     (cannot compare Hu vs detected — missing data)")


def print_summary(results: list[CloneEvalResult], verbose: bool = False) -> None:
    valid = [r for r in results if r.error is None]
    n = len(valid)

    if n == 0:
        print("No valid results to summarize.")
        return

    print(f"\n{'='*65}")
    print(f"EVALUATION SUMMARY  ({n} clones evaluated, "
          f"{len(results)-n} skipped due to errors)")
    print(f"{'='*65}")

    # ── Germline selection ─────────────────────────────────────────────────
    print("\nGERMLINE SELECTION")
    print("-" * 40)
    for chain, attr in [("VH", "vh"), ("VL", "vl")]:
        r1 = sum(getattr(r, f"{attr}_recall_at_1") or False for r in valid)
        r3 = sum(getattr(r, f"{attr}_recall_at_3") or False for r in valid)
        r5 = sum(getattr(r, f"{attr}_recall_at_5") or False for r in valid)
        ranks = [getattr(r, f"{attr}_true_germline_rank")
                 for r in valid if getattr(r, f"{attr}_true_germline_rank")]
        avg_rank = sum(ranks) / len(ranks) if ranks else None
        not_found = sum(1 for r in valid
                        if getattr(r, f"{attr}_true_germline_rank") is None)

        print(f"  {chain}:")
        print(f"    Recall@1:  {r1}/{n} ({r1/n:.0%})")
        print(f"    Recall@3:  {r3}/{n} ({r3/n:.0%})")
        print(f"    Recall@5:  {r5}/{n} ({r5/n:.0%})")
        if avg_rank:
            print(f"    Avg rank of true germline: {avg_rank:.1f}")
        if not_found:
            print(
                f"    Not in top {n}: {not_found} clone(s) — germline outside ranked list")

    # ── Grafting accuracy ──────────────────────────────────────────────────
    grafting_ran = any(r.vh_grafted_seq is not None for r in valid)
    print("\nGRAFTING ACCURACY")
    print("-" * 40)
    if not grafting_ran:
        print("  (skipped — run with --mode full to include grafting accuracy)")
    else:
        for chain, attr in [("VH", "vh"), ("VL", "vl")]:
            exact = sum(getattr(r, f"{attr}_seq_match")
                        or False for r in valid)
            ids = [getattr(r, f"{attr}_seq_identity")
                   for r in valid if getattr(r, f"{attr}_seq_identity") is not None]
            diffs = [getattr(r, f"{attr}_n_differences")
                     for r in valid
                     if getattr(r, f"{attr}_n_differences") is not None
                     and getattr(r, f"{attr}_n_differences") >= 0]
            avg_id = sum(ids) / len(ids) if ids else None
            avg_diff = sum(diffs) / len(diffs) if diffs else None

            print(f"  {chain}:")
            print(f"    Exact match:      {exact}/{n} ({exact/n:.0%})")
            if avg_id:
                print(f"    Avg seq identity: {avg_id:.1%}")
            if avg_diff is not None:
                print(f"    Avg differences:  {avg_diff:.1f} positions")

    # ── Position frequency analysis ────────────────────────────────────────
    print("\nPOSITION FREQUENCY ANALYSIS")
    print("  (positions that differ across multiple clones — systematic issues)")
    print("-" * 55)

    for chain, attr in [("VH", "vh"), ("VL", "vl")]:
        all_diffs = []
        for r in valid:
            all_diffs.extend(getattr(r, f"{attr}_diffs") or [])

        if not all_diffs:
            print(f"  {chain}: no differences found")
            continue

        pos_counter = Counter(d["position"] for d in all_diffs)
        # Only show positions that appear in more than one clone
        systematic = [(pos, count) for pos, count in pos_counter.most_common()
                      if count > 1]

        if not systematic:
            print(f"  {chain}: no systematic (multi-clone) differences found")
            continue

        print(f"  {chain} — positions differing in multiple clones:")
        print(f"    {'Pos':>5}  {'Clones':>7}  {'Region':<8}  {'Boundary?':>10}  "
              f"{'Pipeline→Lab changes'}")
        print(f"    {'-'*5}  {'-'*7}  {'-'*8}  {'-'*10}  {'-'*25}")

        for pos, count in systematic:
            # Get region and boundary info from first occurrence
            sample = next(d for d in all_diffs if d["position"] == pos)
            # Summarise the amino acid changes seen
            changes = Counter(
                f"{d['pipeline_aa']}→{d['lab_aa']}"
                for d in all_diffs if d["position"] == pos
            )
            changes_str = ", ".join(
                f"{ch}({cnt})" for ch, cnt in changes.most_common())
            boundary_mark = "⚠ YES" if sample["is_boundary"] else "no"
            print(f"    {pos:>5}  {count:>7}  {sample['region']:<8}  "
                  f"{boundary_mark:>10}  {changes_str}")

        boundary_count = sum(1 for pos, _ in systematic
                             if any(d["is_boundary"] for d in all_diffs
                                    if d["position"] == pos))
        if boundary_count > 0:
            print(
                f"\n  ⚠ {boundary_count} systematic difference(s) at CDR/FR boundaries")
            print(f"    → Strong signal of CDR definition mismatch (IMGT vs Kabat)")
            print(f"    → Try re-running with cdr_definition='kabat' in grafting step")

    # ── Per-clone detail ───────────────────────────────────────────────────
    print(f"\nPER-CLONE DETAIL")
    print("-" * 65)
    print(f"{'Clone':<10} {'VH Germline':>12} {'Det%':>5} {'Rank':>5} {'Graft%':>7} "
          f"{'VL Germline':>12} {'Det%':>5} {'Rank':>5} {'Graft%':>7}")
    print("-" * 80)
    for r in results:
        if r.error:
            print(f"{r.clone_id:<10} ERROR: {r.error}")
            continue
        vh_rank = str(
            r.vh_true_germline_rank) if r.vh_true_germline_rank else "N/F"
        vl_rank = str(
            r.vl_true_germline_rank) if r.vl_true_germline_rank else "N/F"
        vh_graft = f"{r.vh_seq_identity:.1%}" if r.vh_seq_identity is not None else "N/A"
        vl_graft = f"{r.vl_seq_identity:.1%}" if r.vl_seq_identity is not None else "N/A"
        vh_det = f"{r.vh_true_germline_fr_identity:.0%}" if r.vh_true_germline_fr_identity is not None else "N/A"
        vl_det = f"{r.vl_true_germline_fr_identity:.0%}" if r.vl_true_germline_fr_identity is not None else "N/A"
        print(f"{r.clone_id:<10} "
              f"{(r.vh_true_germline or 'N/A'):>12} {vh_det:>5} {vh_rank:>5} {vh_graft:>7} "
              f"{(r.vl_true_germline or 'N/A'):>12} {vl_det:>5} {vl_rank:>5} {vl_graft:>7}")

    # ── Per-clone deep debug (verbose mode) ──────────────────────────────
    if verbose:
        print("\nPER-CLONE DEBUG — CDR/FR split, germline detection, germline comparison")
        print("=" * 65)
        for r in results:
            if r.error:
                continue
            print(f"\nClone: {r.clone_id}")
            print_clone_debug(r)
            if r.vh_diffs or r.vl_diffs:
                print(
                    "\n  [4] Grafted sequence differences (pipeline vs lab Hu):")
                print_position_diffs(r)


# ── Export ────────────────────────────────────────────────────────────────────

def export_results(results: list[CloneEvalResult], output_path: str) -> None:
    """Export per-clone summary results to CSV (excludes diff details)."""
    fields = [f for f in CloneEvalResult.__dataclass_fields__
              if not f.endswith("_diffs") and not f.endswith("_seq")]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            row = {field: getattr(r, field) for field in fields}
            for key in ["vh_seq_identity", "vl_seq_identity"]:
                if key in row and row[key] is not None:
                    row[key] = f"{row[key]:.4f}"
            writer.writerow(row)
    print(f"\nResults exported to: {output_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate pipeline against lab ground truth")
    parser.add_argument("--csv",     required=True,
                        help="Path to benchmark CSV")
    parser.add_argument("--top-n",   type=int, default=10,
                        help="Germline candidates to rank per chain (default: 10)")
    parser.add_argument("--output",  default="outputs/evaluation_results.csv",
                        help="Output CSV path")
    parser.add_argument("--verbose", action="store_true",
                        help="Print per-clone position-level diff breakdown")
    parser.add_argument("--mode",    choices=["fr", "diagnose", "full"], default="full",
                        help=(
                            "fr       — FR analysis + summary table only. "
                            "diagnose — FR analysis + germline match table, "
                            "           expands mismatches with FR comparison. "
                            "full     — complete output including verbose per-clone detail (default)."
    ))
    args = parser.parse_args()

    run_grafting = args.mode == "full"
    mode_label = {
        "fr":       "FR analysis only",
        "diagnose": "FR analysis + germline match table (mismatches expanded)",
        "full":     "FR analysis + grafting + full verbose detail",
    }[args.mode]
    print(f"Loading: {args.csv}")
    print(f"Mode: {args.mode} — {mode_label}")
    print(f"Top-N germlines: {args.top_n}\n")

    with open(args.csv, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, dialect="excel")
        rows = [row for row in reader if row.get("mouse_vh", "").strip()]

    print(f"Found {len(rows)} clones with sequences\n")

    results = []
    for row in rows:
        clone_id = row.get("clone", "unknown").strip()
        print(f"Evaluating {clone_id}...", end=" ", flush=True)
        result = evaluate_clone(
            clone_id, row, top_n=args.top_n, run_grafting=run_grafting)
        if result.error:
            print(f"ERROR: {result.error}")
        else:
            vh_id = f"{result.vh_seq_identity:.1%}" if result.vh_seq_identity is not None else "N/A"
            vl_id = f"{result.vl_seq_identity:.1%}" if result.vl_seq_identity is not None else "N/A"
            print(f"VH rank={result.vh_true_germline_rank}, VH id={vh_id}, "
                  f"VL rank={result.vl_true_germline_rank}, VL id={vl_id}")
        results.append(result)

    if args.mode == "diagnose":
        print_summary(results, verbose=False)
        print_germline_match_table(results)
    else:
        print_summary(results, verbose=(args.mode == "full"))
    export_results(results, args.output)


if __name__ == "__main__":
    main()
