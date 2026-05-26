"""
generate_sequences.py

Generates all 9 evaluation sequences for a given antibody clone.

Sequences produced:
  0.  sapiens_on_mouse             — mouse → Sapiens → CDR restored (pure model output)
  0r. sapiens_on_mouse_raw         — mouse → Sapiens raw (before CDR restoration)
  1.  pipeline_grafted             — mouse CDRs + pipeline top-1 germline FRs
  2.  pipeline_humanized           — seq 1 + Sapiens → CDR restored
  2r. pipeline_humanized_raw       — seq 1 + Sapiens raw
  3.  lab_grafted                  — lab's Hu sequence (ground truth, from CSV)
  4.  detected_germline_grafted    — mouse CDRs + detected lab germline FRs
  5.  lab_final                    — lab's final humanized sequence (ground truth, from CSV)
  6.  detected_humanized           — seq 4 + Sapiens → CDR restored
  6r. detected_humanized_raw       — seq 4 + Sapiens raw
  7.  detected_direct_backmut      — seq 4 + back-mutations where seq3 != seq5
  8.  lab_stated_germline_grafted  — mouse CDRs + lab's stated germline FRs (from database)
  9.  lab_stated_germline_humanized — seq 8 + Sapiens → CDR restored
  9r. lab_stated_germline_humanized_raw — seq 8 + Sapiens raw

Usage (from project root):
    python3 pipeline/generate_sequences.py \\
        --csv data/benchmarks/humanization_benchmark.csv \\
        --output outputs/all_sequences.csv
"""

# isort: skip_file
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E402
from abnumber import Chain
from pipeline.step_b_germline_scoring import rank_germlines, normalize_germline_name, print_normalization_report
from pipeline.step_a_numbering import number_sequence, IMGT_REGIONS
from evaluation.evaluate import detect_germline, get_germline_fr_by_region
from typing import Optional
from dataclasses import dataclass, field
import argparse
import csv
import json
import re


# ── Germline sequence retrieval ───────────────────────────────────────────────

def get_germline_seq_by_region(germline_name: str, chain_type: str) -> dict:
    """
    Retrieve the germline sequence split by all regions (FR + CDR).
    Wraps evaluate.py's get_germline_fr_by_region (which returns FR only)
    and extends it with CDR regions using the same ANARCI lookup.

    Returns {region: {str(imgt_pos): aa}} for FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4.
    Stored in all_sequences.csv so the report can render side-by-side
    sequence vs germline comparisons without needing database access at render time.
    """
    try:
        from anarci.germlines import all_germlines
        human_germlines = all_germlines["V"][chain_type]["human"]

        gene = germline_name.split("*")[0]
        allele01 = f"{gene}*01"
        matched_seq = None
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

        # Parse aligned IMGT string into {pos: aa} — same as evaluate.py
        all_residues = {
            pos_idx + 1: aa
            for pos_idx, aa in enumerate(matched_seq)
            if aa != "-"
        }

        # Split into all regions using IMGT_REGIONS (same dict used in evaluate.py)
        result = {}
        for region, positions in IMGT_REGIONS.items():
            region_res = {str(p): all_residues[p]
                          for p in positions if p in all_residues}
            if region_res:
                result[region] = region_res
        return result

    except Exception as e:
        print(
            f"    Warning: could not retrieve germline sequence for {germline_name}: {e}")
        return {}


# ── Lab germline name parsing (same as evaluate.py) ──────────────────────────

def parse_lab_germline(name: str) -> Optional[str]:
    """
    Normalize a lab-provided germline name to gene level.
    Strips allele and functional suffixes: IGHV3-23*01F → IGHV3-23
    Fixes transpositions: IGVK4-1 → IGKV4-1
    """
    if not name or not name.strip():
        return None
    name = name.strip()
    name = re.sub(r'^IGVK', 'IGKV', name)
    name = re.sub(r'^IGVL', 'IGLV', name)
    name = re.sub(r'^IGVH', 'IGHV', name)
    name = re.sub(r'[*]\d+[A-Za-z]*$', '', name)
    return name.strip() or None


# ── Sequence container ────────────────────────────────────────────────────────

@dataclass
class CloneSequences:
    """All 9 evaluation sequences for one antibody clone."""
    clone_id: str

    # VH sequences
    vh_0_sapiens_on_mouse:              Optional[str] = None
    vh_0_sapiens_on_mouse_raw:          Optional[str] = None
    vh_1_pipeline_grafted:              Optional[str] = None
    vh_2_pipeline_humanized:            Optional[str] = None
    vh_3_lab_grafted:                   Optional[str] = None
    vh_4_detected_grafted:              Optional[str] = None
    vh_5_lab_final:                     Optional[str] = None
    vh_6_detected_humanized:            Optional[str] = None
    vh_7_detected_direct_backmut:       Optional[str] = None
    vh_8_stated_germline_grafted:       Optional[str] = None
    vh_9_stated_germline_humanized:     Optional[str] = None

    # VH Sapiens raw output (before CDR restoration)
    vh_2_pipeline_humanized_raw:        Optional[str] = None
    vh_6_detected_humanized_raw:        Optional[str] = None
    vh_9_stated_germline_humanized_raw: Optional[str] = None

    # VL sequences
    vl_0_sapiens_on_mouse:              Optional[str] = None
    vl_0_sapiens_on_mouse_raw:          Optional[str] = None
    vl_1_pipeline_grafted:              Optional[str] = None
    vl_2_pipeline_humanized:            Optional[str] = None
    vl_3_lab_grafted:                   Optional[str] = None
    vl_4_detected_grafted:              Optional[str] = None
    vl_5_lab_final:                     Optional[str] = None
    vl_6_detected_humanized:            Optional[str] = None
    vl_7_detected_direct_backmut:       Optional[str] = None
    vl_8_stated_germline_grafted:       Optional[str] = None
    vl_9_stated_germline_humanized:     Optional[str] = None

    # VL Sapiens raw output (before CDR restoration)
    vl_2_pipeline_humanized_raw:        Optional[str] = None
    vl_6_detected_humanized_raw:        Optional[str] = None
    vl_9_stated_germline_humanized_raw: Optional[str] = None

    # Metadata
    vh_pipeline_germline:  Optional[str] = None
    vl_pipeline_germline:  Optional[str] = None
    vh_detected_germline:  Optional[str] = None
    vl_detected_germline:  Optional[str] = None
    vh_stated_germline:    Optional[str] = None  # from CSV vh_germline column
    vl_stated_germline:    Optional[str] = None  # from CSV vl_germline column
    vl_chain_type:         Optional[str] = None
    # CDR definition selected per sequence group (imgt/kabat/chothia)
    vh_cdr_def_pipeline:   Optional[str] = None  # used for seqs 1, 2
    vl_cdr_def_pipeline:   Optional[str] = None
    vh_cdr_def_detected:   Optional[str] = None  # used for seqs 4, 6, 7
    vl_cdr_def_detected:   Optional[str] = None
    vh_cdr_def_stated:     Optional[str] = None  # used for seqs 8, 9
    vl_cdr_def_stated:     Optional[str] = None
    error:                 Optional[str] = None

    # Germline sequences by region — stored as {region: {imgt_pos: aa}} dicts
    # Retrieved from ANARCI database at generation time using the same mechanism
    # as evaluate.py so the report can render side-by-side comparisons without
    # needing database access.
    vh_germ_pipeline_seq:  Optional[dict] = None  # pipeline top-1 germline
    vl_germ_pipeline_seq:  Optional[dict] = None
    vh_germ_detected_seq:  Optional[dict] = None  # detected lab germline
    vl_germ_detected_seq:  Optional[dict] = None
    vh_germ_stated_seq:    Optional[dict] = None  # lab-stated germline
    vl_germ_stated_seq:    Optional[dict] = None


# ── CDR definition selection ─────────────────────────────────────────────────

def select_cdr_definition(mouse_seq: str, germline_name: str,
                          chain_type: str, scheme: str = "imgt") -> str:
    """
    Select the CDR definition (imgt/kabat/chothia) that produces the longest
    total CDR length when grafting mouse_seq onto germline_name.

    This mirrors the lab's strategy of choosing whichever numbering scheme
    preserves the most CDR residues, minimizing the risk of truncating
    important CDR loop residues during humanization.

    Returns the name of the selected definition.
    """
    best_def = "kabat"  # fallback default
    max_cdr_len = -1

    for cdr_def in ["imgt", "kabat", "chothia"]:
        try:
            normalized = normalize_germline_name(germline_name)
            mouse_chain = Chain(mouse_seq, scheme=scheme,
                                cdr_definition=cdr_def)
            grafted = mouse_chain.graft_cdrs_onto_human_germline(
                v_gene=normalized, backmutate_vernier=False)

            # Count CDR residues using abnumber's built-in CDR dicts
            # These reflect the chosen cdr_definition's boundaries
            cdr_len = sum(
                1 for pos, aa in mouse_chain
                if (mouse_chain.cdr1_dict and pos in mouse_chain.cdr1_dict) or
                   (mouse_chain.cdr2_dict and pos in mouse_chain.cdr2_dict) or
                   (mouse_chain.cdr3_dict and pos in mouse_chain.cdr3_dict)
            )

            if cdr_len > max_cdr_len:
                max_cdr_len = cdr_len
                best_def = cdr_def

        except Exception:
            continue

    return best_def


# ── Core grafting helper ──────────────────────────────────────────────────────

def graft(mouse_seq: str, germline_name: str,
          scheme: str = "imgt", cdr_definition: str = None,
          chain_type: str = None) -> tuple[Optional[str], str]:
    """
    Graft mouse CDRs onto a human germline.
    If cdr_definition is None, selects the definition producing the longest CDRs.
    chain_type ('H', 'K', 'L') must be provided to avoid auto-detection errors.
    Returns (grafted_sequence, cdr_definition_used).
    """
    try:
        normalized = normalize_germline_name(germline_name)

        # Resolve chain type if not provided
        if chain_type is None:
            try:
                from pipeline.step_a_numbering import number_sequence
                num = number_sequence(mouse_seq, chain_type=None)
                chain_type = num["chain_type"]
            except Exception:
                chain_type = "H"

        if cdr_definition is None:
            cdr_definition = select_cdr_definition(
                mouse_seq, germline_name, chain_type, scheme)

        mouse_chain = Chain(mouse_seq, scheme=scheme,
                            cdr_definition=cdr_definition)
        grafted = mouse_chain.graft_cdrs_onto_human_germline(
            v_gene=normalized, backmutate_vernier=False)
        return grafted.seq, cdr_definition

    except Exception as e:
        print(f"    Grafting failed for {germline_name}: {e}")
        return None, cdr_definition or "unknown"


# ── Position map helpers ─────────────────────────────────────────────────────

def _sort_key(pos):
    """Sort key that handles both int positions and tuple (int, str) insertion positions."""
    if isinstance(pos, tuple):
        return (pos[0], pos[1])  # e.g. (111, 'A') sorts after (111, ' ')
    return (pos, ' ')            # int positions sort before their insertions


def _build_position_map(numbered: dict) -> dict:
    """
    Build a complete {position: aa} map including CDR3 insertion positions.
    Includes ALL positions — both integer IMGT positions and tuple insertion positions.
    """
    result = dict(numbered["fr_residues"])
    result.update(numbered["cdr_residues"])  # includes tuple-keyed insertions
    return result


def _reconstruct_sequence(pos_map: dict) -> str:
    """
    Reconstruct sequence from position map in correct IMGT order.
    Handles both integer positions and tuple CDR3 insertion positions.
    """
    return "".join(
        pos_map[pos]
        for pos in sorted(pos_map.keys(), key=_sort_key)
        if pos_map.get(pos)
    )


# ── Sapiens humanization helper ───────────────────────────────────────────────

def humanize_sapiens(grafted_seq: str, chain_type: str) -> tuple[Optional[str], Optional[str]]:
    """
    Apply Sapiens humanization to a grafted sequence.
    Sapiens scores every position independently — including CDRs.
    After humanization, CDR positions are restored from the input sequence
    so that only FR positions are modified.

    Returns:
        (cdr_restored_seq, raw_sapiens_seq)
        cdr_restored_seq: FRs humanized by Sapiens, CDRs preserved from grafted input
        raw_sapiens_seq:  full Sapiens output before CDR restoration (for comparison)

    predict_scores(seq, chain_type) returns a DataFrame — take idxmax per row.
    Requires biophi: pip install biophi
    """
    try:
        from sapiens import predict_scores

        # Run Sapiens on full sequence
        scores_df = predict_scores(grafted_seq, chain_type)
        sapiens_seq = "".join(scores_df.idxmax(axis=1).tolist())

        if len(sapiens_seq) != len(grafted_seq):
            print(f"    Sapiens length mismatch — returning original")
            return grafted_seq, sapiens_seq

        # Number the grafted sequence to identify which STRING INDICES are CDR
        # We do NOT number the Sapiens output — Sapiens may reassign insertion
        # positions differently, causing residue swaps at CDR3 insertions.
        # Instead we identify CDR string positions directly from the grafted seq.
        grafted_num = number_sequence(grafted_seq, chain_type=chain_type)

        # Map IMGT positions to string indices in grafted_seq
        # Build ordered list of (imgt_pos, string_index, aa) for all positions
        grafted_all = _build_position_map(grafted_num)
        ordered_positions = sorted(grafted_all.keys(), key=_sort_key)

        # Identify which string indices correspond to CDR positions
        cdr_string_indices = set()
        for str_idx, pos in enumerate(ordered_positions):
            if pos in grafted_num["cdr_residues"]:
                cdr_string_indices.add(str_idx)

        # Restore CDR positions in Sapiens output at string level
        # This avoids any IMGT renumbering issues with CDR3 insertions
        sapiens_list = list(sapiens_seq)
        grafted_list = list(grafted_seq)

        for str_idx in cdr_string_indices:
            if str_idx < len(sapiens_list) and str_idx < len(grafted_list):
                sapiens_list[str_idx] = grafted_list[str_idx]

        cdr_restored_seq = "".join(sapiens_list)

        # Count changes at string level
        fr_string_indices = set(range(len(grafted_list))) - cdr_string_indices
        fr_changes = sum(
            1 for i in fr_string_indices
            if i < len(grafted_list) and i < len(sapiens_list)
            # compare against raw Sapiens
            and grafted_list[i] != sapiens_seq[i]
        )
        cdr_changes_by_sapiens = sum(
            1 for i in cdr_string_indices
            if i < len(grafted_list) and i < len(sapiens_seq)
            and grafted_list[i] != sapiens_seq[i]
        )
        print(f"    Sapiens: {fr_changes} FR position(s) humanized, "
              f"{cdr_changes_by_sapiens} CDR position(s) restored")

        return cdr_restored_seq, sapiens_seq

    except ImportError:
        print("    Sapiens not available — install biophi: pip install biophi")
        return None, None
    except Exception as e:
        print(f"    Sapiens humanization failed: {e}")
        return None, None


# ── Direct back-mutation helper ───────────────────────────────────────────────

def apply_direct_backmutations(
    base_seq: str, reference_seq: str, source_seq: str, chain_type: str,
) -> Optional[str]:
    """
    Apply back-mutations to base_seq at FR positions where reference_seq != source_seq.
    Sequence 7: base=seq4, reference=seq3, source=seq5
    """
    try:
        numbered_base = number_sequence(base_seq,      chain_type=chain_type)
        numbered_ref = number_sequence(reference_seq, chain_type=chain_type)
        numbered_src = number_sequence(source_seq,    chain_type=chain_type)

        backmut_positions = {}
        for pos in numbered_ref["fr_residues"]:
            ref_aa = numbered_ref["fr_residues"].get(pos)
            src_aa = numbered_src["fr_residues"].get(pos)
            if ref_aa and src_aa and ref_aa != src_aa:
                backmut_positions[pos] = src_aa

        if not backmut_positions:
            print(f"    No back-mutation positions found between seq3 and seq5")
            return base_seq

        print(f"    Found {len(backmut_positions)} back-mutation positions: "
              f"{sorted(backmut_positions.keys())}")

        # Apply back-mutations at string level to preserve CDR3 insertion order
        # Using IMGT position→string index mapping from base sequence
        base_all = _build_position_map(numbered_base)
        ordered_pos = sorted(base_all.keys(), key=_sort_key)

        # Build string-level mutation map: string_index → new_aa
        str_mutations = {}
        for imgt_pos, new_aa in backmut_positions.items():
            if imgt_pos in base_all:
                str_idx = ordered_pos.index(imgt_pos)
                str_mutations[str_idx] = new_aa

        result_list = list(base_seq)
        for str_idx, new_aa in str_mutations.items():
            if str_idx < len(result_list):
                result_list[str_idx] = new_aa

        return "".join(result_list)
    except Exception as e:
        print(f"    Direct back-mutation failed: {e}")
        return None


# ── Main sequence generation function ────────────────────────────────────────

def generate_sequences(
    clone_id:        str,
    mouse_vh:        str,
    mouse_vl:        str,
    lab_hu_vh:       str,
    lab_hu_vl:       str,
    lab_final_vh:    str,
    lab_final_vl:    str,
    lab_vh_germline: Optional[str] = None,   # from CSV vh_germline column
    lab_vl_germline: Optional[str] = None,   # from CSV vl_germline column
    top_n:           int = 20,
    cdr_definition:  str = "kabat",
) -> CloneSequences:
    """Generate all 9 evaluation sequences for one antibody clone."""
    result = CloneSequences(clone_id=clone_id)

    try:
        print(f"\n  Generating sequences for {clone_id}...")

        # ── Step A: number mouse sequences ────────────────────────────────────
        vh_numbered = number_sequence(mouse_vh, chain_type="H")
        vl_numbered = number_sequence(mouse_vl, chain_type=None)
        vl_chain_type = vl_numbered["chain_type"]
        result.vl_chain_type = vl_chain_type

        # ── Step B: rank germlines from mouse FRs ─────────────────────────────
        vh_rankings = rank_germlines(
            vh_numbered["fr_residues"], "H", top_n=top_n)
        vl_rankings = rank_germlines(
            vl_numbered["fr_residues"], vl_chain_type, top_n=top_n)

        if not vh_rankings or not vl_rankings:
            result.error = "No germline candidates found"
            return result

        pipe_vh_germ = vh_rankings[0]["germline"]
        pipe_vl_germ = vl_rankings[0]["germline"]
        result.vh_pipeline_germline = pipe_vh_germ
        result.vl_pipeline_germline = pipe_vl_germ
        print(
            f"    Pipeline VH top-1: {pipe_vh_germ} ({vh_rankings[0]['fr_identity']:.1%})")
        print(
            f"    Pipeline VL top-1: {pipe_vl_germ} ({vl_rankings[0]['fr_identity']:.1%})")

        print_normalization_report(vh_rankings, "H")
        print_normalization_report(vl_rankings, vl_chain_type)

        # ── Detect lab germline from Hu sequence ──────────────────────────────
        det_vh_germ, det_vh_id, _ = detect_germline(lab_hu_vh, "H")
        det_vl_germ, det_vl_id, _ = detect_germline(lab_hu_vl, vl_chain_type)
        result.vh_detected_germline = det_vh_germ
        result.vl_detected_germline = det_vl_germ
        print(f"    Detected VH germline: {det_vh_germ} ({det_vh_id:.1%})")
        print(f"    Detected VL germline: {det_vl_germ} ({det_vl_id:.1%})")

        # ── Parse lab-stated germline from CSV ────────────────────────────────
        stated_vh_germ = parse_lab_germline(
            lab_vh_germline) if lab_vh_germline else None
        stated_vl_germ = parse_lab_germline(
            lab_vl_germline) if lab_vl_germline else None
        result.vh_stated_germline = stated_vh_germ
        result.vl_stated_germline = stated_vl_germ
        if stated_vh_germ:
            print(f"    Lab-stated VH germline: {stated_vh_germ}")
        if stated_vl_germ:
            print(f"    Lab-stated VL germline: {stated_vl_germ}")

        # ── Store germline sequences from ANARCI database ────────────────────
        print(f"    Retrieving germline sequences from database...")
        result.vh_germ_pipeline_seq = get_germline_seq_by_region(
            pipe_vh_germ, "H")
        result.vl_germ_pipeline_seq = get_germline_seq_by_region(
            pipe_vl_germ, vl_chain_type)
        if det_vh_germ:
            result.vh_germ_detected_seq = get_germline_seq_by_region(
                det_vh_germ, "H")
        if det_vl_germ:
            result.vl_germ_detected_seq = get_germline_seq_by_region(
                det_vl_germ, vl_chain_type)
        if stated_vh_germ:
            result.vh_germ_stated_seq = get_germline_seq_by_region(
                stated_vh_germ, "H")
        if stated_vl_germ:
            result.vl_germ_stated_seq = get_germline_seq_by_region(
                stated_vl_germ, vl_chain_type)

        # ── Sequence 0: Sapiens applied directly to mouse (no grafting) ─────
        # This is the pure model output — no human germline selection,
        # no CDR grafting. Sapiens predicts the most human-like residue
        # at every position starting from the raw mouse sequence.
        # CDR positions are restored after Sapiens runs (same as seqs 2, 6, 9).
        # Post-hoc germline detection on seq 0 will identify which human
        # germline Sapiens converged toward (handled in score_sequences.py).
        print(f"    Generating seq 0 (Sapiens directly on mouse — pure model output)...")
        result.vh_0_sapiens_on_mouse, result.vh_0_sapiens_on_mouse_raw = (
            humanize_sapiens(mouse_vh, "H"))
        result.vl_0_sapiens_on_mouse, result.vl_0_sapiens_on_mouse_raw = (
            humanize_sapiens(mouse_vl, vl_chain_type))

        # ── Sequence 1: pipeline grafted ──────────────────────────────────────
        print(
            f"    Generating seq 1 (pipeline grafted — selecting longest CDR definition)...")
        result.vh_1_pipeline_grafted, result.vh_cdr_def_pipeline = graft(
            mouse_vh, pipe_vh_germ, chain_type="H")
        result.vl_1_pipeline_grafted, result.vl_cdr_def_pipeline = graft(
            mouse_vl, pipe_vl_germ, chain_type=vl_chain_type)
        print(f"    VH CDR definition selected: {result.vh_cdr_def_pipeline}")
        print(f"    VL CDR definition selected: {result.vl_cdr_def_pipeline}")

        # ── Sequence 2: pipeline humanized (Sapiens) ──────────────────────────
        print(f"    Generating seq 2 (pipeline humanized via Sapiens)...")
        if result.vh_1_pipeline_grafted:
            result.vh_2_pipeline_humanized, result.vh_2_pipeline_humanized_raw = (
                humanize_sapiens(result.vh_1_pipeline_grafted, "H"))
        if result.vl_1_pipeline_grafted:
            result.vl_2_pipeline_humanized, result.vl_2_pipeline_humanized_raw = (
                humanize_sapiens(result.vl_1_pipeline_grafted, vl_chain_type))

        # ── Sequence 3: lab grafted (from CSV — ground truth) ─────────────────
        print(f"    Setting seq 3 (lab grafted — from CSV)...")
        result.vh_3_lab_grafted = lab_hu_vh
        result.vl_3_lab_grafted = lab_hu_vl

        # ── Sequence 4: detected germline grafted ─────────────────────────────
        print(f"    Generating seq 4 (detected germline grafted — selecting longest CDR definition)...")
        if det_vh_germ:
            result.vh_4_detected_grafted, result.vh_cdr_def_detected = graft(
                mouse_vh, det_vh_germ, chain_type="H")
            print(
                f"    VH CDR definition selected: {result.vh_cdr_def_detected}")
        if det_vl_germ:
            result.vl_4_detected_grafted, result.vl_cdr_def_detected = graft(
                mouse_vl, det_vl_germ, chain_type=vl_chain_type)
            print(
                f"    VL CDR definition selected: {result.vl_cdr_def_detected}")

        # ── Sequence 5: lab final (from CSV — ground truth) ───────────────────
        print(f"    Setting seq 5 (lab final — from CSV)...")
        result.vh_5_lab_final = lab_final_vh
        result.vl_5_lab_final = lab_final_vl

        # ── Sequence 6: detected germline + Sapiens ───────────────────────────
        print(f"    Generating seq 6 (detected germline + Sapiens)...")
        if result.vh_4_detected_grafted:
            result.vh_6_detected_humanized, result.vh_6_detected_humanized_raw = (
                humanize_sapiens(result.vh_4_detected_grafted, "H"))
        if result.vl_4_detected_grafted:
            result.vl_6_detected_humanized, result.vl_6_detected_humanized_raw = (
                humanize_sapiens(result.vl_4_detected_grafted, vl_chain_type))

        # ── Sequence 7: detected germline + direct back-mutations ─────────────
        print(f"    Generating seq 7 (detected germline + direct back-mutations)...")
        if result.vh_4_detected_grafted:
            result.vh_7_detected_direct_backmut = apply_direct_backmutations(
                base_seq=result.vh_4_detected_grafted,
                reference_seq=lab_hu_vh, source_seq=lab_final_vh, chain_type="H")
        if result.vl_4_detected_grafted:
            result.vl_7_detected_direct_backmut = apply_direct_backmutations(
                base_seq=result.vl_4_detected_grafted,
                reference_seq=lab_hu_vl, source_seq=lab_final_vl, chain_type=vl_chain_type)

        # ── Sequence 8: lab-stated germline grafted (from database) ───────────
        # Uses the germline sequence from the ANARCI database, not the actual input
        print(f"    Generating seq 8 (lab-stated germline grafted — selecting longest CDR definition)...")
        if stated_vh_germ:
            result.vh_8_stated_germline_grafted, result.vh_cdr_def_stated = graft(
                mouse_vh, stated_vh_germ, chain_type="H")
            print(
                f"    VH CDR definition selected: {result.vh_cdr_def_stated}")
        else:
            print(f"    Skipping seq 8 VH — no lab-stated germline in CSV")
        if stated_vl_germ:
            result.vl_8_stated_germline_grafted, result.vl_cdr_def_stated = graft(
                mouse_vl, stated_vl_germ, chain_type=vl_chain_type)
            print(
                f"    VL CDR definition selected: {result.vl_cdr_def_stated}")
        else:
            print(f"    Skipping seq 8 VL — no lab-stated germline in CSV")

        # ── Sequence 9: lab-stated germline + Sapiens ─────────────────────────
        print(f"    Generating seq 9 (lab-stated germline + Sapiens)...")
        if result.vh_8_stated_germline_grafted:
            result.vh_9_stated_germline_humanized, result.vh_9_stated_germline_humanized_raw = (
                humanize_sapiens(result.vh_8_stated_germline_grafted, "H"))
        if result.vl_8_stated_germline_grafted:
            result.vl_9_stated_germline_humanized, result.vl_9_stated_germline_humanized_raw = (
                humanize_sapiens(result.vl_8_stated_germline_grafted, vl_chain_type))

    except Exception as e:
        result.error = str(e)

    return result


# ── Output helpers ────────────────────────────────────────────────────────────

SEQ_LABELS = {
    "0":  "sapiens_on_mouse",
    "0r": "sapiens_on_mouse_raw",
    "1": "pipeline_grafted",
    "2": "pipeline_humanized",
    "3": "lab_grafted",
    "4": "detected_grafted",
    "5": "lab_final",
    "6": "detected_humanized",
    "7": "detected_direct_backmut",
    "8": "stated_germline_grafted",
    "9": "stated_germline_humanized",
    # Raw Sapiens output before CDR restoration — for comparison with CDR-restored versions
    "2r": "pipeline_humanized_raw",
    "6r": "detected_humanized_raw",
    "9r": "stated_germline_humanized_raw",
}


def _get_seq(seqs: CloneSequences, num: str, label: str, chain: str) -> Optional[str]:
    """Get sequence from CloneSequences handling the raw (r-suffixed) seq IDs.

    Raw sequences (2r, 6r, 9r) strip the 'r' from the num when building the
    field name since the field is e.g. vh_2_pipeline_humanized_raw not vh_2r_...
    """
    # Strip trailing 'r' from num — the field uses base num + full label
    # e.g. num="2r", label="pipeline_humanized_raw" → field="vh_2_pipeline_humanized_raw"
    base_num = num.rstrip("r") if num.endswith("r") else num
    field = f"{chain}_{base_num}_{label}"
    return getattr(seqs, field, None)


def print_sequences(seqs: CloneSequences) -> None:
    print(f"\n  Clone: {seqs.clone_id}")
    if seqs.error:
        print(f"  ERROR: {seqs.error}")
        return
    print(f"  {'Seq':<4} {'Label':<30} {'VH length':>10} {'VL length':>10}")
    print(f"  {'-'*4} {'-'*30} {'-'*10} {'-'*10}")
    for num, label in SEQ_LABELS.items():
        vh = _get_seq(seqs, num, label, "vh")
        vl = _get_seq(seqs, num, label, "vl")
        print(f"  {num:<4} {label:<30} {str(len(vh) if vh else 'N/A'):>10} "
              f"{str(len(vl) if vl else 'N/A'):>10}")


def export_sequences(all_seqs: list, output_path: str) -> None:
    rows = []
    for seqs in all_seqs:
        if seqs.error:
            continue
        for num, label in SEQ_LABELS.items():
            vh = _get_seq(seqs, num, label, "vh")
            vl = _get_seq(seqs, num, label, "vl")
            # Determine which CDR definition was used for this sequence group
            cdr_def_map = {
                "0":  (None, None),  # Sapiens on mouse — no CDR grafting used
                "0r": (None, None),
                "1": (seqs.vh_cdr_def_pipeline, seqs.vl_cdr_def_pipeline),
                "2": (seqs.vh_cdr_def_pipeline, seqs.vl_cdr_def_pipeline),
                "2r": (seqs.vh_cdr_def_pipeline, seqs.vl_cdr_def_pipeline),
                "3": (None, None),  # direct from CSV
                "4": (seqs.vh_cdr_def_detected, seqs.vl_cdr_def_detected),
                "5": (None, None),  # direct from CSV
                "6": (seqs.vh_cdr_def_detected, seqs.vl_cdr_def_detected),
                "6r": (seqs.vh_cdr_def_detected, seqs.vl_cdr_def_detected),
                "7": (seqs.vh_cdr_def_detected, seqs.vl_cdr_def_detected),
                "8": (seqs.vh_cdr_def_stated, seqs.vl_cdr_def_stated),
                "9": (seqs.vh_cdr_def_stated, seqs.vl_cdr_def_stated),
                "9r": (seqs.vh_cdr_def_stated, seqs.vl_cdr_def_stated),
            }
            vh_cdr_def, vl_cdr_def = cdr_def_map.get(num, (None, None))

            rows.append({
                "clone":              seqs.clone_id,
                "seq_id":             num,
                "seq_label":          label,
                "vh_sequence":        vh or "",
                "vl_sequence":        vl or "",
                "vh_germline":        seqs.vh_pipeline_germline or "",
                "vl_germline":        seqs.vl_pipeline_germline or "",
                "vh_det_germline":    seqs.vh_detected_germline or "",
                "vl_det_germline":    seqs.vl_detected_germline or "",
                "vh_stated_germline": seqs.vh_stated_germline or "",
                "vl_stated_germline": seqs.vl_stated_germline or "",
                "vl_chain_type":      seqs.vl_chain_type or "",
                "vh_cdr_def":         vh_cdr_def or "",
                "vl_cdr_def":         vl_cdr_def or "",
                # Germline sequences by region (JSON-serialised, clone-level attributes)
                # Only populated for seq_ids where the germline is relevant:
                #   pipeline: seqs 0,1,2
                #   detected: seqs 3,4,5,6,7
                #   stated:   seqs 8,9
                "vh_germ_pipeline_seq": json.dumps(seqs.vh_germ_pipeline_seq or {}),
                "vl_germ_pipeline_seq": json.dumps(seqs.vl_germ_pipeline_seq or {}),
                "vh_germ_detected_seq": json.dumps(seqs.vh_germ_detected_seq or {}),
                "vl_germ_detected_seq": json.dumps(seqs.vl_germ_detected_seq or {}),
                "vh_germ_stated_seq":   json.dumps(seqs.vh_germ_stated_seq or {}),
                "vl_germ_stated_seq":   json.dumps(seqs.vl_germ_stated_seq or {}),
            })

    if not rows:
        print("No sequences to export.")
        return

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nExported {len(rows)} sequence records to {output_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate all 9 evaluation sequences for antibody clones")
    parser.add_argument(
        "--csv",          help="Benchmark CSV path (batch mode)")
    parser.add_argument("--clone",        help="Clone ID (single mode)")
    parser.add_argument("--mouse-vh",     help="Mouse VH sequence")
    parser.add_argument("--mouse-vl",     help="Mouse VL sequence")
    parser.add_argument("--lab-hu-vh",    help="Lab Hu VH sequence (seq 3)")
    parser.add_argument("--lab-hu-vl",    help="Lab Hu VL sequence (seq 3)")
    parser.add_argument("--lab-final-vh", help="Lab final VH sequence (seq 5)")
    parser.add_argument("--lab-final-vl", help="Lab final VL sequence (seq 5)")
    parser.add_argument(
        "--vh-germline",  help="Lab-stated VH germline (seq 8/9)")
    parser.add_argument(
        "--vl-germline",  help="Lab-stated VL germline (seq 8/9)")
    parser.add_argument("--top-n",        type=int, default=20)
    parser.add_argument("--cdr-def",      default="kabat",
                        choices=["kabat", "imgt"])
    parser.add_argument("--output",       default="outputs/sequences.csv")
    args = parser.parse_args()

    all_seqs = []

    if args.csv:
        print(f"Loading clones from: {args.csv}")
        with open(args.csv, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, dialect="excel")
            rows = [r for r in reader if r.get("mouse_vh", "").strip()]
        print(f"Found {len(rows)} clones\n")

        for row in rows:
            clone_id = row.get("clone", "unknown").strip()
            seqs = generate_sequences(
                clone_id=clone_id,
                mouse_vh=row["mouse_vh"].strip(),
                mouse_vl=row["mouse_vl"].strip(),
                lab_hu_vh=row["hu_vh"].strip(),
                lab_hu_vl=row["hu_vl"].strip(),
                lab_final_vh=row["final_vh"].strip(),
                lab_final_vl=row["final_vl"].strip(),
                lab_vh_germline=row.get("vh_germline", "").strip() or None,
                lab_vl_germline=row.get("vl_germline", "").strip() or None,
                top_n=args.top_n,
                cdr_definition=args.cdr_def,
            )
            print_sequences(seqs)
            all_seqs.append(seqs)

    elif args.mouse_vh and args.mouse_vl:
        seqs = generate_sequences(
            clone_id=args.clone or "query",
            mouse_vh=args.mouse_vh,
            mouse_vl=args.mouse_vl,
            lab_hu_vh=args.lab_hu_vh or "",
            lab_hu_vl=args.lab_hu_vl or "",
            lab_final_vh=args.lab_final_vh or "",
            lab_final_vl=args.lab_final_vl or "",
            lab_vh_germline=args.vh_germline,
            lab_vl_germline=args.vl_germline,
            top_n=args.top_n,
            cdr_definition=args.cdr_def,
        )
        print_sequences(seqs)
        all_seqs.append(seqs)
    else:
        parser.error("Provide either --csv or --mouse-vh + --mouse-vl")

    if all_seqs:
        export_sequences(all_seqs, args.output)


if __name__ == "__main__":
    main()
