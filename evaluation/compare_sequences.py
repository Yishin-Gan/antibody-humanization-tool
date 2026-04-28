"""
compare_sequences.py

Side-by-side sequence comparisons for all 7 evaluation sequences.

Comparisons performed (each split by VH/VL and CDR/FR regions):
  1 vs 2  — pipeline grafted vs pipeline humanized (Sapiens effect)
  3 vs 5  — lab grafted vs lab final (lab's back-mutations)
  4 vs 6  — detected grafted vs detected humanized (Sapiens on detected germline)
  4 vs 7  — detected grafted vs detected direct back-mutated
  2 vs 5  — pipeline humanized vs lab final (pipeline quality vs lab)
  6 vs 5  — detected humanized vs lab final (detected germline + Sapiens vs lab)
  7 vs 5  — detected direct back-mutated vs lab final (direct approach vs lab)

Input:  outputs/all_sequences.csv (from generate_sequences.py)
Output: printed report + outputs/comparison_report.csv

Usage (from project root):
    python3 pipeline/compare_sequences.py \\
        --csv outputs/all_sequences.csv \\
        --output outputs/comparison_report.csv
"""

# isort: skip_file
from pipeline.step_a_numbering import number_sequence, IMGT_REGIONS
from typing import Optional
from dataclasses import dataclass, field
import csv
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))


# ── Sequence labels ───────────────────────────────────────────────────────────

SEQ_LABELS = {
    "1": "pipeline_grafted",
    "2": "pipeline_humanized",
    "3": "lab_grafted",
    "4": "detected_grafted",
    "5": "lab_final",
    "6": "detected_humanized",
    "7": "detected_direct_backmut",
}

# Comparisons to run: (seq_a_id, seq_b_id, description)
COMPARISONS = [
    ("1", "2", "Pipeline grafted vs Pipeline humanized (Sapiens effect)"),
    ("3", "5", "Lab grafted vs Lab final (lab back-mutations)"),
    ("4", "6", "Detected grafted vs Detected humanized (Sapiens on detected germline)"),
    ("4", "7", "Detected grafted vs Detected direct back-mutated"),
    ("2", "5", "Pipeline humanized vs Lab final (pipeline quality)"),
    ("6", "5", "Detected humanized vs Lab final"),
    ("7", "5", "Detected direct back-mutated vs Lab final"),
]


# ── Region extraction ─────────────────────────────────────────────────────────

def get_regions(sequence: str, chain_type: str) -> Optional[dict]:
    """
    Number a sequence and return FR and CDR regions as {region: {pos: aa}}.
    Returns None if numbering fails.
    """
    try:
        result = number_sequence(sequence, chain_type=chain_type)
        regions = {}
        for region, positions in IMGT_REGIONS.items():
            if region.startswith("FR"):
                residues = result["fr_by_region"].get(region, {})
            else:
                residues = result["cdr_by_region"].get(region, {})
            regions[region] = {p: residues[p]
                               for p in sorted(residues) if residues.get(p)}
        return regions
    except Exception as e:
        print(f"    Numbering failed: {e}")
        return None


# ── Comparison formatting ─────────────────────────────────────────────────────

def format_region_diff(
    region:  str,
    seq_a:   dict,
    seq_b:   dict,
    label_a: str,
    label_b: str,
) -> list[str]:
    """
    Format side-by-side comparison of two region sequences.
    Returns list of lines to print.
    """
    all_positions = sorted(set(seq_a.keys()) | set(seq_b.keys()))
    if not all_positions:
        return [f"    {region}: (empty in both sequences)"]

    lines = [f"    {region}:"]

    a_str = ""
    b_str = ""
    diff_str = ""

    for pos in all_positions:
        aa_a = seq_a.get(pos, "-")
        aa_b = seq_b.get(pos, "-")
        a_str += aa_a
        b_str += aa_b
        diff_str += " " if aa_a == aa_b else "^"

    # Truncate labels to fixed width for alignment
    la = f"{label_a[:12]:<12}"
    lb = f"{label_b[:12]:<12}"

    lines.append(f"      {la}: {a_str}")
    lines.append(f"      {lb}: {b_str}")

    if "^" in diff_str:
        lines.append(f"      {'Diff':<12}  {diff_str}")
        diffs = [
            f"pos{pos}({seq_a.get(pos, '-')}->{seq_b.get(pos, '-')})"
            for pos in all_positions
            if seq_a.get(pos, "-") != seq_b.get(pos, "-")
        ]
        lines.append(f"      {'Changes':<12}: {', '.join(diffs)}")
    else:
        lines.append(f"      {'Diff':<12}  (identical)")

    return lines


def compare_two_sequences(
    seq_a:      str,
    seq_b:      str,
    chain_type: str,
    label_a:    str,
    label_b:    str,
    chain_label: str,
) -> dict:
    """
    Compare two sequences across all CDR and FR regions.
    Returns comparison result dict and prints formatted output.
    """
    print(f"\n    ── {chain_label} ──────────────────────────────────────────")

    regions_a = get_regions(seq_a, chain_type)
    regions_b = get_regions(seq_b, chain_type)

    if regions_a is None or regions_b is None:
        print(f"    Cannot compare — numbering failed for one or both sequences")
        return {"error": "numbering_failed"}

    total_diffs = 0
    total_compared = 0
    cdr_diffs = 0
    fr_diffs = 0
    diff_positions = []

    # CDR regions
    print(f"\n      CDR regions:")
    for region in ["CDR1", "CDR2", "CDR3"]:
        ra = regions_a.get(region, {})
        rb = regions_b.get(region, {})
        for line in format_region_diff(region, ra, rb, label_a, label_b):
            print(line)

        # Count diffs
        all_pos = sorted(set(ra.keys()) | set(rb.keys()))
        for pos in all_pos:
            aa_a = ra.get(pos, "-")
            aa_b = rb.get(pos, "-")
            if aa_a != "-" and aa_b != "-":
                total_compared += 1
                if aa_a != aa_b:
                    total_diffs += 1
                    cdr_diffs += 1
                    diff_positions.append({
                        "position": pos, "region": region,
                        label_a: aa_a, label_b: aa_b
                    })

    # FR regions
    print(f"\n      FR regions:")
    for region in ["FR1", "FR2", "FR3", "FR4"]:
        ra = regions_a.get(region, {})
        rb = regions_b.get(region, {})
        for line in format_region_diff(region, ra, rb, label_a, label_b):
            print(line)

        # Count diffs
        all_pos = sorted(set(ra.keys()) | set(rb.keys()))
        for pos in all_pos:
            aa_a = ra.get(pos, "-")
            aa_b = rb.get(pos, "-")
            if aa_a != "-" and aa_b != "-":
                total_compared += 1
                if aa_a != aa_b:
                    total_diffs += 1
                    fr_diffs += 1
                    diff_positions.append({
                        "position": pos, "region": region,
                        label_a: aa_a, label_b: aa_b
                    })

    # Summary
    identity = (total_compared - total_diffs) / \
        total_compared if total_compared else 0
    print(f"\n      Summary: {total_diffs} differences "
          f"({cdr_diffs} CDR, {fr_diffs} FR) "
          f"— identity {identity:.1%}")

    return {
        "total_diffs":    total_diffs,
        "cdr_diffs":      cdr_diffs,
        "fr_diffs":       fr_diffs,
        "identity":       identity,
        "diff_positions": diff_positions,
    }


# ── Main comparison runner ────────────────────────────────────────────────────

def run_comparisons_for_clone(
    clone_id:  str,
    sequences: dict,
) -> list[dict]:
    """
    Run all 7 comparisons for one clone.

    Args:
        clone_id:  clone identifier
        sequences: dict of {seq_id: {"vh": seq, "vl": seq, "vl_chain_type": ct}}

    Returns:
        list of result dicts for CSV export
    """
    results = []

    print(f"\n{'='*65}")
    print(f"Clone: {clone_id}")
    print(f"{'='*65}")

    for seq_a_id, seq_b_id, description in COMPARISONS:
        seq_a = sequences.get(seq_a_id)
        seq_b = sequences.get(seq_b_id)

        if not seq_a or not seq_b:
            print(f"\n  [{seq_a_id} vs {seq_b_id}] {description}")
            print(f"  SKIPPED — one or both sequences not available")
            print(f"  (seq {seq_a_id}: {'present' if seq_a else 'MISSING'}, "
                  f"seq {seq_b_id}: {'present' if seq_b else 'MISSING'})")
            continue

        label_a = f"seq{seq_a_id}_{SEQ_LABELS[seq_a_id][:12]}"
        label_b = f"seq{seq_b_id}_{SEQ_LABELS[seq_b_id][:12]}"

        print(f"\n  [{seq_a_id} vs {seq_b_id}] {description}")
        print(f"  seq {seq_a_id}: {SEQ_LABELS[seq_a_id]}")
        print(f"  seq {seq_b_id}: {SEQ_LABELS[seq_b_id]}")

        vl_chain_type = seq_a.get("vl_chain_type", "K")

        # VH comparison
        vh_result = {}
        if seq_a.get("vh") and seq_b.get("vh"):
            vh_result = compare_two_sequences(
                seq_a["vh"], seq_b["vh"],
                chain_type="H",
                label_a=label_a, label_b=label_b,
                chain_label="VH",
            )
        else:
            print(f"\n    ── VH — MISSING sequence(s)")

        # VL comparison
        vl_result = {}
        if seq_a.get("vl") and seq_b.get("vl"):
            vl_result = compare_two_sequences(
                seq_a["vl"], seq_b["vl"],
                chain_type=vl_chain_type,
                label_a=label_a, label_b=label_b,
                chain_label=f"VL ({vl_chain_type})",
            )
        else:
            print(f"\n    ── VL — MISSING sequence(s)")

        results.append({
            "clone":        clone_id,
            "comparison":   f"{seq_a_id}_vs_{seq_b_id}",
            "seq_a":        SEQ_LABELS[seq_a_id],
            "seq_b":        SEQ_LABELS[seq_b_id],
            "description":  description,
            "vh_identity":  f"{vh_result.get('identity', ''):.1%}" if vh_result.get('identity') is not None else "",
            "vh_total_diff": vh_result.get("total_diffs", ""),
            "vh_cdr_diff":  vh_result.get("cdr_diffs", ""),
            "vh_fr_diff":   vh_result.get("fr_diffs", ""),
            "vl_identity":  f"{vl_result.get('identity', ''):.1%}" if vl_result.get('identity') is not None else "",
            "vl_total_diff": vl_result.get("total_diffs", ""),
            "vl_cdr_diff":  vl_result.get("cdr_diffs", ""),
            "vl_fr_diff":   vl_result.get("fr_diffs", ""),
        })

    return results


# ── CSV loader ────────────────────────────────────────────────────────────────

def load_sequences_csv(csv_path: str) -> dict:
    """
    Load sequences from generate_sequences.py output CSV.

    Returns:
        {clone_id: {seq_id: {"vh": str, "vl": str, "vl_chain_type": str}}}
    """
    clones = {}

    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, dialect="excel")
        for row in reader:
            clone_id = row["clone"].strip()
            seq_id = row["seq_id"].strip()
            vh = row.get("vh_sequence", "").strip()
            vl = row.get("vl_sequence", "").strip()
            chain_type = row.get("vl_chain_type", "K").strip() or "K"

            if clone_id not in clones:
                clones[clone_id] = {}

            clones[clone_id][seq_id] = {
                "vh":            vh or None,
                "vl":            vl or None,
                "vl_chain_type": chain_type,
            }

    return clones


# ── Summary table ─────────────────────────────────────────────────────────────

def print_summary_table(all_results: list[dict]) -> None:
    """Print a compact summary table of all comparisons across all clones."""
    print(f"\n{'='*65}")
    print("SUMMARY TABLE")
    print(f"{'='*65}")
    print(f"{'Clone':<10} {'Comparison':<12} "
          f"{'VH id%':>7} {'VH diff':>8} "
          f"{'VL id%':>7} {'VL diff':>8}")
    print(f"{'-'*10} {'-'*12} {'-'*7} {'-'*8} {'-'*7} {'-'*8}")

    for r in all_results:
        print(
            f"{r['clone']:<10} {r['comparison']:<12} "
            f"{r['vh_identity']:>7} {str(r['vh_total_diff']):>8} "
            f"{r['vl_identity']:>7} {str(r['vl_total_diff']):>8}"
        )


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Side-by-side sequence comparisons for all 7 evaluation sequences")
    parser.add_argument("--csv",    required=True,
                        help="Path to all_sequences.csv from generate_sequences.py")
    parser.add_argument("--output", default="outputs/comparison_report.csv",
                        help="Output CSV path for summary results")
    parser.add_argument("--clone",  default=None,
                        help="Only compare this clone (optional)")
    args = parser.parse_args()

    print(f"Loading sequences from: {args.csv}")
    clones = load_sequences_csv(args.csv)
    print(f"Found {len(clones)} clone(s)\n")

    all_results = []
    for clone_id, sequences in clones.items():
        if args.clone and clone_id != args.clone:
            continue
        results = run_comparisons_for_clone(clone_id, sequences)
        all_results.extend(results)

    # Summary table
    print_summary_table(all_results)

    # Export CSV
    if all_results:
        fieldnames = list(all_results[0].keys())
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)
        print(f"\nSummary exported to: {args.output}")


if __name__ == "__main__":
    main()
