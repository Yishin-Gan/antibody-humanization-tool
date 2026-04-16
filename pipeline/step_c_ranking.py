"""
Sub-step C: Combine VH and VL germline rankings into a full candidate matrix.
Produces a ranked pair table with scores shown separately, not combined.

Depends on: step_a_numbering.py, step_b_germline_scoring.py
"""

import csv
import io
from itertools import product
from step_b_germline_scoring import rank_germlines


# ── Known therapeutically validated VH/VL family pairings ────────────────────
# Based on frequency in approved humanized antibodies (TheraSAbDab analysis)
VALIDATED_FAMILY_PAIRS = {
    ("IGHV1", "IGKV1"), ("IGHV1", "IGKV3"),
    ("IGHV3", "IGKV1"), ("IGHV3", "IGKV3"),
    ("IGHV5", "IGKV1"), ("IGHV5", "IGKV3"),
    ("IGHV1", "IGLV1"), ("IGHV1", "IGLV2"),
    ("IGHV3", "IGLV1"), ("IGHV3", "IGLV2"),
    ("IGHV5", "IGLV1"), ("IGHV5", "IGLV2"),
}


def build_pair_matrix(
    vh_rankings: list[dict],
    vl_rankings: list[dict],
) -> list[dict]:
    """
    Combine VH and VL ranked lists into a full N×M candidate pair matrix.

    Each pair entry contains scores for both chains separately.
    Pairs are sorted by:
      1. Validated family pairing (True first)
      2. VH FR identity descending
      3. VL FR identity descending

    Args:
        vh_rankings: output of rank_germlines() for VH
        vl_rankings: output of rank_germlines() for VL (kappa or lambda)

    Returns:
        List of pair dicts, one per VH×VL combination.
    """
    pairs = []

    for vh, vl in product(vh_rankings, vl_rankings):
        vh_family = vh["family"]
        vl_family = vl["family"]
        is_validated = (vh_family, vl_family) in VALIDATED_FAMILY_PAIRS

        pairs.append({
            # VH fields
            "vh_germline":   vh["germline"],
            "vh_gene":       vh["gene"],
            "vh_family":     vh_family,
            "vh_fr_identity": vh["fr_identity"],
            "vh_matched":    vh["matched"],
            "vh_comparable": vh["comparable"],
            "vh_preferred":  vh["preferred"],

            # VL fields
            "vl_germline":   vl["germline"],
            "vl_gene":       vl["gene"],
            "vl_family":     vl_family,
            "vl_fr_identity": vl["fr_identity"],
            "vl_matched":    vl["matched"],
            "vl_comparable": vl["comparable"],
            "vl_preferred":  vl["preferred"],

            # Pair-level annotation
            "validated_pairing": is_validated,
        })

    # Sort: validated pairs first, then VH identity, then VL identity
    pairs.sort(
        key=lambda x: (
            x["validated_pairing"],
            x["vh_fr_identity"],
            x["vl_fr_identity"],
        ),
        reverse=True,
    )

    # Add rank after sorting
    for rank, pair in enumerate(pairs, 1):
        pair["rank"] = rank

    return pairs


def print_pair_matrix(pairs: list[dict]) -> None:
    """Print the candidate pair matrix as a formatted table."""
    print(
        f"{'Rank':<5} "
        f"{'VH Germline':<20} {'VH FR%':>7} "
        f"{'VL Germline':<20} {'VL FR%':>7} "
        f"{'Validated':>10}"
    )
    print("-" * 75)
    for p in pairs:
        validated_mark = "✓" if p["validated_pairing"] else ""
        print(
            f"{p['rank']:<5} "
            f"{p['vh_germline']:<20} {p['vh_fr_identity']:>6.1%} "
            f"{p['vl_germline']:<20} {p['vl_fr_identity']:>6.1%} "
            f"{validated_mark:>10}"
        )


def export_to_csv(pairs: list[dict]) -> str:
    """
    Export the pair matrix to a CSV string.
    Call this to write to file:  open('output.csv','w').write(export_to_csv(pairs))
    """
    fieldnames = [
        "rank",
        "vh_germline", "vh_gene", "vh_family",
        "vh_fr_identity", "vh_matched", "vh_comparable", "vh_preferred",
        "vl_germline", "vl_gene", "vl_family",
        "vl_fr_identity", "vl_matched", "vl_comparable", "vl_preferred",
        "validated_pairing",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    for pair in pairs:
        row = {k: pair[k] for k in fieldnames}
        row["vh_fr_identity"] = f"{row['vh_fr_identity']:.4f}"
        row["vl_fr_identity"] = f"{row['vl_fr_identity']:.4f}"
        writer.writerow(row)
    return output.getvalue()


# ── Full pipeline runner ──────────────────────────────────────────────────────

def run_pipeline(
    vh_sequence: str,
    vl_sequence: str,
    vl_chain_type: str = "K",
    top_n_vh: int = 5,
    top_n_vl: int = 5,
    min_identity: float = 0.60,
    export_csv: bool = False,
    csv_path: str = "candidates.csv",
) -> list[dict]:
    """
    Full Sub-step A → B → C pipeline for one antibody.

    Args:
        vh_sequence:   Raw mouse VH amino acid string
        vl_sequence:   Raw mouse VL amino acid string
        vl_chain_type: 'K' for kappa, 'L' for lambda
        top_n_vh:      Number of top VH germlines to include
        top_n_vl:      Number of top VL germlines to include
        min_identity:  Minimum FR identity threshold
        export_csv:    If True, write results to csv_path
        csv_path:      Output CSV file path

    Returns:
        List of ranked pair dicts (the full candidate matrix)
    """
    import sys
    sys.path.insert(0, "/home/claude")
    from step_a_numbering import number_sequence

    # Sub-step A: number both sequences
    print("Step A: Numbering sequences...")
    vh_result = number_sequence(vh_sequence, chain_type="H")
    vl_result = number_sequence(vl_sequence, chain_type=vl_chain_type)

    print(
        f"  VH: {vh_result['chain_type']} — ANARCI top hit: {vh_result['v_gene']} ({vh_result['v_identity']:.1%})")
    print(
        f"  VL: {vl_result['chain_type']} — ANARCI top hit: {vl_result['v_gene']} ({vl_result['v_identity']:.1%})")

    # Sub-step B: score all germlines
    print(
        f"\nStep B: Scoring germlines (top {top_n_vh} VH, top {top_n_vl} VL)...")
    vh_rankings = rank_germlines(
        query_fr=vh_result["fr_residues"],
        chain_type="H",
        top_n=top_n_vh,
        min_identity=min_identity,
    )
    vl_rankings = rank_germlines(
        query_fr=vl_result["fr_residues"],
        chain_type=vl_chain_type,
        top_n=top_n_vl,
        min_identity=min_identity,
    )

    print(f"  VH candidates: {len(vh_rankings)}")
    print(f"  VL candidates: {len(vl_rankings)}")

    # Sub-step C: build pair matrix
    print(
        f"\nStep C: Building {len(vh_rankings)}×{len(vl_rankings)} pair matrix...")
    pairs = build_pair_matrix(vh_rankings, vl_rankings)

    if export_csv:
        csv_str = export_to_csv(pairs)
        with open(csv_path, "w") as f:
            f.write(csv_str)
        print(f"  Exported to {csv_path}")

    return pairs


# ── Smoke test ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    test_vh = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
        "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    )
    test_vl = (
        "DIQMTQSPSSLSASVGDRVTITCKASQDVGTSVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGT"
        "DFTLTISSLQPEDFATYYCQQYYTYPPTFGQGTKVEIK"
    )

    pairs = run_pipeline(
        vh_sequence=test_vh,
        vl_sequence=test_vl,
        vl_chain_type="K",
        top_n_vh=5,
        top_n_vl=5,
        export_csv=True,
        csv_path="/home/claude/candidates.csv",
    )

    print()
    print(f"Total pairs in matrix: {len(pairs)}")
    print()
    print_pair_matrix(pairs)

    print()
    validated = [p for p in pairs if p["validated_pairing"]]
    print(f"Validated family pairings: {len(validated)}/{len(pairs)}")
