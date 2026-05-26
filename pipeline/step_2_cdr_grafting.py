"""
Step 2: CDR grafting and OASis humanness scoring.

For each candidate germline (VH and VL independently):
  1. Graft mouse CDRs onto the human germline framework using abnumber
  2. Score the grafted sequence with OASis humanness (if OASis DB available)
  3. Return a result table per chain with FR identity + grafted sequence + humanness

Depends on: step_a_numbering.py, step_b_germline_scoring.py
"""

from abnumber import Chain
from dataclasses import dataclass, field
from typing import Optional


# ── Data class for a single grafting result ───────────────────────────────────

@dataclass
class GraftResult:
    """Holds the result of grafting one mouse chain onto one human germline."""
    chain_type:    str           # 'H', 'K', or 'L'
    germline:      str           # e.g. 'IGHV3-23*04'
    gene:          str           # e.g. 'IGHV3-23'
    family:        str           # e.g. 'IGHV3'
    fr_identity:   float         # from Step B scoring
    grafted_seq:   str           # full grafted amino acid sequence
    fr1_seq:       str           # grafted FR1
    cdr1_seq:      str           # mouse CDR1 (preserved)
    fr2_seq:       str           # grafted FR2
    cdr2_seq:      str           # mouse CDR2 (preserved)
    fr3_seq:       str           # grafted FR3
    cdr3_seq:      str           # mouse CDR3 (preserved)
    fr4_seq:       str           # grafted FR4
    oasis_score:   Optional[float] = None  # filled in if OASis DB available
    oasis_status:  str = "not_run"         # 'scored', 'unavailable', 'error'


# ── Core grafting function ────────────────────────────────────────────────────

def graft_single(
    mouse_sequence: str,
    germline_name:  str,
    chain_type:     str,
    fr_identity:    float,
    scheme:         str = "imgt",
    cdr_definition: str = "imgt",
) -> GraftResult:
    """
    Graft mouse CDRs onto one human germline framework.

    Args:
        mouse_sequence: Raw mouse VH or VL amino acid string
        germline_name:  Target germline e.g. 'IGHV3-23*04'
        chain_type:     'H', 'K', or 'L'
        fr_identity:    FR identity score from Step B (passed through for output)
        scheme:         Numbering scheme (default: imgt)
        cdr_definition: CDR definition (default: imgt)

    Returns:
        GraftResult dataclass
    """
    # Parse mouse sequence into a Chain object
    mouse_chain = Chain(
        mouse_sequence,
        scheme=scheme,
        cdr_definition=cdr_definition,
        assign_germline=False,  # we already know it's mouse
    )

    # Graft mouse CDRs onto the specified human germline
    grafted_chain = mouse_chain.graft_cdrs_onto_human_germline(
        v_gene=germline_name,
        backmutate_vernier=False,  # straight CDR graft only — no back-mutations
    )

    # Parse gene and family from germline name e.g. 'IGHV3-23*04'
    gene = germline_name.split("*")[0]
    family = gene.split("-")[0]

    return GraftResult(
        chain_type=chain_type,
        germline=germline_name,
        gene=gene,
        family=family,
        fr_identity=fr_identity,
        grafted_seq=grafted_chain.seq,
        fr1_seq=grafted_chain.fr1_seq,
        cdr1_seq=grafted_chain.cdr1_seq,
        fr2_seq=grafted_chain.fr2_seq,
        cdr2_seq=grafted_chain.cdr2_seq,
        fr3_seq=grafted_chain.fr3_seq,
        cdr3_seq=grafted_chain.cdr3_seq,
        fr4_seq=grafted_chain.fr4_seq,
    )


def graft_candidates(
    mouse_sequence: str,
    rankings:       list[dict],
    chain_type:     str,
    mode:           str = "all",
    top_n:          int = 5,
    scheme:         str = "imgt",
    cdr_definition: str = "imgt",
) -> list[GraftResult]:
    """
    Graft mouse CDRs onto a configurable subset of ranked germline candidates.

    Args:
        mouse_sequence: Raw mouse amino acid string (VH or VL)
        rankings:       Output of rank_germlines() from Step B
        chain_type:     'H', 'K', or 'L'
        mode:           Filtering mode:
                          'all'       — graft all germlines in rankings
                          'top_n'     — graft top N by FR identity
                          'validated' — graft only therapeutically validated families
        top_n:          Used when mode='top_n'
        scheme:         Numbering scheme
        cdr_definition: CDR definition

    Returns:
        List of GraftResult, one per germline
    """
    # Filter candidates based on mode
    if mode == "top_n":
        candidates = rankings[:top_n]
    elif mode == "validated":
        candidates = [r for r in rankings if r["preferred"]]
    elif mode == "all":
        candidates = rankings
    else:
        raise ValueError(
            f"mode must be 'all', 'top_n', or 'validated'. Got: {mode}")

    if not candidates:
        raise ValueError(
            f"No candidates to graft with mode='{mode}'. "
            f"Try mode='all' or lower min_identity threshold."
        )

    results = []
    for candidate in candidates:
        try:
            result = graft_single(
                mouse_sequence=mouse_sequence,
                germline_name=candidate["germline"],
                chain_type=chain_type,
                fr_identity=candidate["fr_identity"],
                scheme=scheme,
                cdr_definition=cdr_definition,
            )
            results.append(result)
            print(f"  ✓ Grafted onto {candidate['germline']}")
        except Exception as e:
            print(f"  ✗ Failed for {candidate['germline']}: {e}")

    return results


# ── OASis scoring (optional — requires OASis DB) ──────────────────────────────

# Stringency thresholds — fraction of OAS donors a 9-mer must appear in
# to be counted as "human". Higher = more stringent.
OASIS_THRESHOLDS = {
    "loose":   0.01,   # present in ≥1% of donors
    "relaxed": 0.10,   # present in ≥10% of donors (BioPhi default)
    "medium":  0.50,   # present in ≥50% of donors
    "strict":  0.90,   # present in ≥90% of donors
}


def score_oasis(
    results:                list[GraftResult],
    oasis_db_path:          str,
    stringency:             str = "relaxed",
) -> list[GraftResult]:
    """
    Score grafted sequences with OASis humanness using BioPhi's humanness module.

    OASis chops each sequence into overlapping 9-mer peptides and checks what
    fraction of those peptides appear in natural human antibody repertoires (OAS).
    The result is a score between 0 and 1 — higher means more human-like.

    If the DB is unavailable, results are returned unchanged with
    oasis_status='unavailable' — the pipeline does not fail.

    Args:
        results:       List of GraftResult from graft_candidates()
        oasis_db_path: Path to the OASis_9mers_v1.db file (~22GB)
                       Download: wget https://zenodo.org/record/5164685/files/OASis_9mers_v1.db.gz
                                 gunzip OASis_9mers_v1.db.gz
        stringency:    How strict the humanness threshold is. One of:
                         'loose'   — peptide must appear in ≥1%  of OAS donors
                         'relaxed' — peptide must appear in ≥10% of OAS donors (default)
                         'medium'  — peptide must appear in ≥50% of OAS donors
                         'strict'  — peptide must appear in ≥90% of OAS donors

    Returns:
        Same list with oasis_score and oasis_status filled in.
        oasis_score is the OASis identity at the chosen stringency threshold.
    """
    # Validate stringency argument early
    if stringency not in OASIS_THRESHOLDS:
        raise ValueError(
            f"stringency must be one of {list(OASIS_THRESHOLDS.keys())}. Got: '{stringency}'"
        )
    min_fraction_subjects = OASIS_THRESHOLDS[stringency]

    try:
        from biophi.humanization.methods.humanness import (
            get_chain_humanness,
            OASisParams,
        )
        from abnumber import Chain as AbnumberChain
    except ImportError:
        print("  BioPhi not installed — OASis scoring unavailable.")
        print("  Install: pip install git+https://github.com/Merck/BioPhi.git")
        for r in results:
            r.oasis_status = "unavailable"
        return results

    import os
    if not os.path.exists(oasis_db_path):
        print(f"  OASis DB not found at: {oasis_db_path}")
        print(
            "  Download: wget https://zenodo.org/record/5164685/files/OASis_9mers_v1.db.gz")
        print("            gunzip OASis_9mers_v1.db.gz")
        for r in results:
            r.oasis_status = "unavailable"
        return results

    params = OASisParams(
        oasis_db_path=oasis_db_path,
        min_fraction_subjects=min_fraction_subjects,
    )

    print(
        f"  Stringency: '{stringency}' (min_fraction_subjects={min_fraction_subjects:.0%})")

    for r in results:
        try:
            # BioPhi requires an abnumber Chain object — not a raw string
            chain = AbnumberChain(
                r.grafted_seq, scheme="imgt", cdr_definition="imgt")
            humanness = get_chain_humanness(chain, params=params)

            # OASis identity = fraction of 9-mers found in >= min_fraction_subjects of OAS donors
            r.oasis_score = humanness.get_oasis_identity(min_fraction_subjects)
            r.oasis_status = "scored"
            print(f"  ✓ {r.germline}: OASis identity = {r.oasis_score:.3f}")
        except Exception as e:
            r.oasis_status = "error"
            print(f"  ✗ OASis failed for {r.germline}: {e}")

    return results


# ── Output formatting ─────────────────────────────────────────────────────────

def print_graft_results(results: list[GraftResult], chain_label: str = "") -> None:
    """Print a summary table of grafting results."""
    header = f"Grafting results"
    if chain_label:
        header += f" — {chain_label}"
    print(header)
    print("-" * 75)
    print(f"{'#':<4} {'Germline':<20} {'FR Identity':>12} {'OASis':>8}  Grafted seq (truncated)")
    print("-" * 75)
    for i, r in enumerate(results, 1):
        oasis_str = f"{r.oasis_score:.3f}" if r.oasis_score is not None else r.oasis_status
        print(
            f"{i:<4} {r.germline:<20} {r.fr_identity:>11.1%} "
            f"{oasis_str:>8}  {r.grafted_seq[:40]}..."
        )


def export_graft_results(
    vh_results: list[GraftResult],
    vl_results: list[GraftResult],
    output_path: str = "outputs/grafted_candidates.csv",
) -> None:
    """Export VH and VL graft results to a single CSV."""
    import csv
    fieldnames = [
        "chain_type", "germline", "gene", "family",
        "fr_identity", "oasis_score", "oasis_status",
        "grafted_seq",
        "fr1_seq", "cdr1_seq", "fr2_seq",
        "cdr2_seq", "fr3_seq", "cdr3_seq", "fr4_seq",
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in vh_results + vl_results:
            writer.writerow({
                "chain_type":  r.chain_type,
                "germline":    r.germline,
                "gene":        r.gene,
                "family":      r.family,
                "fr_identity": f"{r.fr_identity:.4f}",
                "oasis_score": r.oasis_score if r.oasis_score is not None else "",
                "oasis_status": r.oasis_status,
                "grafted_seq": r.grafted_seq,
                "fr1_seq":     r.fr1_seq,
                "cdr1_seq":    r.cdr1_seq,
                "fr2_seq":     r.fr2_seq,
                "cdr2_seq":    r.cdr2_seq,
                "fr3_seq":     r.fr3_seq,
                "cdr3_seq":    r.cdr3_seq,
                "fr4_seq":     r.fr4_seq,
            })
    print(f"Exported to {output_path}")


# ── Smoke test ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import os
    import sys
    _PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, _PROJECT_ROOT)
    from pipeline.step_a_numbering import number_sequence
    from pipeline.step_b_germline_scoring import rank_germlines

    test_vh = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
        "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    )
    test_vl = (
        "DIQMTQSPSSLSASVGDRVTITCKASQDVGTSVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGT"
        "DFTLTISSLQPEDFATYYCQQYYTYPPTFGQGTKVEIK"
    )

    print("Step A: Numbering...")
    vh_numbered = number_sequence(test_vh, chain_type="H")
    vl_numbered = number_sequence(test_vl, chain_type="K")

    print("\nStep B: Ranking germlines...")
    vh_rankings = rank_germlines(vh_numbered["fr_residues"], "H", top_n=5)
    vl_rankings = rank_germlines(vl_numbered["fr_residues"], "K", top_n=5)

    print("\nStep 2: Grafting CDRs...")
    print("\nVH grafting (mode='top_n', top_n=3):")
    vh_results = graft_candidates(
        mouse_sequence=test_vh,
        rankings=vh_rankings,
        chain_type="H",
        mode="top_n",
        top_n=3,
    )

    print("\nVL grafting (mode='top_n', top_n=3):")
    vl_results = graft_candidates(
        mouse_sequence=test_vl,
        rankings=vl_rankings,
        chain_type="K",
        mode="top_n",
        top_n=3,
    )

    # OASis scoring — will skip gracefully if DB not available
    # stringency options: 'loose', 'relaxed' (default), 'medium', 'strict'
    print("\nOASis scoring (will skip if DB not present):")
    _OASIS_DB = os.environ.get("OASIS_DB_PATH", os.path.join(_PROJECT_ROOT, "data", "OASis_9mers_v1.db"))
    vh_results = score_oasis(
        vh_results, oasis_db_path=_OASIS_DB, stringency="relaxed")
    vl_results = score_oasis(
        vl_results, oasis_db_path=_OASIS_DB, stringency="relaxed")

    print()
    print_graft_results(vh_results, chain_label="VH")
    print()
    print_graft_results(vl_results, chain_label="VL (kappa)")

    export_graft_results(
        vh_results, vl_results,
        output_path=os.path.join(_PROJECT_ROOT, "outputs", "grafted_candidates.csv"),
    )
