"""
Sub-step B: Score all human germlines against the mouse query's FR residues.
Returns a ranked list of top-N candidate germlines by framework identity.

Depends on: step_a_numbering.py (for number_sequence and IMGT_REGIONS)
"""

from anarci.germlines import all_germlines
from step_a_numbering import IMGT_REGIONS, ALL_FR_POSITIONS


# ── Load and parse germline database ─────────────────────────────────────────
# ANARCI ships germlines as pre-aligned strings of exactly 128 characters,
# one character per IMGT position 1-128. Gaps are '-'.
# We parse these once at import time into {germline_name: {pos: aa}} dicts.

CHAIN_TYPE_MAP = {
    "H": "H",   # heavy
    "K": "K",   # kappa light
    "L": "L",   # lambda light
}

# Therapeutically preferred V-gene families (better expression/stability)
PREFERRED_VH_FAMILIES = {"IGHV1", "IGHV3", "IGHV5"}
PREFERRED_VL_FAMILIES = {"IGKV1", "IGKV3", "IGLV1", "IGLV2"}


def _parse_germline_aligned_string(aligned: str) -> dict:
    """
    Convert a 128-character IMGT-aligned string into {imgt_position: aa}.
    Gap characters ('-') are excluded — they mean no residue at that position.
    """
    residues = {}
    for pos_idx, aa in enumerate(aligned):
        imgt_pos = pos_idx + 1  # IMGT positions are 1-indexed
        if aa != "-":
            residues[imgt_pos] = aa
    return residues


def _build_germline_fr_db(chain_type: str) -> dict:
    """
    Build a dict of {germline_name: {imgt_pos: aa}} containing only
    FR positions, for all human germlines of the given chain type.

    Args:
        chain_type: 'H', 'K', or 'L'

    Returns:
        dict of {germline_name: fr_residues_dict}
    """
    human_germlines = all_germlines["V"][chain_type]["human"]

    germline_fr_db = {}
    for name, aligned_seq in human_germlines.items():
        all_residues = _parse_germline_aligned_string(aligned_seq)
        # Keep only FR positions
        fr_residues = {
            pos: aa
            for pos, aa in all_residues.items()
            if pos in ALL_FR_POSITIONS
        }
        germline_fr_db[name] = fr_residues

    return germline_fr_db


# Pre-build databases for all three chain types at import time
_GERMLINE_FR_DB = {
    chain: _build_germline_fr_db(chain)
    for chain in ["H", "K", "L"]
}


# ── Core scoring function ─────────────────────────────────────────────────────

def compute_fr_identity(query_fr: dict, germline_fr: dict) -> tuple[float, int, int]:
    """
    Compute framework sequence identity between query and one germline.

    Only positions where BOTH sequences have a residue are counted.
    Gaps on either side are excluded from the denominator.

    Args:
        query_fr:    {imgt_pos: aa} from the mouse query (FR only)
        germline_fr: {imgt_pos: aa} from the human germline (FR only)

    Returns:
        (identity, matched_positions, comparable_positions)
        identity = matched / comparable  (float 0-1)
    """
    comparable = 0
    matched = 0

    for pos in ALL_FR_POSITIONS:
        query_aa = query_fr.get(pos)
        germline_aa = germline_fr.get(pos)

        if query_aa is None or germline_aa is None:
            continue  # gap on either side — skip

        comparable += 1
        if query_aa == germline_aa:
            matched += 1

    if comparable == 0:
        return 0.0, 0, 0

    return matched / comparable, matched, comparable


def rank_germlines(
    query_fr:   dict,
    chain_type: str,
    top_n:      int = 10,
    min_identity: float = 0.60,
) -> list[dict]:
    """
    Score all human germlines against query FR residues and return top N.

    Args:
        query_fr:     {imgt_pos: aa} FR residues from the mouse query sequence
        chain_type:   'H', 'K', or 'L'
        top_n:        how many top candidates to return
        min_identity: minimum FR identity threshold to include (default 60%)

    Returns:
        List of dicts, sorted by fr_identity descending:
            germline      - germline name e.g. 'IGHV3-23*01'
            gene          - gene-level name without allele e.g. 'IGHV3-23'
            family        - V-gene family e.g. 'IGHV3'
            fr_identity   - float 0-1
            matched       - number of identical FR positions
            comparable    - total comparable FR positions
            preferred     - bool, True if from therapeutically preferred family
    """
    if chain_type not in _GERMLINE_FR_DB:
        raise ValueError(
            f"chain_type must be 'H', 'K', or 'L', got: {chain_type}")

    germline_db = _GERMLINE_FR_DB[chain_type]
    preferred_set = PREFERRED_VH_FAMILIES if chain_type == "H" else PREFERRED_VL_FAMILIES

    scores = []
    for germline_name, germline_fr in germline_db.items():
        identity, matched, comparable = compute_fr_identity(
            query_fr, germline_fr)

        if identity < min_identity:
            continue

        # Parse gene and family from allele name e.g. 'IGHV3-23*01'
        gene = germline_name.split("*")[0]          # 'IGHV3-23'
        family = "-".join(gene.split("-")[:1])         # 'IGHV3'

        scores.append({
            "germline":    germline_name,
            "gene":        gene,
            "family":      family,
            "fr_identity": identity,
            "matched":     matched,
            "comparable":  comparable,
            "preferred":   family in preferred_set,
        })

    # Sort by identity descending, then preferred families first on ties
    scores.sort(key=lambda x: (x["fr_identity"], x["preferred"]), reverse=True)

    return scores[:top_n]


# ── Pretty printer ────────────────────────────────────────────────────────────

def print_rankings(rankings: list[dict], chain_label: str = "") -> None:
    header = f"Top {len(rankings)} germline candidates"
    if chain_label:
        header += f" for {chain_label}"
    print(header)
    print("-" * 60)
    print(f"{'Rank':<5} {'Germline':<20} {'FR Identity':>12} {'Match':>8} {'Preferred':>10}")
    print("-" * 60)
    for rank, entry in enumerate(rankings, 1):
        preferred_mark = "✓" if entry["preferred"] else ""
        print(
            f"{rank:<5} {entry['germline']:<20} "
            f"{entry['fr_identity']:>11.1%} "
            f"{entry['matched']:>3}/{entry['comparable']:<3} "
            f"{preferred_mark:>10}"
        )


# ── Smoke test ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys
    sys.path.insert(0, "/home/claude")
    from step_a_numbering import number_sequence

    test_vh = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
        "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    )
    test_vl = (
        "DIQMTQSPSSLSASVGDRVTITCKASQDVGTSVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGT"
        "DFTLTISSLQPEDFATYYCQQYYTYPPTFGQGTKVEIK"
    )

    print("Running Sub-step A...")
    vh_result = number_sequence(test_vh, chain_type="H")
    vl_result = number_sequence(test_vl, chain_type="K")

    print()
    print("Running Sub-step B: Germline scoring...")
    print()

    vh_rankings = rank_germlines(
        query_fr=vh_result["fr_residues"],
        chain_type="H",
        top_n=5,
    )
    print_rankings(vh_rankings, chain_label="VH")

    print()

    vl_rankings = rank_germlines(
        query_fr=vl_result["fr_residues"],
        chain_type="K",
        top_n=5,
    )
    print_rankings(vl_rankings, chain_label="VL (kappa)")

    print()
    print("Note: ANARCI's single best-guess VH was:", vh_result["v_gene"],
          f"({vh_result['v_identity']:.1%})")
    print("Our top-ranked VH germline is:", vh_rankings[0]["germline"],
          f"({vh_rankings[0]['fr_identity']:.1%})")
