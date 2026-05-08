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

    # Flag ties — computed across ALL scores before truncation to top_n
    # so the tie count reflects the true number of tied germlines,
    # not just those within the top_n window
    if len(scores) > 1:
        top_score = scores[0]["fr_identity"]
        top_preferred = scores[0]["preferred"]
        # Count all tied germlines in the full list (not just top_n)
        all_tied = [
            s for s in scores
            if s["fr_identity"] == top_score and s["preferred"] == top_preferred
        ]
        is_tie = len(all_tied) > 1
        n_tied = len(all_tied)
        # allele-level, deduplicated
        tied_genes = list(dict.fromkeys(s["germline"] for s in all_tied))
        for s in scores:
            if s["fr_identity"] == top_score and s["preferred"] == top_preferred:
                s["is_tie"] = is_tie
                s["n_tied"] = n_tied
                s["tied_genes"] = tied_genes
            else:
                s["is_tie"] = False
                s["n_tied"] = 0
                s["tied_genes"] = []
    else:
        for s in scores:
            s["is_tie"] = False
            s["n_tied"] = 0
            s["tied_genes"] = []

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


# ── Germline name normalization ───────────────────────────────────────────────
# ANARCI and abnumber ship different versions of the IMGT germline database.
# ANARCI has 23 kappa and 26 heavy germlines that abnumber does not recognise.
# This function resolves any ANARCI germline name to the closest abnumber name
# so grafting never fails due to a name mismatch.

def _build_normalization_map() -> dict:
    """
    Build a mapping from ANARCI germline names → closest abnumber germline name.
    Called once at import time and cached.

    Resolution strategy (in order):
      1. Exact match → use as-is
      2. Same gene, different allele → use *01 of that gene
      3. No gene match → find abnumber germline with most similar sequence
    """
    from anarci.germlines import all_germlines
    from abnumber.germlines import get_imgt_v_chains

    mapping = {}

    for chain_type in ["H", "K", "L"]:
        anarci_db = all_germlines["V"][chain_type]["human"]
        abnumber_db = get_imgt_v_chains(chain_type)
        abnumber_names = set(abnumber_db.keys())

        for anarci_name in anarci_db:
            if anarci_name in abnumber_names:
                mapping[anarci_name] = anarci_name  # exact match
                continue

            gene = anarci_name.split("*")[0]

            # Strategy 2: same gene, use *01 allele
            allele01 = f"{gene}*01"
            if allele01 in abnumber_names:
                mapping[anarci_name] = allele01
                continue

            # Strategy 3: any allele of same gene
            same_gene = [n for n in abnumber_names if n.startswith(gene + "*")]
            if same_gene:
                mapping[anarci_name] = sorted(same_gene)[0]
                continue

            # Strategy 4: sequence similarity — find abnumber germline whose
            # aligned sequence is closest to the ANARCI germline's sequence
            anarci_seq = anarci_db[anarci_name].replace("-", "")
            best_name = None
            best_score = -1
            for ab_name, ab_chain in abnumber_db.items():
                # abnumber returns Chain objects — extract sequence via iteration
                try:
                    ab_seq = "".join(aa for pos, aa in ab_chain)
                except Exception:
                    try:
                        ab_seq = ab_chain.seq if hasattr(
                            ab_chain, "seq") else str(ab_chain)
                    except Exception:
                        continue
                ab_seq = ab_seq.replace("-", "")
                # Simple overlap identity
                min_len = min(len(anarci_seq), len(ab_seq))
                if min_len == 0:
                    continue
                matches = sum(a == b for a, b in zip(anarci_seq, ab_seq))
                score = matches / min_len
                if score > best_score:
                    best_score = score
                    best_name = ab_name
            if best_name:
                mapping[anarci_name] = best_name

    return mapping


# Build map at import time
_GERMLINE_NAME_MAP = _build_normalization_map()


def normalize_germline_name(anarci_name: str) -> str:
    """
    Resolve an ANARCI germline name to the closest abnumber-compatible name.
    Handles both allele-level (IGKV1D-7-1*01) and gene-level (IGKV1D-7-1) names.
    Returns the original name if no mapping found.
    """
    # Try exact match first
    if anarci_name in _GERMLINE_NAME_MAP:
        return _GERMLINE_NAME_MAP[anarci_name]
    # Try with *01 allele suffix added (for gene-level names without allele)
    if "*" not in anarci_name:
        allele01 = f"{anarci_name}*01"
        if allele01 in _GERMLINE_NAME_MAP:
            # Return the mapped value but strip allele if input had no allele
            mapped = _GERMLINE_NAME_MAP[allele01]
            return mapped.split("*")[0] if "*" not in anarci_name else mapped
    return anarci_name

# ── Normalization comparison utility ─────────────────────────────────────────


def compare_normalized_germlines(anarci_name: str, chain_type: str) -> None:
    """
    Compare FR residues between an ANARCI germline name and its abnumber-normalized
    equivalent. Always prints the comparison.

    Shows 'identical' when names differ but FR sequences are the same — meaning
    normalization is purely a naming convention change with no biological impact.
    Shows per-position differences when FR sequences actually differ.
    """
    import sys as _sys

    # normalize_germline_name is defined in this same module
    normalized_name = normalize_germline_name(anarci_name)

    print(f"\n  Normalization: {anarci_name} -> {normalized_name}")

    if anarci_name == normalized_name:
        print(f"  Names are identical — no normalization applied")
        return

    anarci_fr = _GERMLINE_FR_DB[chain_type].get(anarci_name)
    normalized_fr = _GERMLINE_FR_DB[chain_type].get(normalized_name)

    if anarci_fr is None:
        print(f"  '{anarci_name}' not found in ANARCI database")
        return
    if normalized_fr is None:
        print(
            f"  '{normalized_name}' not found in ANARCI database — cannot compare")
        return

    all_positions = sorted(set(anarci_fr.keys()) | set(normalized_fr.keys()))
    diffs = []
    for pos in all_positions:
        aa_a = anarci_fr.get(pos)
        aa_n = normalized_fr.get(pos)
        if aa_a != aa_n:
            region = next(
                (name for name, positions in IMGT_REGIONS.items() if pos in positions),
                "unknown"
            )
            diffs.append({"pos": pos, "anarci": aa_a,
                         "normalized": aa_n, "region": region})

    if not diffs:
        print(f"  FR sequences are IDENTICAL despite different names")
        print(f"  -> Normalization has no biological impact for this germline")
        return

    print(f"  FR differences: {len(diffs)} position(s)")
    print(
        f"  {'Pos':>5}  {'Region':<8}  {anarci_name[:15]:>15}  {normalized_name[:15]:>15}")
    print(f"  {'-'*5}  {'-'*8}  {'-'*15}  {'-'*15}")
    for d in diffs:
        aa_a = d["anarci"] or "-"
        aa_n = d["normalized"] or "-"
        print(f"  {d['pos']:>5}  {d['region']:<8}  {aa_a:>15}  {aa_n:>15}")

    if len(diffs) <= 2:
        print(f"  -> Minor difference — normalization is a safe substitution")
    elif len(diffs) <= 5:
        print(f"  -> Moderate difference — verify normalization is acceptable")
    else:
        print(f"  Warning: Significant difference ({len(diffs)} positions) "
              f"— normalization may affect results")


def print_normalization_report(rankings: list[dict], chain_type: str) -> None:
    """
    Print normalization report — only for germlines where ANARCI and abnumber
    names differ. Skips identical names silently.
    """
    mismatches = [
        entry for entry in rankings
        if normalize_germline_name(entry["germline"]) != entry["germline"]
    ]

    if not mismatches:
        return  # all names identical — nothing to report

    print(f"\n{chr(61)*65}")
    print(f"NORMALIZATION REPORT ({chain_type} chain) "
          f"— {len(mismatches)}/{len(rankings)} names required normalization")
    print(f"{chr(61)*65}")
    for entry in mismatches:
        compare_normalized_germlines(entry["germline"], chain_type)


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
