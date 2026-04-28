"""
Sub-step A: Parse and number an antibody sequence using ANARCI (IMGT scheme).
Separates FR and CDR positions for downstream germline scoring and CDR grafting.
"""

from anarci import anarci


# ── IMGT position boundaries ──────────────────────────────────────────────────
IMGT_REGIONS = {
    "FR1":  set(range(1,   27)),   # positions 1-26
    "CDR1": set(range(27,  39)),   # positions 27-38
    "FR2":  set(range(39,  56)),   # positions 39-55
    "CDR2": set(range(56,  66)),   # positions 56-65
    "FR3":  set(range(66,  105)),  # positions 66-104
    "CDR3": set(range(105, 118)),  # positions 105-117
    "FR4":  set(range(118, 129)),  # positions 118-128
}

ALL_FR_POSITIONS = set().union(
    *(v for k, v in IMGT_REGIONS.items() if k.startswith("FR")))
ALL_CDR_POSITIONS = set().union(
    *(v for k, v in IMGT_REGIONS.items() if k.startswith("CDR")))


def number_sequence(sequence: str, chain_type: str = None,
                    cdr_definition: str = "imgt") -> dict:
    """
    Number a raw antibody sequence with ANARCI under the IMGT scheme
    and split residues into framework (FR) and CDR dictionaries.

    Args:
        sequence:   Raw amino acid string (single-letter codes, no gaps).
        chain_type: 'H' for heavy chain, 'K' for kappa, 'L' for lambda.
                    If None, ANARCI infers it automatically.

    Returns:
        dict:
            chain_type   - detected chain type ('H', 'K', or 'L')
            v_gene       - top human germline V gene hit, e.g. 'IGHV3-23*01'
            v_identity   - sequence identity to that V gene (0-1 float)
            j_gene       - top human germline J gene hit
            numbered     - raw ANARCI output: list of ((pos, ins_code), aa)
            fr_residues  - {imgt_pos: aa} for FR positions only
            cdr_residues - {imgt_pos: aa} for CDR positions only
            fr_by_region - {'FR1': {pos: aa}, 'FR2': ..., 'FR3': ..., 'FR4': ...}
            cdr_by_region- {'CDR1': {pos: aa}, 'CDR2': ..., 'CDR3': ...}
    """
    allow_set = {chain_type} if chain_type else {"H", "K", "L"}

    results, numbered_meta, details = anarci(
        [("query", sequence)],
        scheme="imgt",
        assign_germline=True,
        allowed_species=["human", "mouse"],  # accept mouse input sequence
        allow=allow_set,
    )

    if results[0] is None:
        raise ValueError(
            "ANARCI could not number this sequence. "
            "Verify it is a valid VH or VL amino acid string."
        )

    # results[0] is a list of domain hits; take the top hit
    numbered_positions, chain_id, _ = results[0][0]

    # Germline info lives in numbered_meta
    meta = numbered_meta[0][0]  # top hit metadata dict
    detected_chain = meta["chain_type"]

    v_gene_info = meta.get("germlines", {}).get("v_gene", [None, None])
    j_gene_info = meta.get("germlines", {}).get("j_gene", [None, None])

    # e.g. 'IGHV3-23*01'
    v_gene = v_gene_info[0][1] if v_gene_info[0] else None
    v_identity = v_gene_info[1] if v_gene_info[1] else None  # e.g. 0.824
    j_gene = j_gene_info[0][1] if j_gene_info[0] else None

    # ── Split into FR and CDR dicts ────────────────────────────────────────────
    fr_residues = {}
    cdr_residues = {}

    for (pos, ins_code), aa in numbered_positions:
        if aa == "-":
            continue  # gap position — no residue here

        # CDR3 can have insertions encoded as (111, 'A'), (111, 'B') etc.
        # Treat any insertion-code position as CDR regardless of numeric pos.
        if ins_code != " ":
            cdr_residues[(pos, ins_code)] = aa
        elif pos in ALL_FR_POSITIONS:
            fr_residues[pos] = aa
        elif pos in ALL_CDR_POSITIONS:
            cdr_residues[pos] = aa
        # positions outside 1-128 are rare edge cases — silently skip

    # ── Group by named region ──────────────────────────────────────────────────
    fr_by_region = {
        name: {p: fr_residues[p] for p in bounds if p in fr_residues}
        for name, bounds in IMGT_REGIONS.items()
        if name.startswith("FR")
    }

    cdr_by_region = {
        name: {p: cdr_residues.get(p) for p in bounds if p in cdr_residues}
        for name, bounds in IMGT_REGIONS.items()
        if name.startswith("CDR")
    }

    return {
        "chain_type":    detected_chain,
        "v_gene":        v_gene,
        "v_identity":    v_identity,
        "j_gene":        j_gene,
        "numbered":      numbered_positions,
        "fr_residues":   fr_residues,
        "cdr_residues":  cdr_residues,
        "fr_by_region":  fr_by_region,
        "cdr_by_region": cdr_by_region,
    }


def summarise(result: dict) -> None:
    """Print a human-readable summary of a numbering result."""
    print(f"  Chain type  : {result['chain_type']}")
    print(
        f"  V gene      : {result['v_gene']}  (identity: {result['v_identity']:.1%})")
    print(f"  J gene      : {result['j_gene']}")
    print()
    for region in ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]:
        if region.startswith("FR"):
            residues = result["fr_by_region"].get(region, {})
            seq = "".join(residues[p] for p in sorted(residues))
        else:
            residues = result["cdr_by_region"].get(region, {})
            # CDR3 may have tuple keys for insertions — sort carefully
            seq = "".join(residues[p] for p in sorted(residues) if residues[p])
        label = "FR " if region.startswith("FR") else "CDR"
        print(f"  {region} ({len(residues):>2} residues): {seq}")


# ── Smoke test ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Trastuzumab parental mouse VH and VL sequences
    test_vh = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
        "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    )
    test_vl = (
        "DIQMTQSPSSLSASVGDRVTITCKASQDVGTSVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGT"
        "DFTLTISSLQPEDFATYYCQQYYTYPPTFGQGTKVEIK"
    )

    print("=" * 60)
    print("VH NUMBERING")
    print("=" * 60)
    vh_result = number_sequence(test_vh, chain_type="H")
    summarise(vh_result)

    print()
    print("=" * 60)
    print("VL NUMBERING (kappa)")
    print("=" * 60)
    vl_result = number_sequence(test_vl, chain_type="K")
    summarise(vl_result)

    print()
    print("FR residue count  VH:", len(vh_result["fr_residues"]))
    print("CDR residue count VH:", len(vh_result["cdr_residues"]))
    print("FR residue count  VL:", len(vl_result["fr_residues"]))
    print("CDR residue count VL:", len(vl_result["cdr_residues"]))
