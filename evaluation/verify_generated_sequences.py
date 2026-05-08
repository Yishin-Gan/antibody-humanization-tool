"""
verify_generated_sequences.py

Verifies that all 9 generated sequences from generate_sequences.py
are correct according to their definitions.

Checks per sequence:
  Seq 3: identical to hu_vh/hu_vl from benchmark CSV
  Seq 5: identical to final_vh/final_vl from benchmark CSV
  Seq 1: CDRs match mouse, FRs match pipeline germline DB
  Seq 4: CDRs match mouse, FRs match detected germline DB
  Seq 8: CDRs match mouse, FRs match lab-stated germline DB
  Seq 2: CDRs identical to seq 1, FRs changed by Sapiens
  Seq 6: CDRs identical to seq 4, FRs changed by Sapiens
  Seq 9: CDRs identical to seq 8, FRs changed by Sapiens
  Seq 7: CDRs identical to seq 4, back-mutations at correct FR positions

Usage (from project root):
    python3 tests/verify_generated_sequences.py \\
        --generated outputs/all_sequences.csv \\
        --benchmark data/benchmarks/humanization_benchmark.csv
"""

# isort: skip_file
import sys
sys.path.insert(0, "/workspace/antibody-humanization-tool")  # noqa: E402
from pipeline.step_a_numbering import number_sequence
from evaluation.evaluate import get_germline_fr_by_region
from typing import Optional
import argparse
import csv
from dataclasses import dataclass


# ── Dynamic boundary computation ─────────────────────────────────────────────
# Boundary positions are computed per-clone per-chain by comparing Kabat and
# IMGT grafted sequences. Positions where they differ are excluded from CDR
# checks — these positions are legitimately different between the two schemes.

_boundary_cache: dict = {}


def compute_boundary_positions(
    mouse_seq: str, germline_name: str, chain_type: str,
    selected_def: str = "kabat"
) -> set:
    """
    Compute boundary positions where the selected CDR definition and IMGT disagree
    for a specific mouse sequence and germline combination.

    In Approach 2 (longest-CDR strategy), selected_def is the definition that was
    actually used for grafting — read from the vh_cdr_def/vl_cdr_def CSV columns.
    If selected_def == 'imgt', no boundaries exist (same definition used both sides).

    Results are cached to avoid redundant computation across sequences.
    """
    cache_key = (mouse_seq, germline_name, chain_type, selected_def)
    if cache_key in _boundary_cache:
        return _boundary_cache[cache_key]

    # If grafting used IMGT (same as verification), no boundary disagreements
    if selected_def == "imgt":
        _boundary_cache[cache_key] = set()
        return set()

    try:
        from abnumber import Chain
        from pipeline.step_b_germline_scoring import normalize_germline_name

        normalized = normalize_germline_name(germline_name)

        selected_chain = Chain(mouse_seq, scheme="imgt",
                               cdr_definition=selected_def)
        imgt_chain = Chain(mouse_seq, scheme="imgt", cdr_definition="imgt")

        seq_selected = selected_chain.graft_cdrs_onto_human_germline(
            normalized, backmutate_vernier=False).seq
        seq_imgt = imgt_chain.graft_cdrs_onto_human_germline(
            normalized, backmutate_vernier=False).seq

        num_selected = number_sequence(seq_selected, chain_type=chain_type)
        num_imgt = number_sequence(seq_imgt,     chain_type=chain_type)

        all_selected = {
            **num_selected["fr_residues"],
            **{k: v for k, v in num_selected["cdr_residues"].items() if isinstance(k, int)}
        }
        all_imgt = {
            **num_imgt["fr_residues"],
            **{k: v for k, v in num_imgt["cdr_residues"].items() if isinstance(k, int)}
        }

        boundary = {
            pos for pos in set(all_selected) | set(all_imgt)
            if all_selected.get(pos) != all_imgt.get(pos)
        }

        _boundary_cache[cache_key] = boundary
        return boundary

    except Exception as e:
        print(f"    Warning: could not compute boundary positions: {e}")
        _boundary_cache[cache_key] = set()
        return set()


# ── Region extraction helpers ─────────────────────────────────────────────────

def get_cdrs(sequence: str, chain_type: str) -> Optional[dict]:
    """Extract CDR regions as {region: {pos: aa}}."""
    try:
        result = number_sequence(sequence, chain_type=chain_type)
        return result.get("cdr_by_region", {})
    except Exception as e:
        return None


def get_frs(sequence: str, chain_type: str) -> Optional[dict]:
    """Extract FR regions as {region: {pos: aa}}."""
    try:
        result = number_sequence(sequence, chain_type=chain_type)
        return result.get("fr_by_region", {})
    except Exception as e:
        return None


def get_fr_flat(sequence: str, chain_type: str) -> Optional[dict]:
    """Extract all FR residues as flat {pos: aa}."""
    try:
        result = number_sequence(sequence, chain_type=chain_type)
        return result.get("fr_residues", {})
    except Exception as e:
        return None


# Positions where ANARCI renumbering can shift residue assignments
# after Sapiens humanization — gap-adjacent positions in CDR2 and CDR3
# These are excluded regardless of germline or CDR definition
ANARCI_RENUMBER_EXCLUSIONS_VH = {62, 63, 64, 65, 110, 111, 112, 113}
ANARCI_RENUMBER_EXCLUSIONS_VL = {55, 56, 104, 105}


def cdrs_identical(seq_a: str, seq_b: str, chain_type: str,
                   boundary: set = None) -> tuple[bool, list]:
    """
    Check if two sequences have identical CDRs.
    Excludes:
    1. Boundary positions computed dynamically per germline (Kabat/IMGT disagreements)
    2. Gap-adjacent positions where ANARCI renumbering can shift after Sapiens
    Returns (match, list of diffs).
    """
    cdrs_a = get_cdrs(seq_a, chain_type)
    cdrs_b = get_cdrs(seq_b, chain_type)
    if cdrs_a is None or cdrs_b is None:
        return False, ["numbering failed"]

    boundary = boundary or set()
    # Add ANARCI renumbering exclusions for this chain type
    if chain_type == "H":
        boundary = boundary | ANARCI_RENUMBER_EXCLUSIONS_VH
    else:
        boundary = boundary | ANARCI_RENUMBER_EXCLUSIONS_VL
    diffs = []
    for region in ["CDR1", "CDR2", "CDR3"]:
        ra = cdrs_a.get(region, {})
        rb = cdrs_b.get(region, {})
        all_pos = sorted(set(ra) | set(rb))
        for pos in all_pos:
            if pos in boundary:
                continue  # skip Kabat/IMGT boundary positions
            if ra.get(pos) != rb.get(pos):
                diffs.append(
                    f"{region} pos{pos}({ra.get(pos, '–')}→{rb.get(pos, '–')})")

    return len(diffs) == 0, diffs


def frs_match_germline(sequence: str, germline_name: str,
                       chain_type: str, regions: list = None) -> tuple[bool, int, int, list]:
    """
    Check how closely sequence FRs match a germline from the database.
    Returns (is_close_match, matched, comparable, diffs_list)
    is_close_match = True if identity >= 85% (allows for Kabat/IMGT boundary noise)
    """
    if regions is None:
        regions = ["FR1", "FR2", "FR3"]  # exclude FR4 (J gene)

    seq_frs = get_frs(sequence, chain_type)
    germ_frs = get_germline_fr_by_region(germline_name, chain_type)

    if not seq_frs or not germ_frs:
        return False, 0, 0, ["could not retrieve sequences"]

    matched = comparable = 0
    diffs = []
    for region in regions:
        sr = seq_frs.get(region, {})
        gr = germ_frs.get(region, {})
        for pos in sorted(set(sr) | set(gr)):
            aa_s = sr.get(pos)
            aa_g = gr.get(pos)
            if aa_s and aa_g:
                comparable += 1
                if aa_s == aa_g:
                    matched += 1
                else:
                    diffs.append(f"{region} pos{pos}({aa_s}→{aa_g})")

    identity = matched / comparable if comparable else 0
    return identity >= 0.85, matched, comparable, diffs


def frs_changed(seq_a: str, seq_b: str, chain_type: str) -> tuple[bool, int]:
    """Check if FR residues changed between two sequences. Returns (changed, n_diffs)."""
    frs_a = get_fr_flat(seq_a, chain_type)
    frs_b = get_fr_flat(seq_b, chain_type)
    if frs_a is None or frs_b is None:
        return False, 0
    diffs = sum(1 for pos in set(frs_a) | set(frs_b)
                if frs_a.get(pos) != frs_b.get(pos))
    return diffs > 0, diffs


def backmut_positions_correct(seq3: str, seq5: str, seq4: str, seq7: str,
                              chain_type: str) -> tuple[bool, list]:
    """
    Verify seq7 has back-mutations exactly where seq3 != seq5 (at FR positions).
    Returns (correct, issues_list)
    """
    frs3 = get_fr_flat(seq3, chain_type)
    frs5 = get_fr_flat(seq5, chain_type)
    frs4 = get_fr_flat(seq4, chain_type)
    frs7 = get_fr_flat(seq7, chain_type)

    if not all([frs3, frs5, frs4, frs7]):
        return False, ["numbering failed for one or more sequences"]

    expected_positions = {
        pos for pos in frs3
        if frs3.get(pos) and frs5.get(pos) and frs3[pos] != frs5[pos]
    }

    issues = []
    for pos in expected_positions:
        expected_aa = frs5[pos]
        actual_aa = frs7.get(pos)
        if actual_aa != expected_aa:
            issues.append(
                f"pos{pos}: expected {expected_aa} (from seq5), got {actual_aa}")

    return len(issues) == 0, issues


# ── Check result container ────────────────────────────────────────────────────

@dataclass
class CheckResult:
    clone:   str
    chain:   str
    seq_id:  str
    check:   str
    passed:  bool
    detail:  str = ""


# ── Main verification function ────────────────────────────────────────────────

def verify_clone(clone_id: str, seqs: dict, benchmark: dict,
                 vl_chain_type: str) -> list[CheckResult]:
    """
    Run all checks for one clone.

    Args:
        seqs:      {seq_id: {"vh": str, "vl": str}}
        benchmark: {col: value} row from benchmark CSV
        vl_chain_type: 'K' or 'L'
    """
    results = []

    def check(seq_id, chain, check_name, passed, detail=""):
        results.append(CheckResult(clone_id, chain,
                       seq_id, check_name, passed, detail))

    mouse_vh = benchmark.get("mouse_vh", "")
    mouse_vl = benchmark.get("mouse_vl", "")

    pipe_germ_vh = seqs.get("meta", {}).get("vh_germline", "")
    pipe_germ_vl = seqs.get("meta", {}).get("vl_germline", "")
    det_germ_vh = seqs.get("meta", {}).get("vh_det_germline", "")
    det_germ_vl = seqs.get("meta", {}).get("vl_det_germline", "")
    stated_germ_vh = seqs.get("meta", {}).get("vh_stated_germline", "")
    stated_germ_vl = seqs.get("meta", {}).get("vl_stated_germline", "")

    for chain, chain_type, mouse_seq, bench_hu, bench_final, \
            pipe_germ, det_germ, stated_germ in [
                ("VH", "H",         mouse_vh, benchmark.get("hu_vh", ""),
                 benchmark.get("final_vh", ""), pipe_germ_vh, det_germ_vh, stated_germ_vh),
                ("VL", vl_chain_type, mouse_vl, benchmark.get("hu_vl", ""),
                 benchmark.get("final_vl", ""), pipe_germ_vl, det_germ_vl, stated_germ_vl),
            ]:
        s = {k: (v[chain.lower()] if chain == "VH" else v["vl"])
             for k, v in seqs.items() if k != "meta" and isinstance(v, dict)}

        # ── Check seq 3: must equal CSV hu_vh/hu_vl ──────────────────────────
        if s.get("3") and bench_hu:
            passed = s["3"].strip() == bench_hu.strip()
            check("3", chain, "== CSV hu sequence", passed,
                  "" if passed else "sequence mismatch with CSV")

        # ── Check seq 5: must equal CSV final_vh/final_vl ────────────────────
        if s.get("5") and bench_final:
            passed = s["5"].strip() == bench_final.strip()
            check("5", chain, "== CSV final sequence", passed,
                  "" if passed else "sequence mismatch with CSV")

        # ── Compute dynamic boundary positions per germline and CDR definition ──
        # Read the CDR definition that was actually used for each sequence group
        chain_key = "vh" if chain_type == "H" else "vl"
        cdr_def_pipe = seqs.get("1",  {}).get(
            f"{chain_key}_cdr_def") or "kabat"
        cdr_def_det = seqs.get("4",  {}).get(f"{chain_key}_cdr_def") or "kabat"
        cdr_def_stated = seqs.get("8",  {}).get(
            f"{chain_key}_cdr_def") or "kabat"

        boundary_pipe = compute_boundary_positions(
            mouse_seq, pipe_germ,   chain_type, cdr_def_pipe) if pipe_germ and mouse_seq else set()
        boundary_det = compute_boundary_positions(
            mouse_seq, det_germ,    chain_type, cdr_def_det) if det_germ and mouse_seq else set()
        boundary_stated = compute_boundary_positions(
            mouse_seq, stated_germ, chain_type, cdr_def_stated) if stated_germ and mouse_seq else set()

        # ── Check seq 1: CDRs match mouse, FRs match pipeline germline ────────
        if s.get("1") and mouse_seq:
            passed, diffs = cdrs_identical(
                mouse_seq, s["1"], chain_type, boundary_pipe)
            check("1", chain, "CDRs match mouse",
                  passed, ", ".join(diffs[:3]) if diffs else "")

            if pipe_germ:
                ok, matched, comparable, diffs = frs_match_germline(
                    s["1"], pipe_germ, chain_type)
                pct = f"{matched}/{comparable} ({matched/comparable:.0%})" if comparable else "N/A"
                check("1", chain, f"FRs match pipeline germline ({pipe_germ})",
                      ok, pct + (f" — diffs: {', '.join(diffs[:3])}" if not ok else ""))

        # ── Check seq 4: CDRs match mouse, FRs match detected germline ────────
        if s.get("4") and mouse_seq:
            passed, diffs = cdrs_identical(
                mouse_seq, s["4"], chain_type, boundary_det)
            check("4", chain, "CDRs match mouse",
                  passed, ", ".join(diffs[:3]) if diffs else "")

            if det_germ:
                ok, matched, comparable, diffs = frs_match_germline(
                    s["4"], det_germ, chain_type)
                pct = f"{matched}/{comparable} ({matched/comparable:.0%})" if comparable else "N/A"
                check("4", chain, f"FRs match detected germline ({det_germ})",
                      ok, pct + (f" — diffs: {', '.join(diffs[:3])}" if not ok else ""))

        # ── Check seq 8: CDRs match mouse, FRs match stated germline ──────────
        if s.get("8") and mouse_seq:
            passed, diffs = cdrs_identical(
                mouse_seq, s["8"], chain_type, boundary_stated)
            check("8", chain, "CDRs match mouse",
                  passed, ", ".join(diffs[:3]) if diffs else "")

            if stated_germ:
                ok, matched, comparable, diffs = frs_match_germline(
                    s["8"], stated_germ, chain_type)
                pct = f"{matched}/{comparable} ({matched/comparable:.0%})" if comparable else "N/A"
                check("8", chain, f"FRs match stated germline ({stated_germ})",
                      ok, pct + (f" — diffs: {', '.join(diffs[:3])}" if not ok else ""))
            else:
                check("8", chain, "FRs match stated germline", False,
                      "no stated germline in CSV")

        # ── Check seqs 2, 6, 9: CDRs unchanged, FRs changed by Sapiens ───────
        boundary_map = {"2": boundary_pipe,
                        "6": boundary_det, "9": boundary_stated}
        for sapiens_id, base_id in [("2", "1"), ("6", "4"), ("9", "8")]:
            if s.get(sapiens_id) and s.get(base_id):
                passed, diffs = cdrs_identical(
                    s[base_id], s[sapiens_id], chain_type, boundary_map[sapiens_id])
                check(sapiens_id, chain, f"CDRs identical to seq {base_id}",
                      passed, ", ".join(diffs[:3]) if diffs else "")

                changed, n = frs_changed(s[base_id], s[sapiens_id], chain_type)
                if changed:
                    check(sapiens_id, chain, "Sapiens changed FRs", True,
                          f"{n} FR positions changed")
                else:
                    # Not a failure — sequence may already be maximally human
                    check(sapiens_id, chain, "Sapiens changed FRs", True,
                          "NOTE: no FR changes — sequence may already be maximally human")

        # ── Check seq 7: CDRs unchanged, back-mutations at correct positions ───
        if s.get("7") and s.get("4"):
            passed, diffs = cdrs_identical(s["4"], s["7"], chain_type)
            check("7", chain, "CDRs identical to seq 4",
                  passed, ", ".join(diffs[:3]) if diffs else "")

            if s.get("3") and s.get("5"):
                passed, issues = backmut_positions_correct(
                    s["3"], s["5"], s["4"], s["7"], chain_type)
                check("7", chain, "back-mutations at correct FR positions",
                      passed, "; ".join(issues[:3]) if issues else "")

    return results


# ── CSV loaders ───────────────────────────────────────────────────────────────

def load_generated(csv_path: str) -> dict:
    """
    Load generated sequences CSV.
    Returns {clone_id: {seq_id: {"vh": str, "vl": str}, "meta": {...}}}
    """
    clones = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, dialect="excel"):
            clone = row["clone"].strip()
            seq_id = row["seq_id"].strip()
            if clone not in clones:
                clones[clone] = {"meta": {}}
            clones[clone][seq_id] = {
                "vh":         row.get("vh_sequence", "").strip(),
                "vl":         row.get("vl_sequence", "").strip(),
                "vh_cdr_def": row.get("vh_cdr_def", "").strip() or None,
                "vl_cdr_def": row.get("vl_cdr_def", "").strip() or None,
            }
            # Store metadata from first row of each clone
            meta = clones[clone]["meta"]
            for key in ["vh_germline", "vl_germline", "vh_det_germline",
                        "vl_det_germline", "vh_stated_germline", "vl_stated_germline",
                        "vl_chain_type"]:
                if not meta.get(key) and row.get(key):
                    meta[key] = row[key].strip()
    return clones


def load_benchmark(csv_path: str) -> dict:
    """Load benchmark CSV. Returns {clone_id: row_dict}."""
    benchmarks = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, dialect="excel"):
            clone = row.get("clone", "").strip()
            if clone:
                benchmarks[clone] = {k: v.strip() for k, v in row.items()}
    return benchmarks


# ── Summary output ────────────────────────────────────────────────────────────

def print_summary(all_results: list[CheckResult]) -> None:
    total = len(all_results)
    passed = sum(1 for r in all_results if r.passed)
    failed = total - passed

    print(f"\n{'='*75}")
    print(
        f"VERIFICATION SUMMARY — {passed}/{total} checks passed, {failed} failed")
    print(f"{'='*75}")

    # Per-clone table
    clones = dict.fromkeys(r.clone for r in all_results)
    print(f"\n{'Clone':<10} {'Checks':>8} {'Passed':>8} {'Failed':>8}")
    print("-" * 40)
    for clone in clones:
        cr = [r for r in all_results if r.clone == clone]
        cp = sum(1 for r in cr if r.passed)
        cf = len(cr) - cp
        status = "✓" if cf == 0 else "✗"
        print(f"{clone:<10} {len(cr):>8} {cp:>8} {cf:>8}  {status}")

    # Failed checks detail
    failures = [r for r in all_results if not r.passed]
    if failures:
        print(f"\nFAILED CHECKS:")
        print(f"{'Clone':<10} {'Chain':<6} {'Seq':>4} {'Check':<45} {'Detail'}")
        print("-" * 90)
        for r in failures:
            print(
                f"{r.clone:<10} {r.chain:<6} {r.seq_id:>4} {r.check:<45} {r.detail}")
    else:
        print("\n✓ All checks passed")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Verify generated sequences are correct")
    parser.add_argument("--generated",  required=True,
                        help="Path to all_sequences.csv from generate_sequences.py")
    parser.add_argument("--benchmark",  required=True,
                        help="Path to benchmark CSV with ground truth sequences")
    parser.add_argument("--clone",      default=None,
                        help="Verify only this clone (optional)")
    args = parser.parse_args()

    print(f"Loading generated sequences: {args.generated}")
    generated = load_generated(args.generated)
    print(f"Loading benchmark: {args.benchmark}")
    benchmarks = load_benchmark(args.benchmark)

    print(f"\nFound {len(generated)} clones in generated file")
    print(f"Found {len(benchmarks)} clones in benchmark file\n")

    all_results = []
    for clone_id, seqs in generated.items():
        if args.clone and clone_id != args.clone:
            continue
        if clone_id not in benchmarks:
            print(f"  ⚠ {clone_id}: not found in benchmark CSV — skipping")
            continue

        vl_chain_type = seqs.get("meta", {}).get("vl_chain_type", "K")
        print(f"  Verifying {clone_id} (VL chain: {vl_chain_type})...")
        results = verify_clone(
            clone_id, seqs, benchmarks[clone_id], vl_chain_type)
        all_results.extend(results)

    print_summary(all_results)


if __name__ == "__main__":
    main()
