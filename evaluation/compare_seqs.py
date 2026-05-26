"""
compare_seqs.py

Pairwise sequence comparison and funnel analysis for the antibody humanization pipeline.
Operates on scores.csv and all_sequences.csv to produce the comparison.csv used by
the analysis report.

Output schema: one row per (clone, comparison) — sparse columns across comparison types.
Comparison types:
    level_0_baseline            seq3 vs mouse — humanization difficulty context
    level_1_seq0_scorecard      seq0 objective metrics (Sapiens full autonomy)
    level_2_seq2_scorecard      seq2 objective metrics (pipeline germline + Sapiens)
    level_2_germline_selection  seq1 vs seq3 — germline selection quality
    level_3_seq6_scorecard      seq6 objective metrics (lab germline + Sapiens)
    level_3_conflict_positions  Sapiens vs lab disagreements at Vernier positions
    level_4_validation          seq7 vs seq5 — mechanical back-mutation check

Usage:
    python3 evaluation/compare_seqs.py \\
        --generated outputs/all_sequences.csv \\
        --scores    outputs/scores.csv \\
        --benchmark data/benchmarks/humanization_benchmark.csv \\
        --output    outputs/comparison.csv
"""

# isort: skip_file
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # noqa: E402
from typing import Optional
import argparse
import ast
import csv


# ── Helper: safe float/int conversion ────────────────────────────────────────

def _f(val) -> Optional[float]:
    try:
        return float(val) if val and val != "None" else None
    except (TypeError, ValueError):
        return None


def _i(val) -> Optional[int]:
    try:
        return int(float(val)) if val and val != "None" else None
    except (TypeError, ValueError):
        return None


# ── Sequence identity helpers (lazy import) ───────────────────────────────────

def _fr_identity(seq_a: str, seq_b: str, chain_type: str) -> Optional[float]:
    if not seq_a or not seq_b:
        return None
    try:
        from pipeline.step_a_numbering import number_sequence
        num_a = number_sequence(seq_a, chain_type=chain_type)
        num_b = number_sequence(seq_b, chain_type=chain_type)
        positions = set(num_a["fr_residues"]) & set(num_b["fr_residues"])
        if not positions:
            return None
        matched = sum(1 for p in positions
                      if num_a["fr_residues"][p] == num_b["fr_residues"][p])
        return round(matched / len(positions), 3)
    except Exception:
        return None


def _seq_identity(seq_a: str, seq_b: str, chain_type: str) -> Optional[float]:
    if not seq_a or not seq_b:
        return None
    try:
        from pipeline.step_a_numbering import number_sequence
        num_a = number_sequence(seq_a, chain_type=chain_type)
        num_b = number_sequence(seq_b, chain_type=chain_type)
        all_a = {**num_a["fr_residues"],
                 **{k: v for k, v in num_a["cdr_residues"].items() if isinstance(k, int)}}
        all_b = {**num_b["fr_residues"],
                 **{k: v for k, v in num_b["cdr_residues"].items() if isinstance(k, int)}}
        positions = set(all_a) & set(all_b)
        if not positions:
            return None
        matched = sum(1 for p in positions if all_a[p] == all_b[p])
        return round(matched / len(positions), 3)
    except Exception:
        return None


# ── Level 0: Baseline difficulty ──────────────────────────────────────────────

def level_0_baseline(clone_id, seq3_scores, seq5_scores, vl_chain_type="K"):
    row = {"clone": clone_id, "comparison": "level_0_baseline"}
    row["vh_lab_germline"] = seq3_scores.get("oasis_vh_v_germline", "")
    row["vl_lab_germline"] = seq3_scores.get("oasis_vl_v_germline", "")
    row["vh_fr_identity_mouse_vs_seq3"] = _f(
        seq3_scores.get("vh_identity_vs_mouse"))
    row["vl_fr_identity_mouse_vs_seq3"] = _f(
        seq3_scores.get("vl_identity_vs_mouse"))
    row["vh_vernier_mutable_count"] = _i(
        seq3_scores.get("vh_vernier_mutable_count"))
    row["vl_vernier_mutable_count"] = _i(
        seq3_scores.get("vl_vernier_mutable_count"))
    row["vh_vernier_backmut_count"] = _i(
        seq5_scores.get("vh_vernier_backmut_count"))
    row["vl_vernier_backmut_count"] = _i(
        seq5_scores.get("vl_vernier_backmut_count"))
    row["vh_non_vernier_backmut_count"] = _i(
        seq5_scores.get("vh_non_vernier_backmut_count"))
    row["vl_non_vernier_backmut_count"] = _i(
        seq5_scores.get("vl_non_vernier_backmut_count"))
    row["total_backmut_count"] = sum(filter(None, [
        _i(seq5_scores.get("vh_fr_backmut_count")),
        _i(seq5_scores.get("vl_fr_backmut_count")),
    ]))
    row["vh_cdr3_length"] = _i(seq3_scores.get("vh_cdr3_length"))
    row["vl_cdr3_length"] = _i(seq3_scores.get("vl_cdr3_length"))
    for chain in ("vh", "vl"):
        for liab in ("deamidation", "oxidation", "isomerization", "glycosylation"):
            row[f"{chain}_cdr_{liab}_count"] = _i(
                seq3_scores.get(f"{chain}_{liab}_cdr_count"))
    row["vh_identity_mouse_vs_seq3"] = _f(
        seq3_scores.get("vh_identity_vs_mouse"))
    row["vl_identity_mouse_vs_seq3"] = _f(
        seq3_scores.get("vl_identity_vs_mouse"))
    return row


# ── Levels 1-3: Sequence scorecard ───────────────────────────────────────────

def sequence_scorecard(clone_id, comparison, seq_scores, seq5_scores):
    row = {"clone": clone_id, "comparison": comparison}
    def pull(k): return _f(seq_scores.get(k))
    def pull_i(k): return _i(seq_scores.get(k))

    # Humanness
    row["oasis_identity"] = pull("oasis_identity")
    row["vh_oasis_identity"] = pull("oasis_vh_identity")
    row["vl_oasis_identity"] = pull("oasis_vl_identity")
    row["vh_oasis_cdr_identity"] = pull("vh_oasis_cdr_identity")
    row["vh_oasis_fr_identity"] = pull("vh_oasis_fr_identity")
    row["vl_oasis_cdr_identity"] = pull("vl_oasis_cdr_identity")
    row["vl_oasis_fr_identity"] = pull("vl_oasis_fr_identity")
    row["vh_germline_identity"] = pull("vh_germline_identity")
    row["vl_germline_identity"] = pull("vl_germline_identity")
    row["vh_detected_germline"] = (seq_scores.get("vh_detected_germline_seq0") or
                                   seq_scores.get("oasis_vh_v_germline", ""))
    row["vl_detected_germline"] = (seq_scores.get("vl_detected_germline_seq0") or
                                   seq_scores.get("oasis_vl_v_germline", ""))

    # Vernier zone
    for chain in ("vh", "vl"):
        for metric in ("vernier_mutable_count", "vernier_backmut_count",
                       "vernier_humanized_count", "non_vernier_mutable_count",
                       "non_vernier_backmut_count", "non_vernier_humanized_count"):
            row[f"{chain}_{metric}"] = pull_i(f"{chain}_{metric}")
        row[f"{chain}_backmut_detail"] = seq_scores.get(
            f"{chain}_backmut_detail", str([]))

    # Structure confidence
    for col in ("conf_mean_vh", "conf_mean_vl", "conf_mean_fv",
                "conf_cdr_mean_vh", "conf_fr_mean_vh",
                "conf_cdr_mean_vl", "conf_fr_mean_vl"):
        row[col] = pull(col)

    # CamSol
    for chain in ("vh", "vl"):
        row[f"{chain}_camsol_score"] = pull(f"{chain}_camsol_score")
        row[f"{chain}_camsol_cdr_score"] = pull(f"{chain}_camsol_cdr_score")
        row[f"{chain}_camsol_fr_score"] = pull(f"{chain}_camsol_fr_score")
        row[f"{chain}_camsol_hotspot_count"] = pull_i(
            f"{chain}_camsol_hotspot_count")

    # FR liabilities
    for chain in ("vh", "vl"):
        for liab in ("deamidation", "oxidation", "isomerization",
                     "glycosylation", "asp_pro"):
            row[f"{chain}_{liab}_fr_count"] = pull_i(
                f"{chain}_{liab}_fr_count")

    # Physicochemical
    row["vh_pi"] = pull("vh_pi")
    row["vl_pi"] = pull("vl_pi")
    row["fv_pi"] = pull("fv_pi")
    row["fv_net_charge_ph7"] = pull("fv_net_charge_ph7")

    # Seq0-specific: count positions Sapiens changed at Vernier (stored in detail)
    # For seq0, backmut/humanized counts are 0 (seq0_mode uses humanized_by_sapiens)
    # The mutable_count correctly shows how many positions Sapiens changed total
    row["vh_fr_sapiens_changed_count"] = pull_i("vh_fr_mutable_count")
    row["vh_vernier_sapiens_changed_count"] = pull_i(
        "vh_vernier_mutable_count")
    row["vl_fr_sapiens_changed_count"] = pull_i("vl_fr_mutable_count")
    row["vl_vernier_sapiens_changed_count"] = pull_i(
        "vl_vernier_mutable_count")

    # Reference vs seq5
    row["vh_identity_vs_seq5"] = pull("vh_identity_vs_lab_final")
    row["vl_identity_vs_seq5"] = pull("vl_identity_vs_lab_final")

    vh5 = _f(seq5_scores.get("oasis_vh_identity"))
    vl5 = _f(seq5_scores.get("oasis_vl_identity"))
    row["vh_oasis_delta_vs_seq5"] = (
        round(row["vh_oasis_identity"] - vh5, 4)
        if row["vh_oasis_identity"] is not None and vh5 is not None else None)
    row["vl_oasis_delta_vs_seq5"] = (
        round(row["vl_oasis_identity"] - vl5, 4)
        if row["vl_oasis_identity"] is not None and vl5 is not None else None)

    seq5_vh_bm = _i(seq5_scores.get("vh_vernier_backmut_count"))
    seq5_vl_bm = _i(seq5_scores.get("vl_vernier_backmut_count"))
    row["vh_vernier_backmut_delta_vs_seq5"] = (
        row["vh_vernier_backmut_count"] - seq5_vh_bm
        if row["vh_vernier_backmut_count"] is not None and seq5_vh_bm is not None else None)
    row["vl_vernier_backmut_delta_vs_seq5"] = (
        row["vl_vernier_backmut_count"] - seq5_vl_bm
        if row["vl_vernier_backmut_count"] is not None and seq5_vl_bm is not None else None)

    return row


# ── Level 2: Germline selection sub-comparison ────────────────────────────────

def level_2_germline_selection(clone_id, seq1_scores, seq3_scores,
                               seq1_vh, seq1_vl, seq3_vh, seq3_vl,
                               vl_chain_type="K"):
    row = {"clone": clone_id, "comparison": "level_2_germline_selection"}

    row["pipeline_vh_germline"] = seq1_scores.get("oasis_vh_v_germline", "")
    row["lab_vh_germline"] = seq3_scores.get("oasis_vh_v_germline", "")
    row["pipeline_vl_germline"] = seq1_scores.get("oasis_vl_v_germline", "")
    row["lab_vl_germline"] = seq3_scores.get("oasis_vl_v_germline", "")

    def family(name): return name.split("-")[0] if name else ""
    def gene(name): return name.split("*")[0] if name else ""

    for chain in ("vh", "vl"):
        pipe = row[f"pipeline_{chain}_germline"]
        lab = row[f"lab_{chain}_germline"]
        row[f"{chain}_germline_family_match"] = (
            family(pipe) == family(lab) if pipe and lab else None)
        row[f"{chain}_germline_gene_match"] = (
            gene(pipe) == gene(lab) if pipe and lab else None)

    row["vh_seq1_vs_seq3_fr_identity"] = _fr_identity(seq1_vh, seq3_vh, "H")
    row["vl_seq1_vs_seq3_fr_identity"] = _fr_identity(
        seq1_vl, seq3_vl, vl_chain_type)
    row["vh_seq1_vs_seq3_identity"] = _seq_identity(seq1_vh, seq3_vh, "H")
    row["vl_seq1_vs_seq3_identity"] = _seq_identity(
        seq1_vl, seq3_vl, vl_chain_type)

    for chain in ("vh", "vl"):
        pipe = _i(seq1_scores.get(f"{chain}_vernier_mutable_count"))
        lab = _i(seq3_scores.get(f"{chain}_vernier_mutable_count"))
        row[f"{chain}_pipeline_vernier_mutable"] = pipe
        row[f"{chain}_lab_vernier_mutable"] = lab
        row[f"{chain}_vernier_mutable_delta"] = (
            pipe - lab if pipe is not None and lab is not None else None)

    return row


# ── Level 3: Sapiens conflict positions ───────────────────────────────────────

def level_3_conflict_positions(clone_id, mouse_vh, mouse_vl,
                               seq3_vh, seq3_vl, seq5_vh, seq5_vl,
                               seq6_vh, seq6_vl, vl_chain_type="K"):
    row = {"clone": clone_id, "comparison": "level_3_conflict_positions"}

    try:
        from pipeline.step_a_numbering import number_sequence

        for chain, m, s3, s5, s6, ct in [
            ("vh", mouse_vh, seq3_vh, seq5_vh, seq6_vh, "H"),
            ("vl", mouse_vl, seq3_vl, seq5_vl, seq6_vl, vl_chain_type),
        ]:
            if not all([m, s3, s5, s6]):
                row[f"{chain}_conflict_count"] = None
                row[f"{chain}_conflict_positions"] = str([])
                row[f"{chain}_agreement_count"] = None
                continue

            fr_m = number_sequence(m,  chain_type=ct)["fr_residues"]
            fr_s3 = number_sequence(s3, chain_type=ct)["fr_residues"]
            fr_s5 = number_sequence(s5, chain_type=ct)["fr_residues"]
            fr_s6 = number_sequence(s6, chain_type=ct)["fr_residues"]

            positions = set(fr_m) & set(fr_s3) & set(fr_s5) & set(fr_s6)
            conflicts = []
            agreements = []

            for pos in sorted(positions):
                maa, s3aa, s5aa, s6aa = (
                    fr_m.get(pos), fr_s3.get(pos),
                    fr_s5.get(pos), fr_s6.get(pos))
                if not all([maa, s3aa, s5aa, s6aa]):
                    continue
                if s3aa == maa:
                    continue  # not mutable

                lab_restored = (s5aa == maa)
                sapiens_restored = (s6aa == maa)

                if lab_restored and not sapiens_restored:
                    # Higher risk: lab said keep mouse, Sapiens kept human
                    conflicts.append({
                        "imgt_pos": pos, "mouse_aa": maa,
                        "seq3_aa": s3aa, "seq5_aa": s5aa, "seq6_aa": s6aa,
                        "lab_decision": "back_mutated",
                        "sapiens_decision": "humanized",
                        "direction": "sap_over_humanized",
                    })
                elif not lab_restored and sapiens_restored:
                    # Lower risk: lab kept human, Sapiens restored to mouse
                    conflicts.append({
                        "imgt_pos": pos, "mouse_aa": maa,
                        "seq3_aa": s3aa, "seq5_aa": s5aa, "seq6_aa": s6aa,
                        "lab_decision": "humanized",
                        "sapiens_decision": "back_mutated",
                        "direction": "sap_over_conservative",
                    })
                elif lab_restored and sapiens_restored:
                    agreements.append(
                        {"imgt_pos": pos, "decision": "both_back_mutated"})
                elif not lab_restored and not sapiens_restored:
                    agreements.append(
                        {"imgt_pos": pos, "decision": "both_humanized"})

            # Split conflicts by direction
            sap_over_humanized = [
                c for c in conflicts if c["direction"] == "sap_over_humanized"]
            sap_over_conservative = [
                c for c in conflicts if c["direction"] == "sap_over_conservative"]

            row[f"{chain}_conflict_count"] = len(conflicts)
            row[f"{chain}_conflict_sap_over_humanized"] = len(
                sap_over_humanized)
            row[f"{chain}_conflict_sap_over_conservative"] = len(
                sap_over_conservative)
            row[f"{chain}_conflict_positions"] = str(conflicts)
            row[f"{chain}_agreement_count"] = len(agreements)

    except Exception as e:
        print(f"    Level 3 failed for {clone_id}: {e}")
        for chain in ("vh", "vl"):
            row[f"{chain}_conflict_count"] = None
            row[f"{chain}_conflict_positions"] = str([])
            row[f"{chain}_agreement_count"] = None

    return row


# ── Level 4: Mechanical validation ───────────────────────────────────────────

def level_4_validation(clone_id, seq5_vh, seq5_vl, seq7_vh, seq7_vl,
                       mouse_vh, mouse_vl, seq7_scores, seq5_scores,
                       vl_chain_type="K"):
    row = {"clone": clone_id, "comparison": "level_4_validation"}

    row["vh_seq7_vs_seq5_identity"] = _f(
        seq7_scores.get("vh_identity_vs_lab_final"))
    row["vl_seq7_vs_seq5_identity"] = _f(
        seq7_scores.get("vl_identity_vs_lab_final"))
    row["vh_seq7_vs_seq5_fr_identity"] = _fr_identity(seq7_vh, seq5_vh, "H")
    row["vl_seq7_vs_seq5_fr_identity"] = _fr_identity(
        seq7_vl, seq5_vl, vl_chain_type)

    for chain in ("vh", "vl"):
        row[f"{chain}_seq7_backmut_count"] = _i(
            seq7_scores.get(f"{chain}_fr_backmut_count"))
        row[f"{chain}_seq5_backmut_count"] = _i(
            seq5_scores.get(f"{chain}_fr_backmut_count"))

    # Germline difference between seq4 and seq3 explains most discrepancies:
    # seq7 starts from seq4 (detected germline), seq5 starts from seq3 (lab germline).
    # When germlines differ, residues in seq7 will differ from seq5 at non-back-mutated
    # positions even if the back-mutation logic is correct.
    row["vh_seq4_vs_seq3_fr_identity"] = _f(
        seq7_scores.get("vh_identity_vs_lab_hu"))
    row["vl_seq4_vs_seq3_fr_identity"] = _f(
        seq7_scores.get("vl_identity_vs_lab_hu"))

    try:
        from pipeline.step_a_numbering import number_sequence

        for chain, m, s5, s7, ct in [
            ("vh", mouse_vh, seq5_vh, seq7_vh, "H"),
            ("vl", mouse_vl, seq5_vl, seq7_vl, vl_chain_type),
        ]:
            if not all([m, s5, s7]):
                row[f"{chain}_discrepancy_count"] = None
                row[f"{chain}_discrepancy_positions"] = str([])
                continue

            fr_m = number_sequence(m,  chain_type=ct)["fr_residues"]
            fr_s5 = number_sequence(s5, chain_type=ct)["fr_residues"]
            fr_s7 = number_sequence(s7, chain_type=ct)["fr_residues"]
            positions = set(fr_m) & set(fr_s5) & set(fr_s7)
            discrepancies = []

            for pos in sorted(positions):
                maa, s5aa, s7aa = fr_m.get(pos), fr_s5.get(pos), fr_s7.get(pos)
                if not all([maa, s5aa, s7aa]):
                    continue
                if s5aa != maa:
                    continue  # seq5 didn't back-mutate here
                if s7aa != maa:
                    discrepancies.append({
                        "imgt_pos": pos, "mouse_aa": maa,
                        "seq5_aa": s5aa, "seq7_aa": s7aa,
                    })

            row[f"{chain}_discrepancy_count"] = len(discrepancies)
            row[f"{chain}_discrepancy_positions"] = str(discrepancies)

    except Exception as e:
        print(f"    Level 4 failed for {clone_id}: {e}")
        for chain in ("vh", "vl"):
            row[f"{chain}_discrepancy_count"] = None
            row[f"{chain}_discrepancy_positions"] = str([])

    return row


# ── Master orchestrator ───────────────────────────────────────────────────────

def run_all_comparisons(clone_id, clone_seqs, clone_scores, vl_chain_type="K"):
    results = []

    def gs(sid, chain): return clone_seqs.get(
        sid, {}).get(f"{chain}_sequence", "") or ""

    def sc(sid): return clone_scores.get(sid, {})

    if "3" in clone_scores and "5" in clone_scores:
        results.append(level_0_baseline(
            clone_id, sc("3"), sc("5"), vl_chain_type))

    if "0" in clone_scores and "5" in clone_scores:
        results.append(sequence_scorecard(clone_id, "level_1_seq0_scorecard",
                                          sc("0"), sc("5")))

    if "2" in clone_scores and "5" in clone_scores:
        results.append(sequence_scorecard(clone_id, "level_2_seq2_scorecard",
                                          sc("2"), sc("5")))

    if all(k in clone_scores for k in ("1", "3")) and all(k in clone_seqs for k in ("1", "3")):
        results.append(level_2_germline_selection(
            clone_id, sc("1"), sc("3"),
            gs("1", "vh"), gs("1", "vl"), gs("3", "vh"), gs("3", "vl"),
            vl_chain_type))

    if "6" in clone_scores and "5" in clone_scores:
        results.append(sequence_scorecard(clone_id, "level_3_seq6_scorecard",
                                          sc("6"), sc("5")))

    if all(k in clone_seqs for k in ("3", "5", "6")):
        results.append(level_3_conflict_positions(
            clone_id,
            clone_seqs.get("_mouse_vh", ""), clone_seqs.get("_mouse_vl", ""),
            gs("3", "vh"), gs("3", "vl"), gs("5", "vh"), gs("5", "vl"),
            gs("6", "vh"), gs("6", "vl"), vl_chain_type))

    if all(k in clone_seqs for k in ("5", "7")) and all(k in clone_scores for k in ("7", "5")):
        results.append(level_4_validation(
            clone_id,
            gs("5", "vh"), gs("5", "vl"), gs("7", "vh"), gs("7", "vl"),
            clone_seqs.get("_mouse_vh", ""), clone_seqs.get("_mouse_vl", ""),
            sc("7"), sc("5"), vl_chain_type))

    return results


# ── CSV loaders ───────────────────────────────────────────────────────────────

def load_generated(csv_path):
    result = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            clone = row.get("clone", "").strip()
            seq_id = row.get("seq_id", "").strip()
            if clone not in result:
                result[clone] = {}
            result[clone][seq_id] = {k: v.strip() for k, v in row.items()}
    return result


def load_scores(csv_path):
    result = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            clone = row.get("clone", "").strip()
            seq_id = row.get("seq_id", "").strip()
            if clone not in result:
                result[clone] = {}
            result[clone][seq_id] = {k: v.strip() for k, v in row.items()}
    return result


def load_benchmark(csv_path):
    result = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            clone = row.get("clone", "").strip()
            if clone:
                result[clone] = {k: v.strip() for k, v in row.items()}
    return result


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generated",  required=True)
    parser.add_argument("--scores",     required=True)
    parser.add_argument("--benchmark",  required=True)
    parser.add_argument("--output",     default="outputs/comparison.csv")
    parser.add_argument("--clone",      default=None)
    args = parser.parse_args()

    print(f"Loading sequences:  {args.generated}")
    generated = load_generated(args.generated)
    print(f"Loading scores:     {args.scores}")
    scores = load_scores(args.scores)
    print(f"Loading benchmark:  {args.benchmark}")
    benchmarks = load_benchmark(args.benchmark)
    print(f"Found {len(generated)} clones\n")

    all_results = []

    for clone_id, clone_seqs in generated.items():
        if args.clone and clone_id != args.clone:
            continue
        print(f"  Comparing {clone_id}...")
        bench = benchmarks.get(clone_id, {})
        clone_seqs["_mouse_vh"] = bench.get("mouse_vh", "")
        clone_seqs["_mouse_vl"] = bench.get("mouse_vl", "")
        vl_chain_type = next(
            (row.get("vl_chain_type", "K")
             for row in clone_seqs.values()
             if isinstance(row, dict) and row.get("vl_chain_type")),
            "K")
        clone_scores = scores.get(clone_id, {})
        comparisons = run_all_comparisons(
            clone_id, clone_seqs, clone_scores, vl_chain_type)
        all_results.extend(comparisons)

    if not all_results:
        print("No comparisons generated.")
        return

    all_keys = []
    seen = set()
    for row in all_results:
        for k in row:
            if k not in seen:
                all_keys.append(k)
                seen.add(k)

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys, extrasaction="ignore")
        writer.writeheader()
        for row in all_results:
            writer.writerow(row)

    print(
        f"\nWrote {len(all_results)} comparisons ({len(all_keys)} columns) -> {args.output}")
    types = {}
    for r in all_results:
        t = r.get("comparison", "unknown")
        types[t] = types.get(t, 0) + 1
    for t, n in types.items():
        print(f"  {t}: {n} rows")


if __name__ == "__main__":
    main()
