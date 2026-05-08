"""
generate_sequences.py

Generates all 9 evaluation sequences for a given antibody clone.

Sequences produced:
  1. pipeline_grafted              — mouse CDRs + pipeline top-1 germline FRs
  2. pipeline_humanized            — seq 1 + Sapiens back-mutations
  3. lab_grafted                   — lab's Hu sequence (ground truth, read from CSV)
  4. detected_germline_grafted     — mouse CDRs + detected lab germline FRs
  5. lab_final                     — lab's final humanized sequence (ground truth, read from CSV)
  6. detected_humanized            — seq 4 + Sapiens back-mutations
  7. detected_direct_backmut       — seq 4 + back-mutations where seq3 != seq5
  8. lab_stated_germline_grafted   — mouse CDRs + lab's stated germline FRs (from database)
  9. lab_stated_germline_humanized — seq 8 + Sapiens back-mutations

Usage (from project root):
    python3 pipeline/generate_sequences.py \\
        --csv data/benchmarks/humanization_benchmark.csv \\
        --output outputs/all_sequences.csv
"""

# isort: skip_file
import sys
sys.path.insert(0, "/workspace/antibody-humanization-tool")  # noqa: E402
from abnumber import Chain
from pipeline.step_b_germline_scoring import rank_germlines, normalize_germline_name, print_normalization_report
from pipeline.step_a_numbering import number_sequence, IMGT_REGIONS
from evaluation.evaluate import detect_germline
from typing import Optional
from dataclasses import dataclass, field
import argparse
import csv
import re


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
    vh_1_pipeline_grafted:              Optional[str] = None
    vh_2_pipeline_humanized:            Optional[str] = None
    vh_3_lab_grafted:                   Optional[str] = None
    vh_4_detected_grafted:              Optional[str] = None
    vh_5_lab_final:                     Optional[str] = None
    vh_6_detected_humanized:            Optional[str] = None
    vh_7_detected_direct_backmut:       Optional[str] = None
    vh_8_stated_germline_grafted:       Optional[str] = None
    vh_9_stated_germline_humanized:     Optional[str] = None

    # VL sequences
    vl_1_pipeline_grafted:              Optional[str] = None
    vl_2_pipeline_humanized:            Optional[str] = None
    vl_3_lab_grafted:                   Optional[str] = None
    vl_4_detected_grafted:              Optional[str] = None
    vl_5_lab_final:                     Optional[str] = None
    vl_6_detected_humanized:            Optional[str] = None
    vl_7_detected_direct_backmut:       Optional[str] = None
    vl_8_stated_germline_grafted:       Optional[str] = None
    vl_9_stated_germline_humanized:     Optional[str] = None

    # Metadata
    vh_pipeline_germline:  Optional[str] = None
    vl_pipeline_germline:  Optional[str] = None
    vh_detected_germline:  Optional[str] = None
    vl_detected_germline:  Optional[str] = None
    vh_stated_germline:    Optional[str] = None  # from CSV vh_germline column
    vl_stated_germline:    Optional[str] = None  # from CSV vl_germline column
    vl_chain_type:         Optional[str] = None
    error:                 Optional[str] = None


# ── Core grafting helper ──────────────────────────────────────────────────────

def graft(mouse_seq: str, germline_name: str,
          scheme: str = "imgt", cdr_definition: str = "kabat") -> Optional[str]:
    """Graft mouse CDRs onto a human germline. Returns grafted sequence or None."""
    try:
        normalized = normalize_germline_name(germline_name)
        mouse_chain = Chain(mouse_seq, scheme=scheme,
                            cdr_definition=cdr_definition)
        grafted = mouse_chain.graft_cdrs_onto_human_germline(
            v_gene=normalized, backmutate_vernier=False)
        return grafted.seq
    except Exception as e:
        print(f"    Grafting failed for {germline_name}: {e}")
        return None


# ── Sapiens humanization helper ───────────────────────────────────────────────

def humanize_sapiens(grafted_seq: str, chain_type: str) -> Optional[str]:
    """
    Apply Sapiens humanization to a grafted sequence.
    Sapiens scores every position independently — including CDRs.
    After humanization, CDR positions are restored from the input sequence
    so that only FR positions are modified.

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
            return grafted_seq

        # Number both sequences to identify CDR positions
        grafted_num = number_sequence(grafted_seq, chain_type=chain_type)
        sapiens_num = number_sequence(sapiens_seq, chain_type=chain_type)

        # Get all CDR positions from the grafted sequence
        cdr_positions = set()
        for pos_key in grafted_num["cdr_residues"]:
            if isinstance(pos_key, int):
                cdr_positions.add(pos_key)

        # Build position → residue map for the Sapiens sequence
        # then overwrite CDR positions with original grafted residues
        grafted_all = {
            **grafted_num["fr_residues"],
            **{k: v for k, v in grafted_num["cdr_residues"].items()
               if isinstance(k, int)}
        }
        sapiens_all = {
            **sapiens_num["fr_residues"],
            **{k: v for k, v in sapiens_num["cdr_residues"].items()
               if isinstance(k, int)}
        }

        # Restore CDR positions from grafted sequence
        for pos in cdr_positions:
            if pos in grafted_all:
                sapiens_all[pos] = grafted_all[pos]

        # Reconstruct sequence in position order
        result_seq = "".join(
            sapiens_all[pos]
            for pos in sorted(sapiens_all.keys())
            if sapiens_all.get(pos)
        )

        # Count FR changes made by Sapiens
        fr_changes = sum(
            1 for pos in grafted_num["fr_residues"]
            if grafted_all.get(pos) != sapiens_all.get(pos)
        )
        print(
            f"    Sapiens: {fr_changes} FR position(s) humanized, CDRs preserved")

        return result_seq

    except ImportError:
        print("    Sapiens not available — install biophi: pip install biophi")
        return None
    except Exception as e:
        print(f"    Sapiens humanization failed: {e}")
        return None


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

        base_fr = dict(numbered_base["fr_residues"])
        base_cdr = {k: v for k, v in numbered_base["cdr_residues"].items()
                    if isinstance(k, int)}

        for pos, aa in backmut_positions.items():
            if pos in base_fr:
                base_fr[pos] = aa

        all_residues = {**base_fr, **base_cdr}
        return "".join(
            all_residues[pos]
            for pos in sorted(all_residues.keys())
            if all_residues.get(pos)
        )
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

        # ── Sequence 1: pipeline grafted ──────────────────────────────────────
        print(f"    Generating seq 1 (pipeline grafted)...")
        result.vh_1_pipeline_grafted = graft(
            mouse_vh, pipe_vh_germ, cdr_definition=cdr_definition)
        result.vl_1_pipeline_grafted = graft(
            mouse_vl, pipe_vl_germ, cdr_definition=cdr_definition)

        # ── Sequence 2: pipeline humanized (Sapiens) ──────────────────────────
        print(f"    Generating seq 2 (pipeline humanized via Sapiens)...")
        if result.vh_1_pipeline_grafted:
            result.vh_2_pipeline_humanized = humanize_sapiens(
                result.vh_1_pipeline_grafted, "H")
        if result.vl_1_pipeline_grafted:
            result.vl_2_pipeline_humanized = humanize_sapiens(
                result.vl_1_pipeline_grafted, vl_chain_type)

        # ── Sequence 3: lab grafted (from CSV — ground truth) ─────────────────
        print(f"    Setting seq 3 (lab grafted — from CSV)...")
        result.vh_3_lab_grafted = lab_hu_vh
        result.vl_3_lab_grafted = lab_hu_vl

        # ── Sequence 4: detected germline grafted ─────────────────────────────
        print(f"    Generating seq 4 (detected germline grafted)...")
        if det_vh_germ:
            result.vh_4_detected_grafted = graft(
                mouse_vh, det_vh_germ, cdr_definition=cdr_definition)
        if det_vl_germ:
            result.vl_4_detected_grafted = graft(
                mouse_vl, det_vl_germ, cdr_definition=cdr_definition)

        # ── Sequence 5: lab final (from CSV — ground truth) ───────────────────
        print(f"    Setting seq 5 (lab final — from CSV)...")
        result.vh_5_lab_final = lab_final_vh
        result.vl_5_lab_final = lab_final_vl

        # ── Sequence 6: detected germline + Sapiens ───────────────────────────
        print(f"    Generating seq 6 (detected germline + Sapiens)...")
        if result.vh_4_detected_grafted:
            result.vh_6_detected_humanized = humanize_sapiens(
                result.vh_4_detected_grafted, "H")
        if result.vl_4_detected_grafted:
            result.vl_6_detected_humanized = humanize_sapiens(
                result.vl_4_detected_grafted, vl_chain_type)

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
        print(f"    Generating seq 8 (lab-stated germline grafted)...")
        if stated_vh_germ:
            result.vh_8_stated_germline_grafted = graft(
                mouse_vh, stated_vh_germ, cdr_definition=cdr_definition)
        else:
            print(f"    Skipping seq 8 VH — no lab-stated germline in CSV")
        if stated_vl_germ:
            result.vl_8_stated_germline_grafted = graft(
                mouse_vl, stated_vl_germ, cdr_definition=cdr_definition)
        else:
            print(f"    Skipping seq 8 VL — no lab-stated germline in CSV")

        # ── Sequence 9: lab-stated germline + Sapiens ─────────────────────────
        print(f"    Generating seq 9 (lab-stated germline + Sapiens)...")
        if result.vh_8_stated_germline_grafted:
            result.vh_9_stated_germline_humanized = humanize_sapiens(
                result.vh_8_stated_germline_grafted, "H")
        if result.vl_8_stated_germline_grafted:
            result.vl_9_stated_germline_humanized = humanize_sapiens(
                result.vl_8_stated_germline_grafted, vl_chain_type)

    except Exception as e:
        result.error = str(e)

    return result


# ── Output helpers ────────────────────────────────────────────────────────────

SEQ_LABELS = {
    "1": "pipeline_grafted",
    "2": "pipeline_humanized",
    "3": "lab_grafted",
    "4": "detected_grafted",
    "5": "lab_final",
    "6": "detected_humanized",
    "7": "detected_direct_backmut",
    "8": "stated_germline_grafted",
    "9": "stated_germline_humanized",
}


def print_sequences(seqs: CloneSequences) -> None:
    print(f"\n  Clone: {seqs.clone_id}")
    if seqs.error:
        print(f"  ERROR: {seqs.error}")
        return
    print(f"  {'Seq':<4} {'Label':<30} {'VH length':>10} {'VL length':>10}")
    print(f"  {'-'*4} {'-'*30} {'-'*10} {'-'*10}")
    for num, label in SEQ_LABELS.items():
        vh = getattr(seqs, f"vh_{num}_{label}", None)
        vl = getattr(seqs, f"vl_{num}_{label}", None)
        print(f"  {num:<4} {label:<30} {str(len(vh) if vh else 'N/A'):>10} "
              f"{str(len(vl) if vl else 'N/A'):>10}")


def export_sequences(all_seqs: list, output_path: str) -> None:
    rows = []
    for seqs in all_seqs:
        if seqs.error:
            continue
        for num, label in SEQ_LABELS.items():
            vh = getattr(seqs, f"vh_{num}_{label}", None)
            vl = getattr(seqs, f"vl_{num}_{label}", None)
            rows.append({
                "clone":           seqs.clone_id,
                "seq_id":          num,
                "seq_label":       label,
                "vh_sequence":     vh or "",
                "vl_sequence":     vl or "",
                "vh_germline":     seqs.vh_pipeline_germline or "",
                "vl_germline":     seqs.vl_pipeline_germline or "",
                "vh_det_germline": seqs.vh_detected_germline or "",
                "vl_det_germline": seqs.vl_detected_germline or "",
                "vh_stated_germline": seqs.vh_stated_germline or "",
                "vl_stated_germline": seqs.vl_stated_germline or "",
                "vl_chain_type":   seqs.vl_chain_type or "",
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
