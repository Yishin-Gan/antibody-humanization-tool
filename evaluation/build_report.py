#!/usr/bin/env python3
"""
build_report.py — Antibody Humanization Analysis Report Builder

Usage:
    python3 build_report.py \
        --scores    outputs/scores.csv \
        --comparison outputs/comparison.csv \
        --sequences outputs/all_sequences.csv \
        --output    outputs/report.html

Generates a self-contained HTML report with five tabs:
  - Features & Metrics  (OASis, CamSol, structure, liabilities)
  - Germline Identity   (per-sequence germline reference table + FR comparison)
  - Pipeline Funnel     (Levels 0–4 from comparison.csv)
  - Sequence Viewer     (seq2 / seq5 / seq6 aligned, CDR highlighted)
  - Mutable Positions   (all FR mutable positions: Vernier + non-Vernier)
"""

import argparse
import ast
import csv
import json
import sys
from pathlib import Path


# ── IMGT Vernier positions ────────────────────────────────────────────────────
VERNIER_VH = {2, 27, 29, 30, 47, 48, 67, 69, 71, 78, 80, 93, 94}
VERNIER_VL = {2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71}

# Which germline column to use for germline_identity per seq_id
GERM_COL_MAP = {
    "0": "pipeline", "0r": "pipeline",
    "1": "pipeline", "2": "pipeline", "2r": "pipeline",
    "3": "detected", "4": "detected", "5": "detected",
    "6": "detected", "6r": "detected", "7": "detected",
    "8": "stated",   "9": "stated",   "9r": "stated",
}

CLONES = ['8C11', '25A1', '2E8', 'Ab21', '10H5', '28E07', '3C3', '1G8', '2B6']


# ── helpers ───────────────────────────────────────────────────────────────────
def fv(v):
    try:
        return round(float(v), 6) if v and v != 'None' else None
    except Exception:
        return None


def iv(v):
    try:
        return int(float(v)) if v and v != 'None' else None
    except Exception:
        return None


def find_regions(seq, cdr1, cdr2, cdr3):
    """Return list of {t, start, end} dicts for FR/CDR regions."""
    if not seq:
        return []
    regions = []
    pos = 0
    for cs, cn in [(cdr1, 'CDR1'), (cdr2, 'CDR2'), (cdr3, 'CDR3')]:
        if not cs:
            continue
        idx = seq.find(cs, pos)
        if idx == -1:
            # Fuzzy fallback: find best match with <= 3 mismatches
            best_i, best_d = -1, len(cs)
            for i in range(pos, len(seq) - len(cs) + 1):
                d = sum(1 for a, b in zip(seq[i:i+len(cs)], cs) if a != b)
                if d < best_d:
                    best_i, best_d = i, d
            if best_d <= 3:
                idx = best_i
            else:
                continue
        if idx > pos:
            regions.append({'t': 'FR', 'start': pos, 'end': idx})
        regions.append({'t': cn, 'start': idx, 'end': idx + len(cs)})
        pos = idx + len(cs)
    if pos < len(seq):
        regions.append({'t': 'FR', 'start': pos, 'end': len(seq)})
    return regions


def extract_per_pos(sc, chain):
    """Extract per-position OASis detail for a chain."""
    key = f'{chain}_oasis_per_position_detail'
    try:
        detail = ast.literal_eval(sc.get(key, '[]'))
        return {str(d['imgt_pos']): {'aa': d['nmer'][0], 'region': d['region']}
                for d in detail}
    except Exception:
        return {}


def extract_germ_db(sc, chain):
    """Extract DB germline FR residues from vh/vl_germline_fr_detail column."""
    key = f'{chain}_germline_fr_detail'
    try:
        items = ast.literal_eval(sc.get(key, '[]'))
        return {str(pos): {'aa': aa, 'region': 'FR'} for pos, aa in items}
    except Exception:
        return {}


def get_germ_ref(clone, sid, sc, sr):
    """Get germline reference names for VH and VL based on seq_id."""
    if sid == '0':
        return sc.get('vh_detected_germline_seq0', ''), sc.get('vl_detected_germline_seq0', '')
    elif sid in ('1', '2', '2r'):
        return sr.get('vh_germline', ''), sr.get('vl_germline', '')
    elif sid in ('3', '4', '5', '6', '6r', '7'):
        return sr.get('vh_det_germline', ''), sr.get('vl_det_germline', '')
    else:
        return sr.get('vh_stated_germline', ''), sr.get('vl_stated_germline', '')


# ── data extraction ───────────────────────────────────────────────────────────
def extract_data(scores_path, comparison_path, sequences_path):
    with open(scores_path) as f:
        score_rows = list(csv.DictReader(f))
    with open(sequences_path) as f:
        seq_rows = list(csv.DictReader(f))
    with open(comparison_path) as f:
        comp_rows = list(csv.DictReader(f))

    scores = {(r['clone'], r['seq_id']): r for r in score_rows}
    seqs   = {(r['clone'], r['seq_id']): r for r in seq_rows}
    comp   = {(r['clone'], r['comparison']): r for r in comp_rows}

    col_map  = {'pipeline': 'vh_germline',     'detected': 'vh_det_germline',
                'stated':   'vh_stated_germline'}
    vcol_map = {'pipeline': 'vl_germline',     'detected': 'vl_det_germline',
                'stated':   'vl_stated_germline'}

    # ── DATA ─────────────────────────────────────────────────────────────
    DATA = {}
    for clone in CLONES:
        DATA[clone] = {}
        for sid in ['0', '2', '5', '6']:
            sr = seqs.get((clone, sid), {})
            sc = scores.get((clone, sid), {})
            vh = sr.get('vh_sequence', '')
            vl = sr.get('vl_sequence', '')
            vh_ref, vl_ref = get_germ_ref(clone, sid, sc, sr)

            DATA[clone][sid] = {
                'vh': vh, 'vl': vl,
                'vh_r': find_regions(vh, sc.get('vh_cdr1_sequence', ''),
                                     sc.get('vh_cdr2_sequence', ''),
                                     sc.get('vh_cdr3_sequence', '')),
                'vl_r': find_regions(vl, sc.get('vl_cdr1_sequence', ''),
                                     sc.get('vl_cdr2_sequence', ''),
                                     sc.get('vl_cdr3_sequence', '')),
                # OASis
                'oa_vh':     fv(sc.get('oasis_vh_identity')),
                'oa_vl':     fv(sc.get('oasis_vl_identity')),
                'oa_fr_vh':  fv(sc.get('vh_oasis_fr_identity')),
                'oa_fr_vl':  fv(sc.get('vl_oasis_fr_identity')),
                'oa_cdr_vh': fv(sc.get('vh_oasis_cdr_identity')),
                'oa_cdr_vl': fv(sc.get('vl_oasis_cdr_identity')),
                # Germline identity
                'g_vh':     fv(sc.get('vh_germline_identity')),
                'g_vl':     fv(sc.get('vl_germline_identity')),
                'g_vh_ref': vh_ref,
                'g_vl_ref': vl_ref,
                # Structure confidence
                'cf_fv':      fv(sc.get('conf_mean_fv')),
                'cf_vh':      fv(sc.get('conf_mean_vh')),
                'cf_vl':      fv(sc.get('conf_mean_vl')),
                'cf_min_vh':  fv(sc.get('conf_min_vh')),
                'cf_min_vl':  fv(sc.get('conf_min_vl')),
                'cf_cdr_vh':  fv(sc.get('conf_cdr_mean_vh')),
                'cf_cdr_vl':  fv(sc.get('conf_cdr_mean_vl')),
                'cf_fr_vh':   fv(sc.get('conf_fr_mean_vh')),
                'cf_fr_vl':   fv(sc.get('conf_fr_mean_vl')),
                'cf_cdr1_vh': fv(sc.get('conf_cdr1_mean_vh')),
                'cf_cdr2_vh': fv(sc.get('conf_cdr2_mean_vh')),
                'cf_cdr3_vh': fv(sc.get('conf_cdr3_mean_vh')),
                'cf_cdr1_vl': fv(sc.get('conf_cdr1_mean_vl')),
                'cf_cdr2_vl': fv(sc.get('conf_cdr2_mean_vl')),
                'cf_cdr3_vl': fv(sc.get('conf_cdr3_mean_vl')),
                # CamSol
                'cs_vh':     fv(sc.get('vh_camsol_score')),
                'cs_vl':     fv(sc.get('vl_camsol_score')),
                'cs_cdr_vh': fv(sc.get('vh_camsol_cdr_score')),
                'cs_fr_vh':  fv(sc.get('vh_camsol_fr_score')),
                'cs_cdr_vl': fv(sc.get('vl_camsol_cdr_score')),
                'cs_fr_vl':  fv(sc.get('vl_camsol_fr_score')),
                'hs_vh':     iv(sc.get('vh_camsol_hotspot_count')),
                'hs_vl':     iv(sc.get('vl_camsol_hotspot_count')),
                'hs_cdr_vh': iv(sc.get('vh_camsol_hotspot_cdr_count')),
                'hs_fr_vh':  iv(sc.get('vh_camsol_hotspot_fr_count')),
                'hs_cdr_vl': iv(sc.get('vl_camsol_hotspot_cdr_count')),
                'hs_fr_vl':  iv(sc.get('vl_camsol_hotspot_fr_count')),
                # Post-Sapiens germline drift (seq2 and seq6 only)
                'ps_germ_vh':  sc.get(f'vh_detected_germline_post_sap_seq{sid}', ''),
                'ps_germ_vl':  sc.get(f'vl_detected_germline_post_sap_seq{sid}', ''),
                'ps_id_vh':    fv(sc.get(f'vh_detected_germline_post_sap_seq{sid}_identity')),
                'ps_id_vl':    fv(sc.get(f'vl_detected_germline_post_sap_seq{sid}_identity')),
                'ps_drift_vh': sc.get(f'vh_germline_drifted_seq{sid}', ''),
                'ps_drift_vl': sc.get(f'vl_germline_drifted_seq{sid}', ''),
                'ps_delta_vh': fv(sc.get(f'vh_germline_identity_delta_seq{sid}')),
                'ps_delta_vl': fv(sc.get(f'vl_germline_identity_delta_seq{sid}')),
                # Identity vs lab final
                'id5_vh': fv(sc.get('vh_identity_vs_lab_final')),
                'id5_vl': fv(sc.get('vl_identity_vs_lab_final')),
                # Physicochemical
                'pi': fv(sc.get('fv_pi')),
                'ch': fv(sc.get('fv_net_charge_ph7')),
                # Liabilities
                'l': {
                    'vh_dc': iv(sc.get('vh_deamidation_cdr_count'))  or 0,
                    'vh_df': iv(sc.get('vh_deamidation_fr_count'))   or 0,
                    'vh_oc': iv(sc.get('vh_oxidation_cdr_count'))    or 0,
                    'vh_of': iv(sc.get('vh_oxidation_fr_count'))     or 0,
                    'vh_if': iv(sc.get('vh_isomerization_fr_count')) or 0,
                    'vl_dc': iv(sc.get('vl_deamidation_cdr_count'))  or 0,
                    'vl_df': iv(sc.get('vl_deamidation_fr_count'))   or 0,
                    'vl_oc': iv(sc.get('vl_oxidation_cdr_count'))    or 0,
                    'vl_of': iv(sc.get('vl_oxidation_fr_count'))     or 0,
                    'vl_if': iv(sc.get('vl_isomerization_fr_count')) or 0,
                },
                # Per-position OASis data
                'vh_pp': extract_per_pos(sc, 'vh'),
                'vl_pp': extract_per_pos(sc, 'vl'),
                # Germline DB FR residues
                'vh_germ_db': extract_germ_db(sc, 'vh'),
                'vl_germ_db': extract_germ_db(sc, 'vl'),
            }

        # Store germline FR sequences from scaffold sequences
        for germ_sid, germ_key in [('1', 'germ_pipe'), ('3', 'germ_lab'),
                                    ('4', 'germ_det'), ('8', 'germ_stated')]:
            sc_g = scores.get((clone, germ_sid), {})
            DATA[clone][germ_key] = {
                'vh': extract_per_pos(sc_g, 'vh'),
                'vl': extract_per_pos(sc_g, 'vl'),
            }

        s5 = scores.get((clone, '5'), {})
        s3 = scores.get((clone, '3'), {})
        DATA[clone]['_m'] = {
            'gvh':       s3.get('oasis_vh_v_germline', ''),
            'gvl':       s3.get('oasis_vl_v_germline', ''),
            's5_oa_vh':  fv(s5.get('oasis_vh_identity')),
            's5_oa_vl':  fv(s5.get('oasis_vl_identity')),
            's5_cs_vh':  fv(s5.get('vh_camsol_score')),
            's5_cf_vh':  fv(s5.get('conf_mean_vh')),
            's5_cf_vl':  fv(s5.get('conf_mean_vl')),
        }

    # ── VERN (Vernier-only positions) ────────────────────────────────────
    VERN = {}
    for clone in CLONES:
        VERN[clone] = {'VH': [], 'VL': []}
        for chain, vset in [('VH', VERNIER_VH), ('VL', VERNIER_VL)]:
            ck = 'vh' if chain == 'VH' else 'vl'
            all_pos = set()
            sd = {}
            mm = {}
            for sid in ['0', '2', '3', '4', '5', '6']:
                r = scores.get((clone, sid), {})
                try:
                    detail = ast.literal_eval(r.get(f'{ck}_backmut_detail', '[]'))
                    sd[sid] = {}
                    for d in detail:
                        pos = d.get('imgt_pos')
                        if pos in vset:
                            all_pos.add(pos)
                            sd[sid][pos] = {
                                'aa': d['query_aa'],
                                't':  d.get('status', ''),
                                'm':  d.get('mouse_aa'),
                            }
                            if d.get('mouse_aa'):
                                mm[pos] = d['mouse_aa']
                except Exception:
                    sd[sid] = {}
            for pos in sorted(all_pos):
                mouse = mm.get(pos, '?')
                row = {'pos': pos, 'mouse': mouse}
                for sid in ['0', '2', '6', '5']:
                    if pos in sd.get(sid, {}):
                        x = sd[sid][pos]
                        row[f's{sid}'] = x['aa']
                        row[f's{sid}_t'] = x['t']
                    else:
                        row[f's{sid}'] = mouse
                        row[f's{sid}_t'] = ('sapiens_kept_mouse' if sid == '0'
                                            else 'not_mutable')
                VERN[clone][chain].append(row)

    # ── MUTABLE (all mutable FR positions, split vernier/non-vernier) ────
    MUTABLE = {}
    for clone in CLONES:
        MUTABLE[clone] = {'VH': {'vern': [], 'non': []},
                          'VL': {'vern': [], 'non': []}}
        for chain, vset in [('VH', VERNIER_VH), ('VL', VERNIER_VL)]:
            ck = 'vh' if chain == 'VH' else 'vl'
            all_pos = set()
            sd = {}
            mm = {}
            for sid in ['0', '1', '2', '3', '4', '5', '6']:
                r = scores.get((clone, sid), {})
                try:
                    detail = ast.literal_eval(r.get(f'{ck}_backmut_detail', '[]'))
                    sd[sid] = {}
                    for d in detail:
                        pos = d.get('imgt_pos')
                        all_pos.add(pos)
                        sd[sid][pos] = {
                            'aa': d['query_aa'],
                            't':  d.get('status', ''),
                            'm':  d.get('mouse_aa'),
                        }
                        if d.get('mouse_aa'):
                            mm[pos] = d['mouse_aa']
                except Exception:
                    sd[sid] = {}
            for pos in sorted(all_pos):
                mouse = mm.get(pos, '?')
                row = {'pos': pos, 'mouse': mouse, 'is_vern': pos in vset}
                for sid in ['2', '5', '6']:
                    if pos in sd.get(sid, {}):
                        x = sd[sid][pos]
                        row[f's{sid}'] = x['aa']
                        row[f's{sid}_t'] = x['t']
                    else:
                        row[f's{sid}'] = mouse
                        row[f's{sid}_t'] = 'not_mutable'
                bucket = 'vern' if pos in vset else 'non'
                MUTABLE[clone][chain][bucket].append(row)

    # ── RPT (Pipeline funnel from comparison.csv) ────────────────────────
    RPT = {}
    for clone in CLONES:
        RPT[clone] = {}
        r0 = comp.get((clone, 'level_0_baseline'), {})
        RPT[clone]['l0'] = {
            'vh_germ': r0.get('vh_lab_germline', ''),
            'vl_germ': r0.get('vl_lab_germline', ''),
            'vh_fr':   fv(r0.get('vh_fr_identity_mouse_vs_seq3')),
            'vl_fr':   fv(r0.get('vl_fr_identity_mouse_vs_seq3')),
            'vh_vm':   iv(r0.get('vh_vernier_mutable_count')),
            'vl_vm':   iv(r0.get('vl_vernier_mutable_count')),
            'vh_cdr3': iv(r0.get('vh_cdr3_length')),
            'vl_cdr3': iv(r0.get('vl_cdr3_length')),
            'bm':      iv(r0.get('total_backmut_count')),
        }
        for level, ctype in [
            ('l1', 'level_1_seq0_scorecard'),
            ('l2', 'level_2_seq2_scorecard'),
            ('l3', 'level_3_seq6_scorecard'),
        ]:
            r = comp.get((clone, ctype), {})
            RPT[clone][level] = {
                'vh_oa': fv(r.get('vh_oasis_identity')),
                'vl_oa': fv(r.get('vl_oasis_identity')),
                'id5':   fv(r.get('vh_identity_vs_seq5')),
                'oa_d':  fv(r.get('vh_oasis_delta_vs_seq5')),
                'frc':   iv(r.get('vh_fr_sapiens_changed_count')),
                'vc':    iv(r.get('vh_vernier_sapiens_changed_count')),
            }
        rc = comp.get((clone, 'level_3_conflict_positions'), {})
        RPT[clone]['l3'].update({
            'vh_oh':    iv(rc.get('vh_conflict_sap_over_humanized')),
            'vh_oc':    iv(rc.get('vh_conflict_sap_over_conservative')),
            'vl_oh':    iv(rc.get('vl_conflict_sap_over_humanized')),
            'vl_oc':    iv(rc.get('vl_conflict_sap_over_conservative')),
            'vh_agree': iv(rc.get('vh_agreement_count')),
            'vl_agree': iv(rc.get('vl_agreement_count')),
        })
        rg = comp.get((clone, 'level_2_germline_selection'), {})
        RPT[clone]['l2'].update({
            'pg':     rg.get('pipeline_vh_germline', ''),
            'lg':     rg.get('lab_vh_germline', ''),
            'fam':    rg.get('vh_germline_family_match', '') == 'True',
            'fri':    fv(rg.get('vh_seq1_vs_seq3_fr_identity')),
            'vd':     iv(rg.get('vh_vernier_mutable_delta')),
            'vl_pg':  rg.get('pipeline_vl_germline', ''),
            'vl_lg':  rg.get('lab_vl_germline', ''),
            'vl_fam': rg.get('vl_germline_family_match', '') == 'True',
            'vl_fri': fv(rg.get('vl_seq1_vs_seq3_fr_identity')),
            'vl_vd':  iv(rg.get('vl_vernier_mutable_delta')),
        })
        r4 = comp.get((clone, 'level_4_validation'), {})
        RPT[clone]['l4'] = {
            'vi': fv(r4.get('vh_seq7_vs_seq5_identity')),
            'li': fv(r4.get('vl_seq7_vs_seq5_identity')),
            'vd': iv(r4.get('vh_discrepancy_count')),
            'ld': iv(r4.get('vl_discrepancy_count')),
        }
        s0v = fv(scores.get((clone, '0'), {}).get('oasis_vh_identity'))
        s2v = fv(scores.get((clone, '2'), {}).get('oasis_vh_identity'))
        s6v = fv(scores.get((clone, '6'), {}).get('oasis_vh_identity'))
        s5v = fv(scores.get((clone, '5'), {}).get('oasis_vh_identity'))
        RPT[clone]['gaps'] = {
            'gA': round(s2v - s0v, 4) if s2v and s0v else None,
            'gB': round(s6v - s2v, 4) if s6v and s2v else None,
            'gC': round(s6v - s5v, 4) if s6v and s5v else None,
        }
        RPT[clone]['s5'] = {
            'vh_oa': s5v,
            'vl_oa': fv(scores.get((clone, '5'), {}).get('oasis_vl_identity')),
            'cs_vh': fv(scores.get((clone, '5'), {}).get('vh_camsol_score')),
        }
        germ = {}
        for sid in ['0', '2', '5', '6']:
            sc = scores.get((clone, sid), {})
            sr = seqs.get((clone, sid), {})
            vh_ref, vl_ref = get_germ_ref(clone, sid, sc, sr)
            germ[sid] = {
                'vh':     fv(sc.get('vh_germline_identity')),
                'vl':     fv(sc.get('vl_germline_identity')),
                'vh_ref': vh_ref,
                'vl_ref': vl_ref,
            }
        RPT[clone]['germ'] = germ

    # ── DRIFT (post-Sapiens germline drift for seq2 and seq6) ────────────
    DRIFT = {}
    for clone in CLONES:
        DRIFT[clone] = {}
        for sid in ['2', '6']:
            sc = scores.get((clone, sid), {})
            DRIFT[clone][sid] = {
                'vh_post':    sc.get(f'vh_detected_germline_post_sap_seq{sid}', ''),
                'vh_post_id': fv(sc.get(f'vh_detected_germline_post_sap_seq{sid}_identity')),
                'vh_drifted': sc.get(f'vh_germline_drifted_seq{sid}', ''),
                'vh_delta':   fv(sc.get(f'vh_germline_identity_delta_seq{sid}')),
                'vl_post':    sc.get(f'vl_detected_germline_post_sap_seq{sid}', ''),
                'vl_post_id': fv(sc.get(f'vl_detected_germline_post_sap_seq{sid}_identity')),
                'vl_drifted': sc.get(f'vl_germline_drifted_seq{sid}', ''),
                'vl_delta':   fv(sc.get(f'vl_germline_identity_delta_seq{sid}')),
            }

    return {'DATA': DATA, 'VERN': VERN, 'MUTABLE': MUTABLE,
            'RPT': RPT, 'DRIFT': DRIFT}


# ── HTML template ─────────────────────────────────────────────────────────────
# Uses __ED_JSON__ as placeholder (replaced via str.replace, no brace escaping)
HTML_TEMPLATE_FILE = Path(__file__).with_name('report_template.html')


def get_html_template() -> str:
    """Load the HTML template from report_template.html next to this script."""
    if HTML_TEMPLATE_FILE.exists():
        return HTML_TEMPLATE_FILE.read_text(encoding='utf-8')
    print(f"ERROR: template not found: {HTML_TEMPLATE_FILE}", file=sys.stderr)
    sys.exit(1)


def build_html(data_json: str) -> str:
    template = get_html_template()
    return template.replace('__ED_JSON__', data_json)


# ── entry point ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='Build antibody humanization report HTML')
    parser.add_argument('--scores',      required=True, help='scores.csv path')
    parser.add_argument('--comparison',  required=True, help='comparison.csv path')
    parser.add_argument('--sequences',   required=True, help='all_sequences.csv path')
    parser.add_argument('--output',      default='report.html', help='output HTML path')
    args = parser.parse_args()

    for p in [args.scores, args.comparison, args.sequences]:
        if not Path(p).exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            sys.exit(1)

    print("Extracting data...")
    data = extract_data(args.scores, args.comparison, args.sequences)
    data_json = json.dumps(data)
    print(f"  DATA: {len(data['DATA'])} clones")
    print(f"  MUTABLE: {sum(len(data['MUTABLE'][c]['VH']['vern']) + len(data['MUTABLE'][c]['VH']['non']) for c in CLONES)} VH rows")
    print(f"  DRIFT: {len(data['DRIFT'])} clones")
    print(f"  JSON: {len(data_json):,} chars")

    print("Building HTML...")
    html = build_html(data_json)

    Path(args.output).write_text(html, encoding='utf-8')
    print(f"  Written: {len(html):,} chars -> {args.output}")


if __name__ == '__main__':
    main()
