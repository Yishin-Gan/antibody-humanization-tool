# Report Data Mapping

Maps the JSON data structures in the report to their CSV sources.

## Data Structures

`build_report.py` extracts 5 top-level structures injected as `const ED = {DATA, VERN, MUTABLE, RPT, DRIFT}`:

| Structure | Description |
|-----------|-------------|
| `DATA[clone][sid]` | Per-sequence metrics for sid in {0, 2, 5, 6} |
| `DATA[clone][germ_*]` | Scaffold germline FR residues from seq1/seq3/seq4/seq8 |
| `VERN[clone][VH/VL]` | Vernier-only mutable positions |
| `MUTABLE[clone][VH/VL][vern/non]` | ALL mutable FR positions split by Vernier membership |
| `RPT[clone]` | Pipeline funnel data (l0/l1/l2/l3/l4/gaps/s5/germ) |
| `DRIFT[clone][sid]` | Post-Sapiens drift for sid in {2, 6} |

## DATA[clone][sid] Fields

| JS Key | scores.csv Column | Type |
|--------|------------------|------|
| oa_vh | oasis_vh_identity | float |
| oa_vl | oasis_vl_identity | float |
| oa_fr_vh | vh_oasis_fr_identity | float |
| oa_fr_vl | vl_oasis_fr_identity | float |
| oa_cdr_vh | vh_oasis_cdr_identity | float |
| oa_cdr_vl | vl_oasis_cdr_identity | float |
| g_vh | vh_germline_identity | float |
| g_vl | vl_germline_identity | float |
| g_vh_ref | (from all_sequences.csv, per GERM_COL_MAP) | str |
| g_vl_ref | (from all_sequences.csv, per GERM_COL_MAP) | str |
| cf_vh | conf_mean_vh | float |
| cf_vl | conf_mean_vl | float |
| cf_min_vh | conf_min_vh | float |
| cf_min_vl | conf_min_vl | float |
| cf_cdr_vh | conf_cdr_mean_vh | float |
| cf_cdr_vl | conf_cdr_mean_vl | float |
| cf_fr_vh | conf_fr_mean_vh | float |
| cf_fr_vl | conf_fr_mean_vl | float |
| cf_cdr1_vh..cf_cdr3_vl | conf_cdr{N}_mean_{chain} | float |
| cs_vh | vh_camsol_score | float |
| cs_vl | vl_camsol_score | float |
| cs_cdr_vh | vh_camsol_cdr_score | float |
| cs_fr_vh | vh_camsol_fr_score | float |
| cs_cdr_vl | vl_camsol_cdr_score | float |
| cs_fr_vl | vl_camsol_fr_score | float |
| hs_vh | vh_camsol_hotspot_count | int |
| hs_vl | vl_camsol_hotspot_count | int |
| hs_cdr_vh | vh_camsol_hotspot_cdr_count | int |
| hs_fr_vh | vh_camsol_hotspot_fr_count | int |
| hs_cdr_vl | vl_camsol_hotspot_cdr_count | int |
| hs_fr_vl | vl_camsol_hotspot_fr_count | int |
| id5_vh | vh_identity_vs_lab_final | float |
| id5_vl | vl_identity_vs_lab_final | float |
| pi | fv_pi | float |
| ch | fv_net_charge_ph7 | float |
| l.vh_dc | vh_deamidation_cdr_count | int |
| l.vh_df | vh_deamidation_fr_count | int |
| l.vh_oc | vh_oxidation_cdr_count | int |
| l.vh_of | vh_oxidation_fr_count | int |
| l.vh_if | vh_isomerization_fr_count | int |
| l.vl_dc..vl_if | (VL equivalents) | int |
| vh_pp | vh_oasis_per_position_detail | dict |
| vl_pp | vl_oasis_per_position_detail | dict |
| vh_germ_db | vh_germline_fr_detail | dict |
| vl_germ_db | vl_germline_fr_detail | dict |
| ps_germ_vh | vh_detected_germline_post_sap_seq{sid} | str |
| ps_drift_vh | vh_germline_drifted_seq{sid} | str |
| ps_delta_vh | vh_germline_identity_delta_seq{sid} | float |

## DATA[clone][germ_*] Scaffold Sources

| Key | Source seq_id | Purpose |
|-----|--------------|---------|
| germ_pipe | seq1 | Pipeline germline scaffold (seq2 compared against this in Panel 1) |
| germ_lab | seq3 | Lab germline scaffold (seq6 compared against this in Panel 1) |
| germ_det | seq4 | Detected germline scaffold |
| germ_stated | seq8 | Stated germline scaffold (seq5 compared against this in Panel 1) |

## MUTABLE[clone][chain] Structure

```
{
  vern: [{pos, mouse, is_vern, s2, s2_t, s5, s5_t, s6, s6_t}, ...],
  non:  [{pos, mouse, is_vern, s2, s2_t, s5, s5_t, s6, s6_t}, ...]
}
```

Status values: back_mutated, humanized, other_substitution, not_mutable, sapiens_changed_non_mutable

## DRIFT[clone][sid] Structure

```
{
  vh_post: str,      // post-Sapiens detected germline name
  vh_post_id: float, // identity vs post-Sapiens germline
  vh_drifted: str,   // "True" or "False"
  vh_delta: float,   // identity delta
  vl_post: str, vl_post_id: float, vl_drifted: str, vl_delta: float
}
```

## Report Tabs

| Tab | Section ID | Render Function | Data Sources |
|-----|-----------|-----------------|-------------|
| Features & Metrics | features | renderFeatures() | DATA |
| Germline Identity | germline | renderGermline() + renderFRComparison() | DATA, DRIFT, germ_pipe/lab/stated |
| Pipeline Funnel | funnel | renderFunnel() | RPT |
| Sequence Viewer | sequences | renderSequences() | DATA |
| Mutable Positions | vernier | renderVernier() + renderMutableNonVernier() | MUTABLE, VERN, DATA |
