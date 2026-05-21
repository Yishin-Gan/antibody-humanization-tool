# Score Variables Reference

All columns output by `score_sequences.py` into `scores.csv`.

## Key Columns

| Column | Type | Description |
|--------|------|-------------|
| clone | str | Clone name (e.g. 8C11) |
| seq_id | str | Sequence ID (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, etc.) |

## OASis Humanness

| Column | Type | Description |
|--------|------|-------------|
| oasis_vh_identity | float | Overall VH OASis 9-mer prevalence (FR+CDR) |
| oasis_vl_identity | float | Overall VL OASis |
| vh_oasis_fr_identity | float | VH FR-only OASis |
| vl_oasis_fr_identity | float | VL FR-only OASis |
| vh_oasis_cdr_identity | float | VH CDR-only OASis |
| vl_oasis_cdr_identity | float | VL CDR-only OASis |
| vh_oasis_per_position_detail | list[dict] | Per-position OASis data (imgt_pos, nmer, region) |
| vl_oasis_per_position_detail | list[dict] | VL per-position OASis data |

## Germline Identity

| Column | Type | Description |
|--------|------|-------------|
| vh_germline_identity | float | VH FR identity vs canonical germline DB |
| vl_germline_identity | float | VL FR identity vs canonical germline DB |
| vh_germline_fr_detail | list[tuple] | DB germline residues used: [(int_pos, str_aa), ...] in number_sequence coords |
| vl_germline_fr_detail | list[tuple] | VL DB germline residues |

## Post-Sapiens Germline Drift (seq2 and seq6 only)

| Column | Type | Description |
|--------|------|-------------|
| vh_detected_germline_post_sap_seq{N} | str | Post-Sapiens detected VH germline name |
| vh_detected_germline_post_sap_seq{N}_identity | float | Identity vs post-Sapiens germline |
| vh_germline_drifted_seq{N} | str | "True"/"False" — did Sapiens drift to different germline? |
| vh_germline_identity_delta_seq{N} | float | Identity delta (post - pre Sapiens germline) |
| vl_detected_germline_post_sap_seq{N} | str | VL equivalents |
| vl_detected_germline_post_sap_seq{N}_identity | float | |
| vl_germline_drifted_seq{N} | str | |
| vl_germline_identity_delta_seq{N} | float | |

## Structure Confidence (ABodyBuilder2)

| Column | Type | Description |
|--------|------|-------------|
| conf_mean_fv | float | Mean Fv confidence (VH+VL) |
| conf_mean_vh / conf_mean_vl | float | Mean per chain |
| conf_min_vh / conf_min_vl | float | Worst single residue per chain |
| conf_cdr_mean_vh / conf_cdr_mean_vl | float | CDR pooled mean |
| conf_fr_mean_vh / conf_fr_mean_vl | float | FR mean |
| conf_cdr1_mean_vh / conf_cdr2_mean_vh / conf_cdr3_mean_vh | float | Per-CDR VH |
| conf_cdr1_mean_vl / conf_cdr2_mean_vl / conf_cdr3_mean_vl | float | Per-CDR VL |

## CamSol Solubility

| Column | Type | Description |
|--------|------|-------------|
| vh_camsol_score / vl_camsol_score | float | Overall CamSol intrinsic |
| vh_camsol_cdr_score / vl_camsol_cdr_score | float | CDR region CamSol |
| vh_camsol_fr_score / vl_camsol_fr_score | float | FR region CamSol |
| vh_camsol_hotspot_count / vl_camsol_hotspot_count | int | Total hotspots (smoothed < -0.5) |
| vh_camsol_hotspot_cdr_count / vl_camsol_hotspot_cdr_count | int | CDR hotspots |
| vh_camsol_hotspot_fr_count / vl_camsol_hotspot_fr_count | int | FR hotspots |

## Physicochemical

| Column | Type | Description |
|--------|------|-------------|
| fv_pi | float | Isoelectric point of Fv |
| fv_net_charge_ph7 | float | Net charge at pH 7 |

## Liabilities

| Column | Type | Description |
|--------|------|-------------|
| vh_deamidation_cdr_count / vh_deamidation_fr_count | int | NG/NS/NT motifs |
| vh_oxidation_cdr_count / vh_oxidation_fr_count | int | M/W residues |
| vh_isomerization_fr_count | int | DG/DS motifs |
| vl_deamidation_cdr_count / vl_deamidation_fr_count | int | VL equivalents |
| vl_oxidation_cdr_count / vl_oxidation_fr_count | int | |
| vl_isomerization_fr_count | int | |

## Backmutation Detail

| Column | Type | Description |
|--------|------|-------------|
| vh_backmut_detail | list[dict] | Per-position mutable FR data (imgt_pos, query_aa, mouse_aa, status) |
| vl_backmut_detail | list[dict] | VL per-position mutable FR data |

**Status values in backmut_detail:**
- `back_mutated` — query matches mouse (reversion)
- `humanized` — query matches grafted humanized sequence
- `other_substitution` — query differs from both mouse and grafted
- `sapiens_changed_non_mutable` — Sapiens changed a non-mutable position (mouse == grafted)

## Sequence Identity

| Column | Type | Description |
|--------|------|-------------|
| vh_identity_vs_lab_final | float | VH identity vs seq5 |
| vl_identity_vs_lab_final | float | VL identity vs seq5 |

## CDR Properties

| Column | Type | Description |
|--------|------|-------------|
| vh_cdr1_sequence / vh_cdr2_sequence / vh_cdr3_sequence | str | CDR sequences |
| vl_cdr1_sequence / vl_cdr2_sequence / vl_cdr3_sequence | str | VL CDR sequences |
