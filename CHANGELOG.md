# Antibody Humanization Pipeline — Changelog

## Branch: feature/longest-cdr-strategy

### Sequence Generation
- Implemented longest-CDR strategy for CDR definition selection
  - For each clone, tries IMGT, Kabat, and Chothia
  - Selects whichever produces the longest total CDR length
  - Mirrors the lab's approach to minimize CDR truncation during humanization
  - Selected definition stored in output CSV (vh_cdr_def, vl_cdr_def columns)

### Bug Fixes
- **CDR3 insertion handling**: Fixed loss of CDR3 insertion-coded positions
  (e.g. `(111,'A')`, `(112,'A')`) during sequence reconstruction
  - Root cause: `isinstance(k, int)` filter dropped tuple-keyed positions
  - Fix: string-level CDR restoration using IMGT position→string index mapping
  - Affected: `humanize_sapiens()` and `apply_direct_backmutations()`

- **Sapiens CDR restoration**: Fixed residue swapping between integer pos112
  and tuple (112,'A') after Sapiens humanization
  - Root cause: Renumbering Sapiens output caused ANARCI to reassign insertion
    positions differently than in the original grafted sequence
  - Fix: Identify CDR string indices from grafted sequence, restore at string
    level without renumbering the Sapiens output

- **Gene-level germline normalization**: Fixed `normalize_germline_name()` for
  gene-level names without allele suffix (e.g. `IGKV1D-7-1`)
  - Root cause: `_GERMLINE_NAME_MAP` keyed by allele-level names only
  - Fix: Try `name*01` if exact gene-level match fails

- **Strategy 4 sequence extraction**: Fixed abnumber chain object iteration
  in normalization map builder
  - Root cause: `ab_chain.seq` and `str(ab_chain)` fail for abnumber Chain
    objects; all similarity scores were 0
  - Fix: `"".join(aa for pos, aa in ab_chain)` for correct sequence extraction

- **Wrong chain type in graft()**: Fixed auto-detection returning `"H"` for
  kappa VL sequences when `chain_type=None`
  - Fix: All `graft()` calls now pass `chain_type` explicitly

- **Absolute imports**: Fixed relative imports in `step_b_germline_scoring.py`
  (`from step_a_numbering` → `from pipeline.step_a_numbering`)

### Verification
- Dynamic boundary computation per clone/germline (replaces hardcoded sets)
  - Computes positions where selected CDR definition and IMGT disagree
  - Results cached to avoid redundant computation
- Added ANARCI renumbering exclusions for gap-adjacent positions
  - VH: {62, 63, 64, 65, 110, 111, 112, 113}
  - VL: {55, 56, 104, 105}
  - These positions shift after Sapiens changes FR sequence

### Verification Results
- **283/284 checks pass**
- 1 known failure: Ab21 VL seq8 — IGKV1-NL1 not in ANARCI database (unfixable)

---

## Branch: feature/sequence-generation

### Sequence Generation
- Fixed baseline using Kabat CDR definition throughout
- Generates 9 evaluation sequences + 3 raw Sapiens outputs per clone:
  - seq 1: pipeline_grafted
  - seq 2: pipeline_humanized (CDR-restored)
  - seq 2r: pipeline_humanized_raw (Sapiens raw)
  - seq 3: lab_grafted (from CSV)
  - seq 4: detected_grafted
  - seq 5: lab_final (from CSV)
  - seq 6: detected_humanized (CDR-restored)
  - seq 6r: detected_humanized_raw
  - seq 7: detected_direct_backmut
  - seq 8: stated_germline_grafted
  - seq 9: stated_germline_humanized (CDR-restored)
  - seq 9r: stated_germline_humanized_raw

### Verification Results
- **277/278 checks pass**
- 1 known failure: Ab21 VL seq8 — IGKV1-NL1 not in ANARCI database (unfixable)

---

## Known Limitations (both branches)

### Database Gaps
| Clone | Chain | Germline   | Issue                            | Affected Seqs |
|-------|-------|------------|----------------------------------|---------------|
| Ab21  | VL    | IGKV1-NL1  | Not in ANARCI database           | seq 8, 9, 9r  |
| 25A1  | VL    | IGKV4-59   | Not in abnumber database         | seq 8, 9, 9r  |

### VL Length Difference vs Lab Ground Truth
- Pipeline-generated VL sequences are 1-2 residues shorter than lab sequences
- Cause: lab sequences contain a J gene tail residue beyond IMGT position 128
- ANARCI does not number this position → FR/CDR comparisons unaffected
- Affected clones: 8C11, 25A1, 2E8, 3C3, 1G8, 2B6, 28E07
- Not affected: Ab21, 10H5

---

## Evaluation Script (evaluate.py)

### Germline Detection Investigation (see germline_detection_investigation.md)
- Root cause of detection mismatches: allele-level inconsistency between
  `rank_germlines()` and `get_germline_fr_by_region()`
  - `rank_germlines` scored against best allele (e.g. IGHV4-30-2*07)
  - `get_germline_fr_by_region` returned first allele found (*01)
  - Fix: `prefer_allele` parameter + Hu-based detection rankings stored separately
- Tie detection: shows warning only when detected and lab germlines share same score
- FR4 excluded from all comparisons (J gene, absent from V gene database)