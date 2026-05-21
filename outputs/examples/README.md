# Reference example outputs

These files are committed snapshots of what the pipeline produces on a
known mouse antibody input. They exist so a new maintainer can see what
"good" output looks like without having to run the pipeline themselves.

| File | What it is |
|---|---|
| `all_sequences.csv` | Every intermediate and final humanized sequence the pipeline generated, with IDs `seq0`–`seq9` |
| `candidates.csv` | Shortlist of grafting candidates considered |
| `grafted_candidates.csv` | The CDR-grafted candidates before back-mutation |
| `comparison.csv` | Side-by-side comparison across modes |
| `evaluation_results.csv` | Per-sequence evaluation metrics |
| `scores.csv` | Per-position humanness / OASis scores (the biggest file) |
| `report.html` | Static HTML report produced by `evaluation/build_report.py` |

**These are reference samples — do not overwrite them.** The live pipeline
writes new runs to `outputs/` at the top level; this folder is left alone.
Regenerate (after intentional pipeline changes) with the scripts in
`pipeline/` and `evaluation/`, then move the new files in here and commit.

For interactive reports produced by the web app, use the in-app Download
buttons (XLSX / FASTA) rather than these files.
