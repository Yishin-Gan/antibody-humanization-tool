# CLAUDE.md — Project context for Claude Code

This file is auto-loaded by Claude Code into every session in this repo.
It captures project-specific conventions that aren't obvious from reading
the code. Keep it short — humans should read [README.md](README.md),
[HANDOFF.md](HANDOFF.md), and [DEVELOPMENT.md](DEVELOPMENT.md) for full
documentation.

## What this project is

A Flask web app that takes a mouse antibody (paired VH + VL sequences)
and runs a humanization pipeline (ANARCI / Sapiens / BioPhi-OASis /
CamSol / ABodyBuilder2), producing a 5-tab interactive report. Runs
locally inside Docker; sequences and results never leave the host.

Three audiences:
- Installer/user → [README.md](README.md)
- Operator running it long-term → [HANDOFF.md](HANDOFF.md)
- Developer editing code → [DEVELOPMENT.md](DEVELOPMENT.md)

## Critical conventions (do NOT violate)

### IMGT numbering must come from `abnumber`, not ANARCI

Both libraries do IMGT numbering but disagree on edge cases (especially
non-standard CDR1 lengths). The OASis per-position humanness data uses
**abnumber's** numbering, so the report MUST also use abnumber for
alignment to line up. See `_linear_to_imgt()` and
`_scaffold_imgt_residues()` in `web/report_data.py`. **Never switch
this to ANARCI for "consistency".** It's been the source of a real bug
already.

### Region boundaries (IMGT, hard-coded, must match across files)

```
FR1   1–26
CDR1  27–38
FR2   39–55
CDR2  56–65
FR3   66–104
CDR3  105–117
FR4   118–128
```

These live in `web/static/report.js` (`buildAlignmentBlock`) AND
`web/download.py` (`_split_by_region`). If you change one, change both.

### Vernier-position / CDR overlap is intentional, not a bug

VH Vernier `{2, 27, 29, 30, 47, 48, 67, 69, 71, 78, 80, 93, 94}` and
VL Vernier `{2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71}` overlap
with IMGT CDR regions in several positions (e.g. VH 27, 29, 30 fall in
CDR1). This is a Kabat-vs-IMGT numbering-history artifact, not a coding
error. The UI shows Vernier with a coloured bracket. **Do not "fix" the
overlap.**

### Sequence IDs (seq0–seq9) — stable, used everywhere

| ID | Meaning |
|---|---|
| seq0 | Sapiens applied directly to mouse |
| seq1 | Pipeline grafted (pre-Sapiens) |
| seq2 | Pipeline humanized (post-Sapiens) — main "pipeline" output |
| seq3 | Lab Hu (intermediate from benchmark CSV) |
| seq4 | Detected germline grafted |
| seq5 | Lab final (ground truth from benchmark) |
| seq6 | Detected humanized — "preferred-germline" mode output |
| seq7 | Direct back-mutation |
| seq8 | Stated germline grafted |
| seq9 | Stated germline humanized |

Mode → seq_id mapping in `web/runner.py`:
- `pipeline` → seq2  |  `preferred` → seq6  |  `sapiens` → seq0  |  `lab` → seq5

### Single-worker by design

Gunicorn is configured `--workers=1 --threads=1` in the Dockerfile.
Sapiens, CamSol, and ABodyBuilder2 hold mutable model state and are
**not** thread-safe. Concurrent user submissions queue. Do not bump
worker/thread counts without making the pipeline thread-safe first.

### Static-asset cache busting

Browsers aggressively cache `/static/*.css|*.js`. Templates include
`?v=N` query strings on every `<link>` and `<script>`. Bump the version
when you change a static file, otherwise users see stale assets.

### OASIS_DB_PATH environment variable

The OASis DB path is read from the `OASIS_DB_PATH` env var (set by
`docker-compose.yml` to `/data/OASis_9mers_v1.db`), with a relative
fallback. Do not hardcode absolute paths anywhere — we did a full pass
removing those (commit `1f4b3c0`).

### JOBS dict is in-memory, lost on restart

Submission results live in `web.app.JOBS`, an in-memory dict. Container
restart loses past report URLs. The `humanization_jobs` named volume in
`docker-compose.yml` is reserved for future persistence but currently
unused. If adding persistence, write `JOBS` to that volume.

## Project layout

```
web/                   Flask app (routes, report builder, downloads, tweak)
  app.py               Routes
  runner.py            Per-submission pipeline orchestration
  report_data.py       Pipeline result → JSON for the report
  rescore.py           Re-scoring for the Tweak tab
  download.py          XLSX + FASTA builders
  templates/           input.html, report.html
  static/              CSS + JS (no build step)
pipeline/              Humanization pipeline (ANARCI, Sapiens, grafting)
evaluation/            Scoring (humanness, structure, solubility) + static report
data/                  OASis DB (gitignored, 23 GB) + committed benchmark CSV
outputs/               Live pipeline outputs (gitignored)
outputs/examples/      Committed reference outputs (do NOT overwrite)
Dockerfile             Production container (gunicorn)
docker-compose.yml     Production deployment (restart: always)
install.sh             One-command installer (per README)
bundle.sh              Builds small/medium/full handoff tarballs
run_web.py             Dev-server launcher (NOT production)
```

## Common operations

```bash
# Day-to-day (operator)
docker compose up -d
docker compose logs -f humanization
docker compose down

# Dev (inside the container or with source bind-mounted)
docker compose restart humanization     # pick up Python changes
docker compose up -d --build            # pick up Dockerfile/requirements changes

# Tests / smoke check
docker compose exec humanization python3 evaluation/build_report.py

# Standalone dev server (not for production)
python3 run_web.py --host 0.0.0.0 --port 5000 [--debug]
```

## Things to avoid

- Hardcoded `/workspace/...`, `/home/yishin/...`, or `/home/claude/...` paths. Use `os.path.dirname(os.path.abspath(__file__))` instead. Cruft from earlier got removed in `1f4b3c0`; don't reintroduce it.
- Switching IMGT numbering to ANARCI (see above).
- Bumping gunicorn workers/threads (see above).
- Adding flask/openpyxl/gunicorn pins to the Dockerfile — they belong in `requirements.txt` (consolidated in `7568a98`).
- Bake the 23 GB OASis DB into the Docker image. Always bind-mount.
- Add backwards-compatibility shims for hardcoded paths "just in case". The audit was complete; nothing should reference the old paths.

## Branch and PR convention

- `main` is releasable; don't commit straight to it.
- Feature branches: `feature/<short-desc>`, `fix/<short-desc>`, `docs/<short-desc>`.
- Squash- or rebase-merge to keep `main` linear.

## External dependencies (where things come from)

- ANARCI — https://github.com/oxpig/ANARCI — IMGT numbering, hmmer-based
- Sapiens — https://github.com/Merck/Sapiens — humanization, HuggingFace weights
- BioPhi / OASis — https://github.com/Merck/BioPhi — humanness; 23 GB DB from Zenodo
- ImmuneBuilder / ABodyBuilder2 — https://github.com/oxpig/ImmuneBuilder — structure confidence
- CamSol — bundled in BioPhi — solubility

If any of these is unavailable at runtime, the pipeline still produces
partial results — only the affected metric goes `n/a`. Search the
codebase for `"n/a"` for the per-metric handling.
