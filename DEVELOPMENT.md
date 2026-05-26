# Developing on the Antibody Humanization Advisor

This document is for someone who wants to **edit the code**, not just run the
tool. If you only want to install and use the app, see [README.md](README.md)
instead. If you operate the deployed app for the team, see
[HANDOFF.md](HANDOFF.md).

You do not need a heavy development environment to contribute to this
project. A text editor, Docker, and git are enough for most changes.

---

## Table of contents

1. [Who this is for and what you'll need](#1-who-this-is-for-and-what-youll-need)
2. [Two dev setups, pick one](#2-two-dev-setups-pick-one)
3. [The dev loop (how editing actually works)](#3-the-dev-loop-how-editing-actually-works)
4. [Project layout](#4-project-layout)
5. [Common-change cookbook](#5-common-change-cookbook)
6. [Project conventions and gotchas](#6-project-conventions-and-gotchas)
7. [Running tests and the sample pipeline](#7-running-tests-and-the-sample-pipeline)
8. [Updating dependencies](#8-updating-dependencies)
9. [External services and large data](#9-external-services-and-large-data)
10. [Branching and pull-request flow](#10-branching-and-pull-request-flow)
11. [When you get stuck](#11-when-you-get-stuck)

---

## 1. Who this is for and what you'll need

You should be comfortable with:

- Using a terminal (running commands, navigating folders).
- Editing files in any text editor.
- Basic `git` (clone, branch, commit, push).
- Reading Python — you don't have to be a Python expert, but you will be reading and changing Python.

You do **not** need to know:

- Docker internals beyond `docker compose up`, `down`, `restart`, `logs`.
- The full antibody humanization pipeline — only the part you're changing.
- Anything about ANARCI, Sapiens, or BioPhi unless you're modifying those steps directly.

**Hardware requirements** are the same as for running the tool — see [README.md §1](README.md#1-what-you-need-before-you-start). The tool runs CPU-only by default; a GPU is only useful if you are doing heavy model R&D, which is outside normal code edits.

---

## 2. Two dev setups, pick one

### Setup A — Edit on host, run in Docker (recommended)

The default and the simplest. The Docker container holds the entire Python
environment; you edit files on your host with any editor and Docker picks them
up. **No Python or scientific packages installed on your host machine.**

1. Install Docker per [README.md §2](README.md#2-install-docker-the-only-prerequisite).
2. Install a text editor. [VSCode](https://code.visualstudio.com/) is a strong choice — its "Dev Containers" extension can open VSCode inside the container so terminal, debugger, and intellisense all see the real runtime. But any editor (Sublime, vim, PyCharm, Notepad++) is fine.
3. Clone the repo:
   ```bash
   git clone https://github.com/Yishin-Gan/antibody-humanization-tool.git
   cd antibody-humanization-tool
   ```
4. Run `./install.sh` once — same install everyone else does. Confirms the image and database are in place.
5. Add bind-mounts for the source code so your edits are visible inside the container. Edit `docker-compose.yml` and add four lines under `volumes:`:
   ```yaml
       volumes:
         - ./data/OASis_9mers_v1.db:/data/OASis_9mers_v1.db:ro
         - humanization_jobs:/app/outputs/jobs
         # ↓ ADD THESE for live source editing
         - ./web:/app/web
         - ./pipeline:/app/pipeline
         - ./evaluation:/app/evaluation
         - ./run_web.py:/app/run_web.py
   ```
6. Restart the container so the bind-mounts take effect:
   ```bash
   docker compose up -d
   ```

That's the entire dev setup. Edit files in your local clone; changes are
visible inside the container instantly.

### Setup B — Host-native Python (only if you want Flask hot-reload)

If you're iterating heavily on Python and want Flask to auto-restart on every
save, install Python deps directly on your host (or in WSL Ubuntu). Slower to
set up but faster per-edit. **Not recommended unless you're doing many
Python-side changes per hour.**

```bash
# In Ubuntu WSL or macOS or Linux:
sudo apt install -y hmmer build-essential python3-pip      # Linux only
pip install -r requirements.txt
pip install flask openpyxl gunicorn
python3 run_web.py --host 0.0.0.0 --debug
```

The `--debug` flag enables Flask's auto-reloader: saving any `.py` file
restarts the server in ~1 second. **Never use this in production** — it is a
development-only server with no security or performance guarantees.

If you go this route, you still need the OASis database at
`data/OASis_9mers_v1.db` (~23 GB) — the app refuses to start without it. The
fastest way to obtain it is to first run `./install.sh` once (which downloads
it), then switch to host-native Python.

---

## 3. The dev loop (how editing actually works)

| You changed... | What you do | How long |
|---|---|---|
| HTML, CSS, or JS in `web/templates/` or `web/static/` | Just refresh the browser. Flask serves static files live. | 0 sec |
| Python in `web/`, `pipeline/`, or `evaluation/` (Setup A) | `docker compose restart humanization` | ~5 sec |
| Python (Setup B, with `--debug`) | Save the file. Flask reloads automatically. | ~1 sec |
| `Dockerfile` or `requirements.txt` | `docker compose up -d --build` (rebuilds the image; first time ~15 min, subsequent rebuilds faster thanks to layer caching) | 1–15 min |
| `docker-compose.yml` | `docker compose up -d` (no `--build` needed) | ~5 sec |

To see what the running container is actually doing while you poke at it:

```bash
docker compose logs -f humanization        # tail logs, Ctrl+C to stop
docker compose exec humanization bash      # open a shell inside the container
```

Cache-busting note: the templates load static assets with versioned query
strings (`?v=11`, etc.). If you change a JS or CSS file and the browser
doesn't pick it up, bump the version number in the `<link>` / `<script>` tag
inside [web/templates/](web/templates/) so browsers refetch.

---

## 4. Project layout

Top-level directories you'll touch most:

| Path | What it does | Touch when... |
|---|---|---|
| [web/](web/) | The Flask web app: routes, report builder, rescoring, downloads, templates, static assets | You're changing the UI or the API |
| [web/app.py](web/app.py) | Flask routes — every endpoint the browser hits | Adding a new endpoint or changing one |
| [web/runner.py](web/runner.py) | Orchestrates the humanization pipeline for one submission | Changing what runs per request, or in what order |
| [web/report_data.py](web/report_data.py) | Converts pipeline output into the JSON the report tab consumes | Adding a metric or column to the report |
| [web/rescore.py](web/rescore.py) | Re-scores a single sequence after the user tweaks residues | Changing what the Tweak tab can recompute |
| [web/download.py](web/download.py) | Builds the XLSX and FASTA downloads | Changing what's in the download files |
| [web/templates/](web/templates/) | Jinja2 HTML — `input.html` for the form, `report.html` for the result page | Changing the page structure |
| [web/static/](web/static/) | CSS + JS for both pages | Changing the look or in-browser behavior |
| [pipeline/](pipeline/) | The humanization pipeline itself — grafting, Sapiens calls, ANARCI numbering | Changing how candidates are generated |
| [evaluation/](evaluation/) | Scoring (humanness, structure, solubility) and the static HTML report builder | Changing how sequences are scored |
| [data/](data/) | The OASis database (gitignored, ~23 GB) and the benchmark CSV | Updating reference data |
| [outputs/](outputs/) | Pipeline run outputs (gitignored at top level) and [outputs/examples/](outputs/examples/) (committed sample outputs) | Don't edit examples by hand — regenerate from a real run |
| [tests/](tests/) | Test scripts (currently gitignored — see Section 7) | You're adding regression tests |

Files at the root you should know about:

- [Dockerfile](Dockerfile) — runtime image definition
- [docker-compose.yml](docker-compose.yml) — runtime deployment
- [run_web.py](run_web.py) — the dev launcher (the production container runs gunicorn instead)
- [requirements.txt](requirements.txt) — pinned Python deps
- [install.sh](install.sh) — one-command installer
- [bundle.sh](bundle.sh) — handoff tarball builder
- [README.md](README.md) / [HANDOFF.md](HANDOFF.md) / DEVELOPMENT.md — docs for users, operators, and developers respectively

---

## 5. Common-change cookbook

A few worked examples for the kinds of changes you'll likely make.

### Add a new metric to the report

1. Add the metric to the data block in [web/report_data.py](web/report_data.py) inside `_scores_block()`. Pull it from the pipeline result and add it to the dict.
2. Add a row in the **Feature Metrics** table in [web/static/report.js](web/static/report.js) (look for `renderMetrics()`).
3. If the metric should appear in the XLSX download, add it to `_metrics_rows()` in [web/download.py](web/download.py).
4. Restart: `docker compose restart humanization`. Submit a sequence, check the Feature Metrics tab.

### Change the input form

Edit [web/templates/input.html](web/templates/input.html). For client-side validation or interactivity, edit [web/static/input.js](web/static/input.js). The page is served by the `/` route in [web/app.py](web/app.py).

### Change report styling

Edit [web/static/report.css](web/static/report.css). Refresh the browser; if the change doesn't show, bump the `?v=` query string on the `<link>` tag in [web/templates/report.html](web/templates/report.html).

### Add a new endpoint

Add it in [web/app.py](web/app.py). Decorate with `@app.get(...)` or `@app.post(...)`. If it returns JSON, return a `dict` — Flask serializes it automatically. Restart the container. Test with `curl http://localhost:5000/your-new-endpoint`.

### Change pipeline behavior

Most pipeline logic is in [pipeline/](pipeline/) and [evaluation/](evaluation/). The web app calls into these via [web/runner.py](web/runner.py) — start there to find what's being called. **Be careful**: many pipeline functions are deterministic and there are committed reference outputs in [outputs/examples/](outputs/examples/) that a successor expects you to match (or update intentionally).

---

## 6. Project conventions and gotchas

These are the rules-of-the-road that are not obvious from reading the code.
Read this section before making non-trivial changes.

### IMGT numbering must come from abnumber, not ANARCI

Both libraries do IMGT numbering, but they disagree on edge-case sequences (especially non-standard CDR1 lengths). The OASis per-position humanness data uses **abnumber's** numbering, so the report MUST also use abnumber for alignment to line up. See `_linear_to_imgt()` and `_scaffold_imgt_residues()` in [web/report_data.py](web/report_data.py). Do not "simplify" this by switching to ANARCI.

### Region boundaries are IMGT, hard-coded

Region buckets use these IMGT positions (inclusive):

```
FR1   1–26
CDR1  27–38
FR2   39–55
CDR2  56–65
FR3   66–104
CDR3  105–117
FR4   118–128
```

These are baked into [web/static/report.js](web/static/report.js) in `buildAlignmentBlock()` and into the download grouping in `web/download.py`. If you change one, change both — they must match.

### Vernier positions differ per chain type

VH Vernier set: `{2, 27, 29, 30, 47, 48, 67, 69, 71, 78, 80, 93, 94}`
VL Vernier set: `{2, 4, 35, 36, 46, 47, 48, 49, 64, 66, 68, 69, 71}`

Some Vernier positions overlap with CDR regions (e.g. VH 27, 29, 30 fall inside IMGT CDR1). This is a numbering-history artifact (Kabat vs IMGT), not a bug — the residue is genuinely both. The report shows Vernier with a coloured bracket box; do not "fix" the overlap.

### Sequence ID convention (seq0–seq9)

The pipeline produces multiple intermediate sequences with stable IDs:

| ID | What it is |
|---|---|
| seq0 | Sapiens applied directly to mouse sequence |
| seq1 | Pipeline grafted (pre-Sapiens) |
| seq2 | Pipeline humanized (post-Sapiens, the main "pipeline" output) |
| seq3 | Lab Hu (the published intermediate, if provided) |
| seq4 | Detected germline grafted |
| seq5 | Lab final (the published humanized antibody, if provided) |
| seq6 | Detected humanized (the "preferred-germline" mode output) |
| seq7 | Direct back-mutation |
| seq8 | Stated germline grafted |
| seq9 | Stated germline humanized |

Mode → seq_id mapping in [web/runner.py](web/runner.py):
- `pipeline` → seq2
- `preferred` → seq6
- `sapiens` → seq0
- `lab` → seq5

### Single-threaded by design

Gunicorn runs with `--workers=1 --threads=1` (see the `CMD` line in [Dockerfile](Dockerfile)). Sapiens, CamSol, and ABodyBuilder2 hold mutable model state and are not safe for concurrent calls. Concurrent user submissions queue — that's intentional. Do not bump worker or thread counts without first making the pipeline thread-safe.

### Dev server default-host is 127.0.0.1

`python3 run_web.py` defaults to localhost-only. For LAN access (e.g. demoing from another machine), you must pass `--host 0.0.0.0`. This trips up everyone the first time.

### JOBS dict is in-memory

The web app stores per-submission job state in an in-memory dict (see `JOBS` in [web/app.py](web/app.py)). A container restart loses past report URLs. The `humanization_jobs` named volume in `docker-compose.yml` is reserved for future persistence work but currently unused. If you add persistence, write the pickled JOBS to that volume on shutdown and reload on startup.

### Static-asset cache busting

Browsers aggressively cache `/static/*.css` and `/static/*.js`. The templates include `?v=N` query strings on every `<link>` and `<script>` for this reason. **Bump the version when you change a static file** or users will see stale assets.

---

## 7. Running tests and the sample pipeline

### Tests

The `tests/` folder is currently gitignored — it has been used as scratch space, not as a maintained test suite. Adding a real test suite is the next thing a new maintainer should do. Suggested minimum:

- A regression test that runs the pipeline on the 4D5 anti-HER2 sequence and asserts key metrics (humanness, region boundaries) match committed reference values.
- A snapshot test that runs the pipeline on the 9-clone benchmark ([data/benchmarks/humanization_benchmark.csv](data/benchmarks/humanization_benchmark.csv)) and diffs the output against [outputs/examples/](outputs/examples/).

When you do add tests, remove `tests/*` from [.gitignore](.gitignore) so they're tracked.

### Smoke test (manual)

Until there's a test suite, the cheap end-to-end check is:

1. `docker compose up -d`
2. Open `http://localhost:5000`
3. Paste the 4D5 sequences from [README.md §7](README.md#7-verifying-the-install-works-end-to-end)
4. Confirm all five report tabs render and the Feature Metrics has no unexpected `n/a` values

### Static (offline) report regeneration

Sometimes useful for diffing against `outputs/examples/report.html`:

```bash
docker compose exec humanization python3 evaluation/build_report.py
```

This regenerates [outputs/report.html](outputs/report.html). Compare against [outputs/examples/report.html](outputs/examples/report.html) to spot regressions.

---

## 8. Updating dependencies

Add or change a line in [requirements.txt](requirements.txt), then rebuild the image:

```bash
docker compose up -d --build
```

The build will pip-install everything fresh. Expect 5–15 minutes depending on
what changed. **Commit the updated `requirements.txt` and explain why in the
commit message** — antibody/structure libraries have known version
incompatibilities, so the reason for the bump matters.

If you add a system package (something installed via `apt-get`), edit the
`RUN apt-get install` line in [Dockerfile](Dockerfile) instead.

---

## 9. External services and large data

The pipeline relies on external models and data — most of which are baked into the image, but the OASis database lives outside.

| Resource | Where it lives | How to refresh |
|---|---|---|
| Sapiens model weights | Pre-warmed into the image at build time (see Dockerfile) — fetched from HuggingFace Hub | Rebuild the image to refetch |
| OASis 9-mer database (~23 GB) | Bind-mounted from host at `data/OASis_9mers_v1.db` | Re-download from https://zenodo.org/records/5164685/files/OASis_9mers_v1.db.gz |
| ANARCI HMM databases | Installed in the image via `hmmer` apt package | Rebuild the image |
| ABodyBuilder2 weights | Downloaded on first call by the ImmuneBuilder package | Cached inside the container's home directory |

If any external source is offline at install time, the pipeline still produces partial results — only the affected metric goes `n/a`. Search the codebase for `"n/a"` to see how each metric handles unavailability.

---

## 10. Branching and pull-request flow

A reasonable default while the team is small:

- `main` is the deployable branch. Don't commit straight to it; always branch.
- Name feature branches descriptively: `feature/add-pdb-export`, `fix/cdr2-bucketing`, `docs/update-handoff`.
- Push branches to GitHub; open a Pull Request against `main`.
- Self-review the diff before merging — until there's a CI suite, you're your own gatekeeper.
- Squash-merge or rebase-merge to keep `main` linear (avoid noisy merge commits).
- Tag releases when you deploy something the team relies on: `git tag v0.2.0 && git push --tags`.

For coordinating with the live deployment, see [HANDOFF.md §3](HANDOFF.md#3-day-to-day-operations) — `git pull && docker compose up -d --build` is the upgrade path.

---

## 11. When you get stuck

In rough order of usefulness:

1. **Read the logs.** `docker compose logs -f humanization` is the first thing to try when anything misbehaves. Most "the app is broken" symptoms have a stack trace in there.
2. **Grep the codebase.** The project is small enough that `grep -r "function_name" web/ pipeline/ evaluation/` will usually find what you need.
3. **Read [HANDOFF.md](HANDOFF.md).** It has a troubleshooting table aimed at operators, but many of the same answers apply to developers.
4. **Use a Python debugger.** `import pdb; pdb.set_trace()` in the code, then `docker compose exec humanization python3 -i` for an interactive session. Or use VSCode's remote debugger via the Dev Containers extension.
5. **Run the static-report builder separately.** It uses the same scoring functions as the web app but with simpler I/O, often easier to debug: `docker compose exec humanization python3 evaluation/build_report.py`.
6. **Check the issue tracker** at https://github.com/Yishin-Gan/antibody-humanization-tool/issues — and open a new issue if you find a real bug. Future-you (or your successor) will thank you for writing it down.

---

*This document evolves as the project does. If you find something missing or
wrong, update it in the same PR as your code change.*
