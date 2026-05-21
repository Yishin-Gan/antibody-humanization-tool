# Antibody Humanization Advisor

A web application that takes a **mouse antibody** (a pair of VH and VL sequences) and runs a humanization pipeline against it. It produces an interactive report covering humanness, structure, solubility, and back-mutation analysis — the information a wet-lab team needs to choose which humanization strategy to take to the bench.

- **All compute runs locally** on the host machine. Sequences and results never leave the workstation.
- **One web app** on port `5000`, reachable in a browser from any machine on the same network.
- **Containerised**: the whole stack (Python, ANARCI, Sapiens, BioPhi/OASis, CamSol, ABodyBuilder2) ships as a Docker image. You do not need to install Python or any scientific package by hand.
- **Single-user-at-a-time** by design. Concurrent submissions queue.

If you only need to *operate* the tool after someone else has installed it, jump straight to [HANDOFF.md](HANDOFF.md). If you are installing it from scratch, keep reading.

---

## Table of contents

1. [What you need before you start](#1-what-you-need-before-you-start)
2. [Install Docker (the only prerequisite)](#2-install-docker-the-only-prerequisite)
3. [Get the project files](#3-get-the-project-files)
4. [Run the installer](#4-run-the-installer)
5. [Open the web app](#5-open-the-web-app)
6. [Day-to-day commands](#6-day-to-day-commands)
7. [Verifying the install works end-to-end](#7-verifying-the-install-works-end-to-end)
8. [Troubleshooting](#8-troubleshooting)
9. [Re-packaging the tool for someone else](#9-re-packaging-the-tool-for-someone-else)
10. [What the moving parts are](#10-what-the-moving-parts-are)

---

## 1. What you need before you start

**A computer you control.** This will host the web app for everyone else who uses it. It does not need to be powerful — but it does need to stay on. Typical setup is a lab workstation that nobody reboots casually.

**Minimum specs:**

| Resource | Minimum | Recommended |
|---|---|---|
| RAM | 12 GB free | 16 GB+ |
| Disk space | 30 GB free | 50 GB+ |
| CPU | 4 cores | 8 cores |
| Operating system | Linux, Windows 10/11, or macOS | Linux |
| Internet | Required for first install (~25 GB download) | Same |

**Time budget for first install:**

| Step | How long |
|---|---|
| Install Docker (one-time) | 5 – 15 min |
| Download the project | < 1 min |
| Build the Docker image | 15 – 20 min |
| Download the OASis humanness database (23 GB) | 20 – 45 min on a fast link |
| First app startup | ~30 sec |
| **Total** | **~1 hour** |

If your internet is slow or the machine has no internet at all, see [section 9](#9-re-packaging-the-tool-for-someone-else) — the original maintainer can hand you a USB drive with everything pre-bundled so first install runs in under 10 minutes with zero downloads.

---

## 2. Install Docker (the only prerequisite)

Docker is what packages and runs the tool. If `docker --version` prints a version number in your terminal, skip to [section 3](#3-get-the-project-files).

### On Linux (Ubuntu, Debian, Fedora, etc.)

Open a terminal and run:

```bash
curl -fsSL https://get.docker.com | sh
```

This is the official Docker install script. It will ask for your sudo password. When it finishes, run:

```bash
sudo systemctl enable --now docker      # start now and on every reboot
sudo usermod -aG docker $USER           # let your user run docker without sudo
```

Then **log out and log back in** so the group change takes effect. Verify:

```bash
docker --version
docker run --rm hello-world
```

The `hello-world` test prints a friendly success message if Docker is working.

### On Windows 10 or 11

1. Download **Docker Desktop for Windows** from https://docs.docker.com/desktop/install/windows-install/
2. Double-click the installer and accept defaults.
3. When it finishes, the installer asks you to log out of Windows. Do so, then log back in.
4. Start Docker Desktop from the Start menu. The Docker whale icon will appear in your system tray. Wait until it says "Docker is running" (about 30 seconds the first time).
5. Open **PowerShell** (not the old "Command Prompt") and verify:

```powershell
docker --version
docker run --rm hello-world
```

> **Important for long-term deployment on Windows:** Docker Desktop is tied to a specific Windows user account. If the workstation will be used by many people over years, **install Docker Desktop while logged in as a shared/service account** (something like `lab-shared`), not as an individual person's account. If that individual leaves and IT deletes their profile, Docker goes with them. This is exactly the scenario [HANDOFF.md](HANDOFF.md) Section 6 was written for — read it before deploying on Windows.

### On macOS

1. Download **Docker Desktop for Mac** from https://docs.docker.com/desktop/install/mac-install/ — pick the Apple Silicon or Intel build that matches your chip.
2. Drag the Docker app into your Applications folder.
3. Open Docker from Applications. The Docker whale icon will appear in your menu bar. Wait until it says "Docker is running".
4. Open Terminal and verify:

```bash
docker --version
docker run --rm hello-world
```

---

## 3. Get the project files

You have three ways to get the files onto your machine. Pick the one that matches your situation.

### Option A — Clone from GitHub (recommended; needs internet)

If you have `git` installed (try `git --version` in your terminal), clone the repo:

```bash
git clone https://github.com/Yishin-Gan/antibody-humanization-tool.git
cd antibody-humanization-tool
```

If you do not have git, install it: on Linux `sudo apt install git`, on Windows download from https://git-scm.com/download/win, on macOS run `xcode-select --install`.

### Option B — Download a ZIP from GitHub (no git needed)

1. Open https://github.com/Yishin-Gan/antibody-humanization-tool in a browser.
2. Click the green **"Code"** button → **"Download ZIP"**.
3. Unzip it somewhere sensible, e.g. `~/projects/` on Linux/Mac, or `C:\projects\` on Windows.
4. Open a terminal and `cd` into the unzipped folder. It will be called `antibody-humanization-tool-main` — rename it to `antibody-humanization-tool` if you like.

### Option C — Extract a handoff bundle (offline install)

If the previous maintainer gave you a USB drive or tarball named something like `humanization-advisor-full-20260521.tar.gz`:

```bash
tar xzf humanization-advisor-full-20260521.tar.gz -C ~/projects/
cd ~/projects/humanization-advisor-full-20260521   # or whatever folder it extracted into
```

A "full" bundle (~25 GB) contains the code, the pre-built Docker image, and the 23 GB OASis database — first install runs entirely offline. A "medium" bundle (~3 GB) skips the database (you download that one online). A "small" bundle (~5 MB) is code only.

---

## 4. Run the installer

From inside the project folder, run:

```bash
chmod +x install.sh        # Linux/macOS only — makes the script executable
./install.sh
```

On Windows, run it from **Git Bash** (which ships with Git for Windows) or from **WSL** — `install.sh` is a bash script and will not run in PowerShell directly.

### What the installer actually does

It is interactive and prints a step-by-step progress log. In short:

1. **Checks Docker is installed and running.** If not, it tells you what to do and exits.
2. **Decides whether to build the Docker image or load a pre-built one.**
   - If a file called `image.tar.gz` is in the same folder (you came from a bundle), it loads that. Takes a few minutes.
   - Otherwise it builds the image from the `Dockerfile`. Takes 15 – 20 minutes the first time, and downloads ~2 GB of Python wheels plus the Sapiens model weights. The installer will ask **"Continue?"** before starting the long build — type `y` and press Enter.
3. **Checks for the OASis humanness database** at `data/OASis_9mers_v1.db`.
   - If it is already there (you came from a "full" bundle, or you placed it manually), it skips this step.
   - Otherwise it offers to download it from https://zenodo.org/records/5164685/files/OASis_9mers_v1.db.gz — that is the official BioPhi/OASis distribution. The file is **23 GB** and takes 20 – 45 minutes on a fast link. Type `y` and press Enter to accept.
4. **Starts the container** with `docker compose up -d`. The container is configured to restart automatically whenever the host machine reboots, so once you finish this section you should not need to run `install.sh` again.
5. **Prints the URL** that the web app is now listening on.

### Non-interactive mode

If you are scripting this or just want to accept all defaults without being prompted, run:

```bash
./install.sh --yes
```

### Using a different source for the OASis database

If you already have the database file somewhere local, **do not download it again.** Place it at `data/OASis_9mers_v1.db` (~23 GB) before running `install.sh`:

```bash
mkdir -p data
cp /path/to/your/OASis_9mers_v1.db data/
./install.sh
```

If you want to download it from a mirror instead of Zenodo:

```bash
OASIS_URL=https://your-mirror.example.com/OASis_9mers_v1.db.gz ./install.sh
```

---

## 5. Open the web app

When `install.sh` finishes, it prints something like this:

```
==================================================
Antibody Humanization Advisor is running.
  Open in a browser:
    Local:  http://localhost:5000
    LAN:    http://192.168.1.42:5000
==================================================
```

- **From the same machine** that you installed on: open http://localhost:5000 in any browser.
- **From a different machine on the same network**: open the LAN URL printed by the installer. It is reachable to any browser on the office Wi-Fi / LAN.

### What to enter on the input page

The web app expects two sequences:

| Field | What it is | Example length |
|---|---|---|
| **VH (heavy chain variable region)** | The variable region of your mouse antibody heavy chain, in single-letter amino acids | ~115 residues |
| **VL (light chain variable region)** | The variable region of your mouse antibody light chain | ~107 residues |

Paste each sequence into its box and click **Run**. The pipeline takes about 1 – 2 minutes. Do not close the browser tab — when it finishes, you are redirected to the interactive report automatically.

If the first run after a container restart hangs for ~30 seconds before the report appears, that is normal: the Sapiens model weights load on the first request and are then cached for subsequent runs.

---

## 6. Day-to-day commands

These are run from inside the project folder.

| Goal | Command |
|---|---|
| Check whether the container is running | `docker compose ps` |
| Stop the service | `docker compose down` |
| Start the service (after stopping) | `docker compose up -d` |
| Tail the logs (when a user reports something broken) | `docker compose logs -f --tail=100 humanization` |
| See live resource usage | `docker stats humanization` |
| Restart after pulling a code update | `git pull && docker compose up -d --build` |

The `--build` flag in the last command rebuilds the image; only needed when code or dependencies changed. For day-to-day stop/start, omit it.

---

## 7. Verifying the install works end-to-end

After `install.sh` finishes, do a quick smoke test before declaring victory.

1. **Container running?**
   ```bash
   docker compose ps
   ```
   The `humanization` service should show state `Up`.

2. **App responds?**
   ```bash
   curl -I http://localhost:5000
   ```
   You should see `HTTP/1.0 200 OK` (or similar). If you get "connection refused", wait 10 seconds for the container to finish booting and try again.

3. **Submit a known sequence pair.** Open http://localhost:5000 in a browser and paste in this test pair (the classic anti-HER2 antibody, 4D5):

   - **VH:** `EVQLQQSGPELVKPGASLKLSCTASGFNIKDTYIHWVKQRPEQGLEWIGRIYPTNGYTRYDPKFQDKATITADTSSNTAYLQVSRLTSEDTAVYYCSRWGGDGFYAMDYWGQGASVTVSS`
   - **VL:** `DIVMTQSHKFMSTSVGDRVSITCKASQDVNTAVAWYQQKPGHSPKLLIYSASFRYTGVPDRFTGSRSGTDFTFTISSVQAEDLAVYYCQQHYTTPPTFGGGTKLEIK`

   Click **Run**. The browser will redirect to an interactive 5-tab report within 1 – 2 minutes. If you see Summary numbers and a Sequence View tab with coloured residues, the install is healthy.

4. **Reboot survival test (optional but important):** restart the host machine. After it boots, run `docker compose ps` again — the container should be `Up` without you doing anything. This is the `restart: always` setting in `docker-compose.yml` doing its job and is the reason the tool stays available across power outages and OS updates.

---

## 8. Troubleshooting

| Symptom | What to do |
|---|---|
| `install.sh` exits with "Docker is not installed" | Go back to [section 2](#2-install-docker-the-only-prerequisite) and install Docker. |
| `install.sh` exits with "Docker is installed but not running" | Start Docker. On Linux: `sudo systemctl start docker`. On Windows/macOS: open Docker Desktop from your apps menu and wait for the whale icon to say "running". |
| Browser shows "This site can't be reached" or "ERR_CONNECTION_REFUSED" | Run `docker compose ps` — is the container `Up`? If not, run `docker compose logs humanization` to see the crash reason. Most common cause: the OASis DB path in `docker-compose.yml` does not match where the file actually is. |
| Page loads but submissions never finish | Run `docker compose logs -f humanization` while submitting. Common causes: out of memory (bump the `12G` limit in `docker-compose.yml`), or OASis DB not actually present at the path `docker-compose.yml` points to. |
| "Permission denied" running `./install.sh` | Run `chmod +x install.sh` first. |
| OASis DB download keeps failing or timing out | Download it manually from a stable machine: https://zenodo.org/records/5164685/files/OASis_9mers_v1.db.gz, decompress with `gunzip`, copy the resulting `OASis_9mers_v1.db` file to `data/`, then re-run `./install.sh`. |
| Disk is filling up | Run `docker system prune -af` to clear unused images and build cache. The OASis DB and your saved jobs are kept. |
| First user submission after a restart hangs ~30 sec | Normal — the Sapiens model weights are warming up. After that first request they are cached and subsequent runs start instantly. |
| The container won't restart on its own after a reboot | On Linux, ensure `sudo systemctl enable docker` was run during install. On Windows with Docker Desktop, the user account that installed Docker must be logged in for Docker to start. See [HANDOFF.md](HANDOFF.md) Section 6. |

If none of the above match, read [HANDOFF.md](HANDOFF.md) — it has a deeper troubleshooting table aimed at the person running this for a team.

---

## 9. Re-packaging the tool for someone else

If you ever need to hand this tool off to a new team or a new machine, especially one without internet, use `bundle.sh`. It produces a single self-contained tarball that the recipient can run `./install.sh` on without any further downloads.

```bash
./bundle.sh small      # ~5 MB    code + installer only; recipient downloads everything on first run
./bundle.sh medium     # ~3 GB    + pre-built Docker image; recipient downloads OASis DB only
./bundle.sh full       # ~25 GB   + OASis DB included; fully offline-installable
```

Each produces a `humanization-advisor-<size>-<date>.tar.gz` file. Hand that file (USB drive, file share, whatever) to the recipient. They run:

```bash
tar xzf humanization-advisor-full-20260521.tar.gz
cd humanization-advisor-full-20260521
./install.sh
```

`install.sh` automatically detects what is bundled and what still needs to be downloaded. The recipient never needs to know which bundle size you sent — the installer figures it out.

**Use `full` if the recipient has no internet** (e.g. an air-gapped lab network). Use `medium` if they have internet but you want to spare them the 20-minute image build. Use `small` if you just want them to pull from GitHub themselves and run the installer.

---

## 10. What the moving parts are

For curious operators, here is what is inside the container.

| Component | What it does | Source |
|---|---|---|
| **Flask web app** (`web/`) | Serves the input form and interactive report; orchestrates the pipeline | This repo |
| **ANARCI** | IMGT antibody numbering | https://github.com/oxpig/ANARCI |
| **Sapiens** | Sequence-based humanization model | https://github.com/Merck/Sapiens |
| **BioPhi / OASis** | Humanness scoring using the 23 GB 9-mer database | https://github.com/Merck/BioPhi |
| **CamSol** | Solubility prediction (bundled inside BioPhi) | (via BioPhi) |
| **ABodyBuilder2** | Structure-confidence prediction for the paired Fv | https://github.com/oxpig/ImmuneBuilder |

If any of these dependencies become unavailable, the pipeline still produces a partial report — only the affected metric goes `n/a`.

For deeper operational topics — adding authentication, persisting results across container restarts, surviving Windows account deletions, hardening for production — read [HANDOFF.md](HANDOFF.md). That document is written for the long-term maintainer.

---

*Repository:* https://github.com/Yishin-Gan/antibody-humanization-tool

*Issues, questions, or suggestions:* open an issue at https://github.com/Yishin-Gan/antibody-humanization-tool/issues
