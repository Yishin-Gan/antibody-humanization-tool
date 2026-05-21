# Antibody Humanization Advisor — Handoff & Runbook

This document is for whoever takes over running this tool for the team.
It assumes Docker is installed on the host machine. If you've never
touched it before, that's fine — every command below is copy-paste.

---

## 1. What this is

A web app that takes a mouse antibody (VH + VL sequence pair) and runs a
humanization pipeline against it, producing humanness, structure, solubility,
and back-mutation analysis. Used by the wet-lab team to decide which
humanization strategy to pursue before committing to bench work.

- **Compute lives on the host machine.** Sequences/results never leave.
- **OASis 9-mer database** (~23 GB SQLite file) is the only large data
  asset. It is bind-mounted into the container, not baked into the image.
- **Single-user-at-a-time** by design (one Sapiens/ABodyBuilder2 process).
  Concurrent submissions queue.

## 2. Files you need to know about

| Path | What it is | Touch? |
|---|---|---|
| `Dockerfile` | How the runtime image is built | Only when upgrading deps |
| `docker-compose.yml` | How the container runs (port, mount, restart) | Edit the OASis DB host path if you move the DB |
| `web/` | Flask app source | Code changes |
| `pipeline/`, `evaluation/` | The humanization pipeline (ANARCI, Sapiens, scoring) | Don't touch unless you know what you're doing |
| `data/OASis_9mers_v1.db` | The big DB. Bind-mounted into the container | Move if disk is tight; update `docker-compose.yml` path |

## 3. Day-to-day operations

### Start the service
```bash
cd <project-root>
docker compose up -d
```

It comes up at **http://<host-ip>:5000**. `restart: always` means it
restarts on host reboot once Docker itself starts at boot.

### Stop the service
```bash
docker compose down
```

### Tail logs (when a user reports something broken)
```bash
docker compose logs -f --tail=100 humanization
```

### Restart after a code or dep change
```bash
git pull        # or however you bring updated code in
docker compose up -d --build   # rebuilds the image and restarts
```

### Inspect what's running
```bash
docker compose ps
docker stats humanization
```

## 4. Network access

The container binds `0.0.0.0:5000` inside; compose publishes the same on
the host. There's currently **no authentication** — assume anyone on the
same network can use it. Mitigations, pick one:

- **Trusted LAN only** (default today): rely on the office firewall.
- **Add HTTP basic auth** by putting nginx in front (see Section 7).
- **VPN-only**: bind the published port to a VPN interface instead of
  `0.0.0.0`. In compose change `"5000:5000"` to `"<vpn-ip>:5000:5000"`.

## 5. Troubleshooting

| Symptom | Try this |
|---|---|
| Page won't load | `docker compose ps` — is the container running? If not, `docker compose logs humanization` for the crash reason |
| First user submission hangs ~30 s | Normal: Sapiens model weights load on the first request after a container restart. After that it's cached. |
| "OASis DB not found" in logs | The bind-mount path is wrong. Edit `docker-compose.yml` to point at the actual file, `docker compose up -d` |
| Disk filling up | `docker system prune -af --volumes` clears unused images. The OASis DB and `humanization_jobs` volume are kept. |
| Out of memory | Bump the `memory: 12G` limit in `docker-compose.yml`. Sapiens + ABodyBuilder2 can use ~8–10 GB peak. |
| Container won't start after reboot | Check `docker info` returns successfully (Docker daemon up). On Docker Desktop hosts, the desktop app must be running — see Section 6. |

## 6. Survival across reboots / account changes

The whole point of containerizing was so the tool doesn't die when an
individual's profile is deleted. To realize that:

- **Linux host (recommended):** Docker is a systemd service. Once
  `systemctl enable docker` is run by IT, the container auto-starts on
  every reboot regardless of which human is logged in. Nothing else to do.
- **Windows host with Docker Desktop:** Docker Desktop ties to a Windows
  user account. Have IT install it under a **service / shared account**
  that survives staff turnover (e.g. `lab-shared`), not an individual's
  account. Enable "Start Docker Desktop when you log in" — then log in
  once as that account.
- **Windows host with native Docker Engine (Windows Server, not Desktop):**
  Docker runs as a Windows Service and survives reboots without anyone
  logging in. Preferred for production on Windows.

If the tool ever needs to move to a new machine: `git pull`, copy the
OASis DB, run `docker compose up -d --build`. That's the whole
migration — no Python env reconstruction.

## 7. Hardening for production (when you have time)

These aren't required day one but are the obvious next steps:

1. **Auth via nginx reverse proxy** — `htpasswd`-protected, HTTPS via
   the org's cert, container only listens on localhost. ~15 min of nginx
   config.
2. **Persistent job results** — current state is in-memory; a container
   restart loses past report URLs. The `humanization_jobs` volume is
   declared but the app doesn't yet write to it. Small change in
   `web/app.py` to pickle/restore the `JOBS` dict on shutdown/startup.
3. **Concurrent jobs** — single-worker today (Sapiens etc. are not
   thread-safe). If the team grows, add a job queue (Redis + RQ) and
   parallelism guards.
4. **Backups** — only the OASis DB matters and it's a static download
   from BioPhi. Don't waste backup quota on it; document the download URL.
5. **Monitoring** — `docker compose logs` is fine for a small team.
   For a real SLA, add Prometheus textfile collector or just a cron job
   that `curl`s the home page and alerts on non-200.

## 8. Contacts & external dependencies

| Component | Source / docs | Notes |
|---|---|---|
| ANARCI | https://github.com/oxpig/ANARCI | Antibody numbering. Pure Python + hmmer |
| Sapiens | https://github.com/Merck/Sapiens | Humanization model. Downloads HF weights on first call (pre-warmed in the Dockerfile) |
| BioPhi / OASis | https://github.com/Merck/BioPhi | Humanness scoring. 23 GB DB is the size cost |
| CamSol | (bundled in BioPhi) | Solubility |
| ABodyBuilder2 | https://github.com/oxpig/ImmuneBuilder | Structure confidence (optional; +60 s per submission) |

If any of these license terms change or models go offline, the pipeline
still produces partial results — only the affected metric goes `n/a`.

---

*Last updated: when you took over. Update this file as you make changes.*
