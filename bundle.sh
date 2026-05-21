#!/usr/bin/env bash
# bundle.sh — package the project as a self-contained handoff.
#
# Three bundle sizes, depending on your --mode flag:
#
#   ./bundle.sh small       # ~5 MB   — code + installer only. First run downloads everything.
#   ./bundle.sh medium      # ~3 GB   — code + pre-built Docker image. First run downloads OASis DB.
#   ./bundle.sh full        # ~25 GB  — code + image + OASis DB. Fully offline-installable.
#
# Output: humanization-advisor-<mode>-<date>.tar.gz in the current directory.
set -e

MODE="${1:-small}"
cd "$(dirname "$0")"
ROOT="$(pwd)"
DATE=$(date +%Y%m%d)
OUT="$ROOT/humanization-advisor-${MODE}-${DATE}.tar.gz"
STAGE=$(mktemp -d)

cleanup() { rm -rf "$STAGE"; }
trap cleanup EXIT

case "$MODE" in
  small|medium|full) ;;
  *) echo "Usage: $0 [small|medium|full]"; exit 1;;
esac

echo "→ Staging source files..."
cp -r pipeline evaluation web "$STAGE/"
cp Dockerfile docker-compose.yml install.sh requirements.txt run_web.py "$STAGE/"
cp HANDOFF.md "$STAGE/" 2>/dev/null || true
# README for the recipient
cat > "$STAGE/README.md" <<'EOF'
# Antibody Humanization Advisor — self-installing bundle

## Install (any host with Docker)

```bash
chmod +x install.sh
./install.sh
```

That's it. The installer detects what's bundled and what needs to be
downloaded, then brings up the web app on port 5000.

- For day-to-day operations, see `HANDOFF.md`.
- For network access from other machines, the app binds to 0.0.0.0
  inside the container and the host port is published on 5000.

## What's in this bundle?

EOF

if [ "$MODE" = "medium" ] || [ "$MODE" = "full" ]; then
  echo "→ Building Docker image (15–20 min on first build)..."
  docker build -t humanization-advisor:latest .

  echo "→ Saving image to image.tar.gz (~3 GB, a few minutes)..."
  docker save humanization-advisor:latest | gzip > "$STAGE/image.tar.gz"
  echo "* image.tar.gz — pre-built Docker image (~$(du -h "$STAGE/image.tar.gz" | cut -f1))" >> "$STAGE/README.md"
fi

if [ "$MODE" = "full" ]; then
  if [ ! -f data/OASis_9mers_v1.db ]; then
    echo "ERROR: data/OASis_9mers_v1.db not found — can't make full bundle."
    exit 1
  fi
  echo "→ Copying OASis DB (~23 GB, ~5 min)..."
  mkdir -p "$STAGE/data"
  cp data/OASis_9mers_v1.db "$STAGE/data/"
  echo "* data/OASis_9mers_v1.db — OASis humanness DB (23 GB)" >> "$STAGE/README.md"
fi

case "$MODE" in
  small)
    cat >> "$STAGE/README.md" <<'EOF'
* source code + Dockerfile + installer
* `install.sh` will build the Docker image and download the OASis DB
  (~23 GB) on first run. Allow ~30 min and a fast internet link.
EOF
    ;;
  medium)
    cat >> "$STAGE/README.md" <<'EOF'
* source code + Dockerfile + installer
* `install.sh` will load the pre-built image and download the OASis DB
  (~23 GB) on first run. Allow ~30 min for the DB download.
EOF
    ;;
  full)
    cat >> "$STAGE/README.md" <<'EOF'
* source code + Dockerfile + installer
* Fully offline. `install.sh` will load the image and use the bundled
  DB. First start in ~10 minutes; no internet required.
EOF
    ;;
esac

echo "→ Compressing bundle to $OUT ..."
tar -czf "$OUT" -C "$STAGE" .

SZ=$(du -h "$OUT" | cut -f1)
echo
echo "=================================================="
echo "  Bundle ready: $OUT  ($SZ)"
echo "=================================================="
echo
echo "Hand off to whoever takes over. They:"
echo "  1. Copy this tarball to the target host"
echo "  2. tar xzf $(basename $OUT)"
echo "  3. cd into the extracted folder and run: ./install.sh"
echo
echo "That's the whole install. App will be at http://<host>:5000"
