#!/usr/bin/env bash
# install.sh — one-command installer for the Antibody Humanization Advisor.
#
# Handles three scenarios automatically:
#   • Offline bundle    — finds image.tar.gz / data/OASis_9mers_v1.db beside this script
#   • Partial bundle    — pre-built image but no DB → offers to download the DB
#   • Source-only       — builds image from Dockerfile, downloads DB if asked
#
# Usage:
#     ./install.sh                # interactive (recommended)
#     ./install.sh --yes          # accept all defaults non-interactively
#     OASIS_URL=https://... ./install.sh   # override DB source
set -e

cd "$(dirname "$0")"
ROOT="$(pwd)"
DATA_DIR="$ROOT/data"
DB_FILE="$DATA_DIR/OASis_9mers_v1.db"
IMAGE_TAR="$ROOT/image.tar.gz"
OASIS_URL="${OASIS_URL:-https://zenodo.org/records/5164685/files/OASis_9mers_v1.db.gz}"
ASSUME_YES="${1:-}"

confirm() {
  [ "$ASSUME_YES" = "--yes" ] && return 0
  read -p "$1 [y/N] " yn
  [[ "$yn" =~ ^[Yy] ]]
}

say()  { printf "\033[1;34m▶ %s\033[0m\n" "$1"; }
ok()   { printf "\033[1;32m✓ %s\033[0m\n" "$1"; }
warn() { printf "\033[1;33m⚠ %s\033[0m\n" "$1"; }
err()  { printf "\033[1;31m✗ %s\033[0m\n" "$1" >&2; }

# ── 1. Docker present? ────────────────────────────────────────────────────────
say "Checking Docker..."
if ! command -v docker >/dev/null 2>&1; then
  err "Docker is not installed."
  echo "  Install Docker first, then re-run this script:"
  echo "    • Linux:   curl -fsSL https://get.docker.com | sh"
  echo "    • Windows: install Docker Desktop from https://docs.docker.com/desktop/"
  echo "    • macOS:   install Docker Desktop from https://docs.docker.com/desktop/"
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  err "Docker is installed but not running. Start it and re-run."
  exit 1
fi
DOCKER_COMPOSE="docker compose"
docker compose version >/dev/null 2>&1 || DOCKER_COMPOSE="docker-compose"
ok "Docker is available ($($DOCKER_COMPOSE version | head -1))"

# ── 2. Image: load from bundle if present, otherwise build ────────────────────
if docker image inspect humanization-advisor:latest >/dev/null 2>&1; then
  ok "Docker image already present"
elif [ -f "$IMAGE_TAR" ]; then
  say "Loading pre-built image from image.tar.gz (a few minutes)..."
  gunzip -c "$IMAGE_TAR" | docker load
  ok "Image loaded"
else
  say "Building image from source. This is slow (~15–20 min on first build)."
  say "It downloads ~2 GB of Python wheels + Sapiens model weights."
  if confirm "Continue?"; then
    $DOCKER_COMPOSE build
    ok "Image built"
  else
    err "Aborted by user."
    exit 1
  fi
fi

# ── 3. OASis DB: present? bundled? need download? ─────────────────────────────
if [ -f "$DB_FILE" ]; then
  ok "OASis DB present ($(du -h "$DB_FILE" | cut -f1))"
else
  warn "OASis DB not found at $DB_FILE"
  echo "  This is a ~23 GB file required for humanness scoring (BioPhi/OASis)."
  echo "  Source: $OASIS_URL"
  echo
  if confirm "Download it now (~30 min on a fast link)?"; then
    mkdir -p "$DATA_DIR"
    say "Downloading..."
    if command -v curl >/dev/null 2>&1; then
      curl -L --progress-bar -o "${DB_FILE}.gz" "$OASIS_URL"
    else
      wget --show-progress -O "${DB_FILE}.gz" "$OASIS_URL"
    fi
    say "Decompressing..."
    gunzip "${DB_FILE}.gz"
    ok "OASis DB ready ($(du -h "$DB_FILE" | cut -f1))"
  else
    err "Aborting. Place the file at $DB_FILE manually, then re-run."
    err "Or set OASIS_URL=<url> and re-run to download from a different source."
    exit 1
  fi
fi

# ── 4. Start the service ──────────────────────────────────────────────────────
say "Bringing up the service..."
$DOCKER_COMPOSE up -d
ok "Container started"

# ── 5. Show how to reach it ───────────────────────────────────────────────────
HOST_IP=$(hostname -I 2>/dev/null | awk '{print $1}')
[ -z "$HOST_IP" ] && HOST_IP="<this-host-ip>"
echo
echo "=================================================="
ok "Antibody Humanization Advisor is running."
echo "  Open in a browser:"
echo "    Local:  http://localhost:5000"
echo "    LAN:    http://$HOST_IP:5000"
echo
echo "  Logs:   $DOCKER_COMPOSE logs -f humanization"
echo "  Stop:   $DOCKER_COMPOSE down"
echo "  Start:  $DOCKER_COMPOSE up -d"
echo "=================================================="
