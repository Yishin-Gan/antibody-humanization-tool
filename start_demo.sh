#!/usr/bin/env bash
# Start the Antibody Humanization Advisor for a team demo.
# Binds to 0.0.0.0 so peers on the same network can reach it.
set -e

PORT="${PORT:-5000}"

# Clean any stale server instances from earlier dev sessions
pkill -f "run_web.py" 2>/dev/null || true
sleep 1

cd "$(dirname "$0")"

# Best-guess workstation IP (first non-loopback IPv4)
HOST_IP=$(ip -4 addr show 2>/dev/null | awk '/inet /{print $2}' | cut -d/ -f1 | grep -v '^127' | head -1)
[ -z "$HOST_IP" ] && HOST_IP="<this-host-ip>"

echo "=================================================="
echo "  Antibody Humanization Advisor — demo mode"
echo "--------------------------------------------------"
echo "  Local:    http://127.0.0.1:$PORT"
echo "  LAN:      http://$HOST_IP:$PORT"
echo "  (no auth — keep on a trusted network)"
echo "  Logs:     this terminal"
echo "  Stop:     Ctrl+C"
echo "=================================================="
echo

exec python3 run_web.py --host 0.0.0.0 --port "$PORT"
