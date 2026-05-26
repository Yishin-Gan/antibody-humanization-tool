#!/usr/bin/env python3
"""Launch the Antibody Humanization Advisor web UI on localhost.

Usage:
    python3 run_web.py            # http://127.0.0.1:5000
    python3 run_web.py --port 8080
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from web.app import app


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=5000)
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()
    print(f"Serving Antibody Humanization Advisor on http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug, threaded=False)


if __name__ == "__main__":
    main()
