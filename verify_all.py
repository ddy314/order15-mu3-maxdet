#!/usr/bin/env python3
"""Run the inherited maximality certificate and the September 2026 checks."""
from __future__ import annotations
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
STEPS = [
    ("inherited 38/38 proof", [sys.executable, str(ROOT / "proof/run_final.py")]),
    ("sharp 23-shell reduction", [sys.executable, str(ROOT / "tools/verify_sharp_reduction.py")]),
    ("extremizer structure and automorphisms", [sys.executable, str(ROOT / "tools/verify_structure.py")]),
]


def main() -> int:
    if not __debug__:
        raise RuntimeError("Assertions must remain enabled; do not use python -O.")
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "OPENBLAS_NUM_THREADS": "1"}
    for label, command in STEPS:
        print(f"\n==> {label}", flush=True)
        subprocess.run(command, cwd=ROOT, env=env, check=True)
    print("\nSUCCESS: maximality, sharp reduction, and extremizer structure all verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
