from pathlib import Path
import subprocess, sys, re, json

ROOT = Path(__file__).resolve().parents[1]
ledger = (ROOT / "proofs" / "final_record_level_ledger.md").read_text(encoding="utf-8")
qs = [87,90,96,99,105,108,114,117,123,126,132,135,141,144,150]

missing = []
for q in qs:
    # Flexible line check for markdown rows such as | 90 | CLOSED | ...
    if not re.search(rf"\|\s*{q}\s*\|\s*CLOSED\b", ledger):
        missing.append(q)

if missing:
    raise SystemExit(f"Ledger integrity check failed; missing CLOSED rows: {missing}")

p = subprocess.run(
    [sys.executable, str(ROOT / "programs" / "verify_q99_final_orbits_exact.py")],
    check=True, capture_output=True, text=True
)
obj = json.loads(p.stdout)
assert obj["A/B residual-compatible"] == 0
assert obj["A/A residual-compatible"] == 0
assert obj["B/B residual-compatible"] == 0

print("PASS")
print("All 15 requested Q values are marked CLOSED in the final ledger.")
print("Q=99 standalone exact certificate replays successfully.")
