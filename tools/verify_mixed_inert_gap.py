#!/usr/bin/env python3
"""Exact controls for the mixed-color inert-prime addendum.

This is a finite audit, not a substitute for the all-order proof in
paper/sections/08_mixed_energy_obstruction.tex. Uses the Python standard
library only; no floating-point determinant computations.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def gram_energy_and_colors(rows: list[list[int]]) -> tuple[int, list[int]]:
    n = len(rows)
    assert n > 0 and all(len(row) == n for row in rows)
    q = 0
    for i in range(n):
        for j in range(i):
            counts = [0, 0, 0]
            for x, y in zip(rows[i], rows[j]):
                counts[(x - y) % 3] += 1
            a, b, c = counts
            norm = a*a + b*b + c*c - a*b - b*c - c*a
            difference = (sum(rows[i]) - sum(rows[j])) % 3
            if n % 3 == 0:
                assert norm % 9 == (0 if difference == 0 else 3)
            q += norm
    return q, [sum(row) % 3 for row in rows]


def inert_odd_prime_factor(n: int) -> int | None:
    x = n
    p = 2
    while p*p <= x:
        exp = 0
        while x % p == 0:
            x //= p
            exp += 1
        if p % 3 == 2 and exp % 2:
            return p
        p += 1
    if x > 1 and x % 3 == 2:
        return x
    return None


def v_p(n: int, p: int) -> int:
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def audit_witness() -> dict[str, object]:
    h = json.loads((ROOT / "data/benchmark_normalized.json").read_text())
    n = len(h)
    assert n == 15
    qrow, colors = gram_energy_and_colors(h)
    qcol, ccolors = gram_energy_and_colors([list(c) for c in zip(*h)])
    assert qrow == qcol == 105
    assert sorted(Counter(colors).values()) == [1, 14]
    assert sorted(Counter(ccolors).values()) == [1, 14]
    assert qrow >= 3*n + 15
    return {
        "n": n, "energy": qrow,
        "row_color_counts": dict(Counter(colors)),
        "column_color_counts": dict(Counter(ccolors)),
    }


def audit_residues() -> list[dict[str, int]]:
    results = []
    # Sample norm-obstructed orders, each divisible by 3.
    for n in (15, 33, 45, 51, 69, 87, 99, 105, 123, 135):
        p = inert_odd_prime_factor(n)
        assert p is not None and p >= 5 and n % 2
        k = v_p(n, p)
        assert k % 2
        # Squared-determinant factors for k=0 and both k=1 star cases.
        assert (n*n - 3*n + 3) % p == 3 % p != 0
        assert (n*n - 3*n - 6) % p == (-6) % p != 0
        for tau in (-30, -18, -9, 0, 9, 18, 30):
            last = n**4 - (3*n+6)*n*n + tau*n + 27*(n-3)
            assert last % p == (-81) % p != 0
        assert (n-2)*k % 2 == 1
        assert (n-4)*k % 2 == 1
        results.append({"n": n, "inert_prime": p, "v_p_n": k})
    return results


def audit_order_three() -> dict[str, int]:
    count = 0
    mixed = 0
    for v in product(range(3), repeat=9):
        h = [list(v[3*i:3*i+3]) for i in range(3)]
        q, colors = gram_energy_and_colors(h)
        qt, _ = gram_energy_and_colors([list(c) for c in zip(*h)])
        assert q == qt
        if len(set(colors)) > 1:
            assert q >= 3*(3-1)
            mixed += 1
        count += 1
    assert count == 3**9
    return {"checked": count, "mixed": mixed}


def main() -> None:
    print(json.dumps({
        "status": "PASS",
        "scope": "finite checks only; general result proved in LaTeX",
        "witness": audit_witness(),
        "inert_prime_residues": audit_residues(),
        "order_three_exhaustive": audit_order_three(),
    }, indent=2))


if __name__ == "__main__":
    main()
