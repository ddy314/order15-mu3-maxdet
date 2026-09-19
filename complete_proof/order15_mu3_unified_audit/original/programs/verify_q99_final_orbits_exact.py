"""Exact local-eigenspace obstruction for the final Q=99 all-same orbits.

This replaces the last solver-dependent simultaneous-Gram check.

The two surviving signed support orbits are

    A = 4 K2^- + K3^+ + C4^+
    B = 3 K2^- + K3^+ + T6,

where T6 is the double-star with positive center edge and four negative leaf
edges.  Both have characteristic polynomial

    (x-21)^2 (x-18)^4 (x-15)^2 (x-12)^6 (x-9).

Eigenvalue 15 is supported only on C4 in A and only on the two leaf-difference
directions of T6 in B.  Hence the corresponding local block Z of any
intertwiner must solve a small linear commutant/intertwiner equation.  The
remaining columns/rows have repeated-coordinate patterns forced by the same
15-eigenspace.  Exhausting the resulting mu3-valued local solutions shows
that none leaves a compatible residual Gram matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product


@dataclass(frozen=True)
class E:
    a: int = 0
    b: int = 0

    def __add__(self, other):
        if isinstance(other, int):
            other = E(other)
        return E(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return E(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, int):
            return E(self.a * other, self.b * other)
        return E(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    __rmul__ = __mul__

    def conj(self):
        return E(self.a - self.b, -self.b)


ZERO = E()
ONE = E(1)
OMEGA = E(0, 1)
OMEGA2 = E(-1, -1)
ROOTS = (ONE, OMEGA, OMEGA2)
ROOT_SET = set(ROOTS)


def gram(block):
    rows = len(block)
    cols = len(block[0])
    return [
        [
            sum((block[i][k] * block[j][k].conj() for k in range(cols)), ZERO)
            for j in range(rows)
        ]
        for i in range(rows)
    ]


def subtract(left, right):
    return [
        [left[i][j] - right[i][j] for j in range(len(left))]
        for i in range(len(left))
    ]


def c4_gram():
    result = [[E(15 if i == j else 0) for j in range(4)] for i in range(4)]
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 0)):
        result[i][j] = result[j][i] = E(3)
    return result


def t6_gram():
    # centers 0,1; leaves 2,3 attached to 0; leaves 4,5 attached to 1.
    result = [[E(15 if i == j else 0) for j in range(6)] for i in range(6)]
    result[0][1] = result[1][0] = E(3)
    for leaf in (2, 3):
        result[0][leaf] = result[leaf][0] = E(-3)
    for leaf in (4, 5):
        result[1][leaf] = result[leaf][1] = E(-3)
    return result


C4 = c4_gram()
T6 = t6_gram()


def is_c4_residual_shape(residual, outside_columns: int) -> bool:
    """Necessary shape from columns (p,q,p,q)^T."""
    d = E(outside_columns)
    if any(residual[i][i] != d for i in range(4)):
        return False
    s = residual[0][1]
    expected = [
        [d, s, d, s],
        [s.conj(), d, s.conj(), d],
        [d, s, d, s],
        [s.conj(), d, s.conj(), d],
    ]
    return residual == expected


def is_t6_residual_shape(residual, outside_columns: int) -> bool:
    """Necessary shape from equality of each leaf pair in every outside column."""
    d = E(outside_columns)
    if any(residual[i][i] != d for i in range(6)):
        return False
    for j in range(6):
        if residual[2][j] != residual[3][j]:
            return False
        if residual[j][2] != residual[j][3]:
            return False
        if residual[4][j] != residual[5][j]:
            return False
        if residual[j][4] != residual[j][5]:
            return False
    return True


def ab_local_blocks():
    """All mu3-valued 4x6 solutions of C4 Z = Z T6.

    The displayed six-parameter form is the exact solution of the linear
    intertwining equation.
    """
    blocks = []
    for z15, z17, z20, z21, z22, z23 in product(ROOTS, repeat=6):
        entries = (
            -z20-z21, -z22-z23, z15, -z15+z22+z23, z17, -z17+z20+z21,
            -z22-z23, -z20-z21, z21, z20, z23, z22,
            -z20-z21, -z22-z23, -z15+z22+z23, z15, -z17+z20+z21, z17,
            -z22-z23, -z20-z21, z20, z21, z22, z23,
        )
        if all(value in ROOT_SET for value in entries):
            blocks.append([list(entries[6*i:6*i+6]) for i in range(4)])
    return blocks


def aa_local_blocks():
    """All mu3-valued 4x4 solutions of C4 Z = Z C4."""
    blocks = []
    for a10, a11, a12, a13, a14, a15 in product(ROOTS, repeat=6):
        entries = (
            a10, a11, -a10+a13+a15, -a11+a12+a14,
            a14, a15, a12, a13,
            -a10+a13+a15, -a11+a12+a14, a10, a11,
            a12, a13, a14, a15,
        )
        if all(value in ROOT_SET for value in entries):
            blocks.append([list(entries[4*i:4*i+4]) for i in range(4)])
    return blocks


def bb_local_blocks():
    """All mu3-valued 6x6 solutions of T6 Z = Z T6."""
    blocks = []
    for b21, b23, b30, b31, b32, b33, b34, b35 in product(ROOTS, repeat=8):
        entries = (
            -b30+b34+b35, -b31+b32+b33, b31, b31, b30, b30,
            -b31+b32+b33, -b30+b34+b35, b30, b30, b31, b31,
            b31, b30, b21, -b21+b34+b35, b23, -b23+b32+b33,
            b31, b30, -b21+b34+b35, b21, -b23+b32+b33, b23,
            b30, b31, b33, b32, b35, b34,
            b30, b31, b32, b33, b34, b35,
        )
        if all(value in ROOT_SET for value in entries):
            blocks.append([list(entries[6*i:6*i+6]) for i in range(6)])
    return blocks


def certificate():
    ab = ab_local_blocks()
    aa = aa_local_blocks()
    bb = bb_local_blocks()

    ab_residual_survivors = sum(
        is_c4_residual_shape(subtract(C4, gram(block)), 9)
        for block in ab
    )
    aa_residual_survivors = sum(
        is_c4_residual_shape(subtract(C4, gram(block)), 11)
        for block in aa
    )
    bb_residual_survivors = sum(
        is_t6_residual_shape(subtract(T6, gram(block)), 9)
        for block in bb
    )

    assert len(ab) == 144
    assert len(aa) == 225
    assert len(bb) == 729
    assert ab_residual_survivors == 0
    assert aa_residual_survivors == 0
    assert bb_residual_survivors == 0

    return {
        "A/B local mu3 intertwiners": len(ab),
        "A/B residual-compatible": ab_residual_survivors,
        "A/A local mu3 commutants": len(aa),
        "A/A residual-compatible": aa_residual_survivors,
        "B/B local mu3 commutants": len(bb),
        "B/B residual-compatible": bb_residual_survivors,
        "B/A": "excluded by conjugate transpose of A/B",
        "theorem": "All four final Q=99 all-same orbit pairings are impossible.",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(certificate(), indent=2))
