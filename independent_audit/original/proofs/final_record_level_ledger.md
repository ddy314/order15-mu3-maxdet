# Order-15 mu3 maximal determinant — record-level proof ledger

Benchmark:

\[
D_0=|\det H_0|^2=277868041444786176=2^{22}3^{20}19.
\]

The repository's published benchmark matrix verifies this value exactly.

## Global reductions

1. Let \(G=HH^*\) and
   \[
   Q=\sum_{i<j}|G_{ij}|^2.
   \]
   Exact trace/variance optimization shows every strict counterexample must have \(Q\le168\); \(Q\ge171\) is below \(D_0\).

2. Row colors are exponent sums modulo 3. For a Gram entry with color difference \(\delta\),
   \[
   |G_{ij}|^2/3\equiv \delta^2\pmod3.
   \]
   Hence same-color norms are \(0\pmod9\), different-color norms are \(3\pmod9\). Since the three color-class sizes sum to 15, the number of cross-color pairs is \(0\) or \(2\pmod3\). Therefore
   \[
   Q\equiv0\text{ or }6\pmod9,
   \]
   and every \(Q\equiv3\pmod9\) shell is empty.

3. Schur/color-partition reduction excludes all strict-counterexample partitions except
   \[
   (15,0,0),\quad(14,1,0),\quad(13,2,0),\quad(13,1,1).
   \]

4. The exact ternary code bound \(B_3(15,10)=12\) is used repeatedly. In particular a same-color class of 13 or more rows cannot be internally orthogonal.

## Shell ledger

### Q == 0 (mod 9)

| Q | record-level status | main mechanism |
|---:|---|---|
| 0, 9, 18 | CLOSED | support-independence / \(B_3(15,10)=12\) |
| 27 | CLOSED | exact all-same low-energy component enumeration |
| 36 | CLOSED | exact all-same low-energy component enumeration |
| 45 | CLOSED | exact support enumeration + Eisenstein norm sieve |
| 54 | CLOSED | kernel-compatible support enumeration + K4 residual PSD obstruction |
| 63 | CLOSED | kernel-compatible enumeration + 5-vertex residual PSD obstruction |
| 72 | CLOSED | no-isolate arithmetic obstruction + kernel residual obstruction |
| 81 | CLOSED | exact support reduction; final K3+6K2 branch has a structural symmetric/antisymmetric block contradiction (repository reduced SAT is redundant) |
| 90 | CLOSED | size-13 12+1 defect / affine-ADE obstruction; all-same unique C5^- + 5K2^+ orbit killed by eigenspace block/Cauchy contradiction |
| 99 | CLOSED | size-13 e=9 defect and e=18 P3/2K2 impossibility; all-same final two cospectral orbits A,B killed by local 15-eigenspace Gram mismatch |
| 108 | CLOSED | all-same exact repository theorem; size-13 e=9 defect, e=18 exact two-hub Schur, e=27 inverse-loss bound |
| 117 | CLOSED | all-same exact support theorem; size-13 exact inverse-loss / Schur bounds |
| 126 | CLOSED | complete record-level all-same support closure (cubic log majorant + exact weighted Fourier); all size-13 allocations exact below record |
| 135 | CLOSED | cubic log majorant + complete high-moment triangle-core enumeration and exact integer Eisenstein Fourier; size-13 below record |
| 144 | CLOSED | exact all-same K6 boundary/arithmetic + size-13 exact Schur/state pairing |
| 153 | CLOSED | exact tail final closure: all-same K6/isolate split + sparse size-13 two-hub Cauchy |
| 162 | CLOSED | exact tail structural reduction / trace-support certificate |

### Q == 6 (mod 9)

| Q | record-level status | main mechanism |
|---:|---|---|
| 6, 15, 24, 33, 42, 51 | CLOSED / impossible for strict counterexample | mixed-color cross-energy plus \(B_3(15,10)=12\); (14,1) needs at least 42 cross + 18 internal energy |
| 60 | CLOSED | twin-leaf PSD obstruction |
| 69 | CLOSED | twin-leaf PSD obstruction |
| 78 | CLOSED | equal-color-sign isolated-orbit obstruction |
| 87 | CLOSED | (13,2) three allocations: orthogonal-class code bound or 12+1 defect; (14,1)x(14,1) isolated-orbit obstruction |
| 96 | CLOSED | size-13 e=9 defect/e=18 P3,2K2 impossibility; remaining (14,1)x(14,1) supports have isolates and are excluded |
| 105 | CLOSED | size-13 e=9/e=18 structural impossibility and e=27 exact bound; (14,1)x(14,1) unique no-isolate 7K2 branch has maximum exactly D0 |
| 114 | CLOSED | exact size-13 residual/inverse-loss bounds + size-14 componentwise Schur bounds |
| 123 | CLOSED | exact size-13 hub losses + size-14 componentwise bounds |
| 132 | CLOSED | exact size-13 refinements + size-14 no-isolate component enumeration |
| 141 | CLOSED | exact size-13 residual + size-14 no-isolate support bounds |
| 150 | CLOSED | size-13 exact state pairing + size-14 componentwise closure, including e=108 K2 + T12 boundary |
| 159 | CLOSED | exact tail final closure: one size-13 residual by two-hub Cauchy; size-14 component configurations and exact P14 tree boundary |
| 168 | CLOSED | exact tail structural/trace-support certificate |

## New structural lemmas replacing solver-heavy branches

### 12+1 defect lemma

If 12 length-15 mu3 rows are pairwise orthogonal and a 13th row has a single nonzero same-color inner product with them, the rank-3 complement tight frame has affine-ADE orthogonality graph \(Om=2m\). Parseval on an affine-A2 basis (and local rank-3 obstructions for the remaining affine-ADE cases) excludes all attainable same-color defect norms needed here, including norm 9. This eliminates the e=9 size-13 branches at Q=87/90/96/99/105/108.

### Energy-18 near-equidistant lemma

For a 13-row same-color block of internal energy 18, the support is P3 or 2K2.

* P3: deleting the middle row gives 12 mutually orthogonal rows; affine-ADE + a 3x3 contingency condition leaves only local C4/C5/E6 or D4 cases. The former violate rank 3; D4 violates the exact ETF projection identity
  \[
  3\sum M_i^2=(\sum M_i)^2.
  \]
* 2K2: a rank-2 complement reduces three sign patterns immediately; the final negative-negative case leaves two exact 3x3 joint tables. Each gives 6579 candidate completion rows. The exact ternary orthogonality graph has 2,436,912 edges and no 9-clique (integer difference-count adjacency; no floating-point dependence).

Thus internal energy 18 is impossible for the relevant 13-row blocks.

### Q=99 final all-same local obstruction

The only pure-11 signed-degree survivors are
\[
A=4K_2^-\sqcup K_3^+\sqcup C_4^+,
\qquad
B=3K_2^-\sqcup K_3^+\sqcup T_6,
\]
with common characteristic polynomial
\[
(x-21)^2(x-18)^4(x-15)^2(x-12)^6(x-9).
\]
For A/B, the unique 15-eigenspace supports force the C4-by-T6 subblock Z to satisfy \(G_{C4}Z=ZG_{T6}\). There are exactly 144 mu3-valued solutions to this local intertwining equation; none has the residual Gram shape forced by the remaining columns. A/A similarly has 225 local mu3 commutant blocks and 0 residual-compatible states. B/B has 729 local mu3 commutant blocks and 0 residual-compatible states. Hence all four orbit pairings are impossible without MILP/Z3.

## Conclusion

Every admissible shell with \(Q\le168\) is record-level closed, and every \(Q\ge171\) lies below the benchmark by the exact trace/variance certificate. Since the benchmark matrix realizes \(D_0\),

\[
\boxed{\max_{H\in\mu_3^{15\times15}} |\det H|^2
=277868041444786176
=2^{22}3^{20}19.}
\]

Equivalently,
\[
\boxed{\max |\det H|=2^{11}3^{10}\sqrt{19}.}
\]

The remaining engineering task is to convert the new lemmas/enumerations into repository scripts/tests and replay them from a clean checkout; there is no remaining mathematical shell in the ledger.
