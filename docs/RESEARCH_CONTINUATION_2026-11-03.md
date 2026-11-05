> Historical version. The current results and regenerated certificates are documented in [the 5 November recovery report](RECOVERY_2026-11-05.md). Claims below describe the earlier state.

# Research continuation after Pádraig Ó Catháin's comments (3 November 2026)

This is a **working research ledger**, not a replacement for the existing
computer-assisted proof of the order-15 maximum. It distinguishes established
results supported by the repository or the previous exact-verification bundle
from newly reported computations whose complete input, program and certificate
still need to be checked in.

## Baseline

Let \(\omega=\exp(2\pi i/3)\), \(\mu_3=\{1,\omega,\omega^2\}\), and
\(Q(H)=\sum_{i<j}|(HH^*)_{ij}|^2\).

The order-15 exact maximum remains

\[
\max_{H\in\mu_3^{15\times15}}|\det H|^2
 =2^{22}3^{20}\cdot19=277868041444786176.
\]

The displayed maximizer has \(Q=105\) and a monomial automorphism group
of order 1008. These statements and the strict upper-bound certificate are
already documented in the baseline paper. The original proof is left intact.

## Research results with derivations or existing exact checks

1. **More conceptual Gram-graph argument.** Writing
   \(HH^*=nI+3T\), \(H^*H=nI+3S\) gives \(TH=HS\).
   Componentwise determinant costs
   \(c(A)=\log(n^{|A|}/\det A)\) make clear why one enumerates cheap
   Hermitian support graphs and then checks whether row/column supports
   can coexist. A previously prepared revised manuscript contains this
   organization and exact audit material; it has not yet replaced
   \`paper/main.tex\` on \`main\`.

2. **Fano/weighing-matrix construction.** The displayed order-15
   extremizer can be constructed from a balanced \(BW(7,4)\) and
   a seven-point Fano-plane support. Its Gram spectrum is
   \(12^{(7)},18^{(6)},(33+\sqrt{177})/2,(33-\sqrt{177})/2\).
   A general constant-border determinant identity is derived in the
   revised manuscript. This describes a *conditional* construction
   family; it does not establish the existence of a weighing matrix
   at every proposed parameter.

3. **An obstruction at the next parameter.** A preceding exact finite
   certificate excludes balanced \(BW(16,10)\), hence the proposed
   construction does not directly extend to order 33. The distinct
   existence problem for balanced \(BW(25,16)\) (which would give an
   order-51 instance) remains open in this work.

4. **General energy gaps.** The preceding revision proves
   \(Q\ge3(n-1)\) for each non-Butson matrix of order \(3\mid n\), with
   equality exactly for a one-entry alteration of \(BH(n,3)\).
   It additionally derives refined bounds using spectral projections
   and the inert-prime obstruction. The independent addendum in this
   PR supplies a self-contained proof that, at norm-obstructed orders
   \(n\), a *mixed row-coloring* already forces
   \(Q\ge3n+15\). This last assertion does not require any finite graph
   catalogue or the uncommitted classification computations.

## Further claims from the continuation run: verification pending

These are **research leads** and must not yet be cited as established
theorems or silently inserted into the abstract.

- **Exact minimum energy at order 15:** the continuation reported
  \(\min_{H\in\mu_3^{15\times15}}Q(H)=105\), with the stronger
  conditional estimate \(Q\ge108\) for the pure-color situation.
  The displayed maximizer shows that 105 is attainable. The still-needed
  item is a complete, reproducible *lower-bound* proof for all matrices,
  independent of the strict-determinant-threshold argument.
- **Equality-case classification:** the continuation reported nine
  monomial classes for one fixed Gram type and nine for its conjugate,
  suggesting **exactly 18** classes under row/column cube-root
  monomial equivalence, *without* identifying conjugation or transpose.
  The counts, class representatives, stabilizers, orbit sizes, exact
  finite enumeration and completeness of the Gram equality boundary
  must all be released before this can replace the current manuscript's
  explicit "not classified" qualification.
- **Orthogonality capacity:** a five-case clique computation was
  reported as a possible self-contained proof of \(M_3(15)\le9\).
  It may remove one external classification input when an independently
  reproducible maximum-clique certificate is included.
- **The order-51 support problem:** a continuation run reported 78
  isomorphism types of symmetric \((25,9,3)\) designs and a reduction
  to 50 binary sign systems with 75 variables each. This is a roadmap,
  not an existence or nonexistence proof for \(BW(25,16)\).

For each item above the acceptance bar is: explicit mathematical
statement and equivalence convention; exact reproducible program;
complete machine-checkable certificate or complete independent proof;
positive controls; and a documented run that reproduces all counts.
Until then the baseline theorem is not modified.

## Immediate next steps

- Bring the earlier full LaTeX revision and its replay suite into a
  separate reviewed commit, avoiding unnecessary changes to proven
  code paths.
- Archive the energy-105 and 18-class scripts together with input
  matrices, checksums, exact branch-coverage data and class witnesses.
- Check all published-source dependencies, especially difference
  matrices and symmetric designs, before upgrading provisional claims.
- Solve or conclusively obstruct the 50 sign systems relevant to
  balanced \(BW(25,16)\).

The present PR deliberately separates *new, checkable mathematics*
from *promising but not yet archived computational conclusions*.
