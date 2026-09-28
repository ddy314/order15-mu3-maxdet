# Structure of the displayed order-15 maximizer

Let `H_c` be the canonical representative stored in `data/benchmark_canonical.json`, with entries `omega^e` for the stored exponents `e`. Its row Gram matrix is

\[
G_c=
\begin{pmatrix}
I_7\otimes C & \delta\mathbf 1_{14}\\
\bar\delta\mathbf 1_{14}^{T} & 15
\end{pmatrix},\qquad
C=\begin{pmatrix}15&3\\3&15\end{pmatrix},\quad \delta=1-\omega.
\]

Thus the nonzero support graph consists of seven triangles meeting in one hub. Among the 105 unordered off-diagonal pairs there are 84 zeros, 14 entries of Eisenstein norm 3, and 7 entries of norm 9, giving `Q=105`.

The spectrum is

\[
12^{(7)},\quad 18^{(6)},\quad \frac{33-\sqrt{177}}2,\quad \frac{33+\sqrt{177}}2,
\]

and therefore

\[
\det G_c=12^7 18^6\cdot 228=277868041444786176.
\]

For this particular matrix, the largest mutually orthogonal row set has size 7. Every maximum set chooses one endpoint from each of the seven matching edges, so there are exactly `2^7 = 128` such sets. This statistic is distinct from the global length-15 capacity `M_3(15)=9` supplied by the difference-matrix classification.

Under the convention that an automorphism is a pair of `mu_3`-monomial matrices `(M,N)` satisfying `M H_c N* = H_c`, without adjoining conjugation or transposition, exact enumeration gives:

- projective monomial automorphism group: order `336`, isomorphic to `C2 x GL(3,2)`;
- full monomial-pair automorphism group: order `1008`, isomorphic to `C6 x GL(3,2)`.

The induced action on the seven matching pairs contains the order-168 Fano-plane automorphism group. The verifier examines all `7! * 2^7 = 645120` permitted row permutations and checks all 225 matrix entries for each retained automorphism. It regenerates the complete list of 336 projective automorphisms in `data/automorphisms.json`.

These statements describe one displayed maximizer. They do not classify all matrices attaining the maximal determinant and do not prove uniqueness of the Hadamard-equivalence class.
