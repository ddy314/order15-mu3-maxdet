# Order-15 maximal determinant over the third roots of unity

[![Exact verification](https://github.com/ddy314/order15-mu3-maxdet/actions/workflows/verify.yml/badge.svg)](https://github.com/ddy314/order15-mu3-maxdet/actions/workflows/verify.yml)

This repository contains a computer-assisted proof of the exact maximal determinant of a $15\times15$ matrix with entries in

$$
\mu_3=\{1,\omega,\omega^2\},\qquad \omega^2+\omega+1=0.
$$

The main theorem is

$$
\max_{H\in\mu_3^{15\times15}} |\det H|^2
=277868041444786176
=2^{22}3^{20}19,
$$

equivalently

$$
\max |\det H|=120932352\sqrt{19}.
$$

> **Status.** The maximal value is proved exactly. The repository does **not** claim that the maximizing Hadamard-equivalence class is unique.

## At a glance

| Item | Exact result |
| --- | --- |
| Maximum squared determinant | $2^{22}3^{20}19$ |
| Maximum determinant | $120932352\sqrt{19}$ |
| Sharp length-15 orthogonality capacity | $M_3(15)=9$ |
| Finite Gram-energy candidates after the sharp reduction | 23 |
| Largest orthogonal row subset in the displayed maximizer | 7 rows |
| Projective monomial automorphism group | order 336, $C_2\times\mathrm{GL}(3,2)$ |
| Full $\mu_3$-monomial-pair automorphism group | order 1008, $C_6\times\mathrm{GL}(3,2)$ |

## What changed in the September 2026 revision

The revised manuscript reorganizes the proof around the Gram energy

$$
Q=\sum_{1\le i<j\le15}|(HH^*)_{ij}|^2,
$$

and credits the classical Gram-matrix search strategy of Moyssiadis--Kounias and later maximal-determinant work. The main new structural input is the published difference-matrix classification of Lampio--Östergård, which gives the exact orthogonality capacity $M_3(15)=9$. This is stronger than the ordinary ternary equidistant-code parameter $B_3(15,10)=12$. The exact orthogonality bound reduces the coarse finite energy list from **38 shells to 23**, and directly eliminates many sparse Gram candidates that previously required heavier decomposition tests.

The displayed extremizer is analyzed explicitly. After dephasing and monomial equivalence, its row Gram matrix has the form

$$
G_c=
\begin{pmatrix}
I_7\otimes\begin{pmatrix}15&3\\3&15\end{pmatrix} & (1-\omega)\mathbf1_{14}\\
(1-\omega^2)\mathbf1_{14}^{T} & 15
\end{pmatrix}.
$$

Its support graph is seven triangles sharing one common vertex. Its spectrum is

$$
12^{(7)},\qquad 18^{(6)},\qquad \frac{33\pm\sqrt{177}}2,
$$

so

$$
\det G_c=12^7\,18^6\,228=277868041444786176.
$$

For this particular maximizer, the largest mutually orthogonal row set has size 7, with exactly $2^7=128$ maximum such subsets. The projective monomial automorphism group has order 336 and is isomorphic to $C_2\times\mathrm{GL}(3,2)$; restoring the common scalar subgroup gives a full $\mu_3$-monomial-pair automorphism group of order 1008, isomorphic to $C_6\times\mathrm{GL}(3,2)$.

## Read the paper

The article is split into a short entry file and focused section files:

- [paper/main.tex](paper/main.tex) — preamble, theorem statement, and section assembly;
- [paper/sections/](paper/sections/) — mathematical sections and appendices;
- [paper/matrix_dephased.tex](paper/matrix_dephased.tex) — the displayed normalized exponent matrix.

The paper now includes the historical Gram-search background, the sharp orthogonality theorem, the 23-shell reduction, a uniform finite-shell framework, the normalized extremizer and Gram spectrum, the automorphism computation, and a discussion of what generalizes beyond order 15.

Compile with a standard LaTeX installation:

    cd paper
    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex

## Verify the result

Python dependencies are pinned in [proof/requirements.txt](proof/requirements.txt). A C++17 compiler is needed only for the complete automorphism enumeration.

    python -m venv .venv
    . .venv/bin/activate
    python -m pip install -r proof/requirements.txt
    python verify_all.py

The unified verifier runs three layers in order:

1. [proof/run_final.py](proof/run_final.py) — the inherited exact 38/38 proof replay;
2. [tools/verify_sharp_reduction.py](tools/verify_sharp_reduction.py) — the sharper 38 → 29 → 23 reduction and explicit independent-set/spectral checks;
3. [tools/verify_structure.py](tools/verify_structure.py) — exact determinant, Gram, orthogonality, and automorphism checks for the displayed extremizer.

Assertions must remain enabled; do not use **python -O**. No GPU or numerical optimizer is required. The GitHub Actions workflow runs the same unified verifier on changes to the mathematical or verification sources.

## Repository layout

| Path | Purpose |
| --- | --- |
| [paper/](paper/) | Revised, split LaTeX manuscript and displayed matrix |
| [proof/](proof/) | Inherited replayable computer-assisted maximality certificate |
| [tools/](tools/) | Sharp-reduction and extremizer-structure verifiers |
| [data/](data/) | Original, dephased, and canonical benchmark exponent matrices |
| [docs/STRUCTURE.md](docs/STRUCTURE.md) | Compact structural summary of the extremizer |
| [docs/REVISION_NOTES_ZH.md](docs/REVISION_NOTES_ZH.md) | Chinese revision and verification notes |
| [verify_all.py](verify_all.py) | Single command for the complete verification stack |

The historical proof code under **proof/** is intentionally kept in its existing layout because its stages import one another through those paths. The cleaner top-level **tools/**, **data/**, and **docs/** directories contain the new results without disturbing the inherited certificate.

## External mathematical inputs

Two finite classification results are cited rather than reimplemented:

- T. Todorov and G. Bogdanova, *Ternary equidistant codes of length 11 ≤ n ≤ 15*, J. Math. Comput. Sci. 10 (2020), 2713--2721. This gives $B_3(15,10)=12$.
- P. H. J. Lampio and P. R. J. Östergård, *Classification of difference matrices over cyclic groups*, J. Statist. Plann. Inference 141 (2011), 1194--1207, doi:10.1016/j.jspi.2010.09.023. This gives the exact orthogonality capacity $M_3(15)=9$.

The revised paper also credits the real maximal-determinant Gram-search literature, including Moyssiadis--Kounias and Orrick.

## Scope

The theorem determines the exact maximum. The finite certificate remains a computer-assisted proof rather than a proof-assistant formalization. The external difference-matrix and code classifications are not rerun here. Equality cases have not been completely classified, and no uniform asymptotic gap from the Hadamard bound is claimed for an infinite congruence class.
