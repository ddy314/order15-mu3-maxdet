# Maximal determinant of a 15 × 15 matrix over the third roots of unity

[![Exact verification](https://github.com/ddy314/order15-mu3-maxdet/actions/workflows/verify.yml/badge.svg)](https://github.com/ddy314/order15-mu3-maxdet/actions/workflows/verify.yml)

How large can the absolute value of a determinant be when every entry of a 15 × 15 matrix is one of the three complex cube roots of 1?

This repository contains a mathematical paper, an explicit matrix attaining the maximum, and programs that check the finite calculations in a computer-assisted proof. **The exact maximum is 120932352√19.** The matrix was already known; the contribution here is the proof that no allowed matrix has a larger determinant magnitude.

## The problem and the result

The three allowed entries are 1, ω, and ω², where ω = exp(2πi/3). They lie equally spaced on the complex unit circle and satisfy 1 + ω + ω² = 0. The notation **μ₃** refers to this set; **order 15** means that the matrix has 15 rows and 15 columns.

Because the determinant can be complex, we maximize its absolute value. The theorem is:

```math
\max_{H\in\{1,\omega,\omega^2\}^{15\times15}} \left|\det H\right|
=120932352\sqrt{19}.
```

For exact arithmetic, it is more convenient to use the square of this value:

```math
\max_{H\in\{1,\omega,\omega^2\}^{15\times15}} \left|\det H\right|^2
=277868041444786176
=2^{22}\cdot3^{20}\cdot19.
```

The proof has two parts: exhibit a matrix reaching this value, then exclude every possible strict improvement. **It does not classify all matrices attaining the maximum or prove that they all belong to one equivalence class.**

## How the proof works

Searching through all matrices directly would mean considering 3²²⁵ possibilities. Instead, the proof studies their **Gram matrices**: the matrices of inner products between rows.

For a matrix H, write G = HH<sup>∗</sup>, where H<sup>∗</sup> is the conjugate transpose. Every diagonal entry of G is 15, and det(G) = |det(H)|². An off-diagonal entry is zero exactly when the corresponding two rows are orthogonal.

The proof measures the total squared magnitude of these off-diagonal entries using a quantity called **Gram energy**, denoted by Q. Each unordered pair of rows is counted once:

```math
Q=\sum_{1\le i<j\le15}\left|G_{ij}\right|^2.
```

The argument proceeds as follows:

1. **Establish the lower bound.** Compute the determinant of the supplied matrix exactly.
2. **Restrict any hypothetical improvement.** Determinant bounds and arithmetic constraints limit the possible Gram energies.
3. **Use an orthogonality bound.** A published classification shows that at most nine length-15 vectors with entries in μ₃ can be mutually orthogonal. This leaves **23 possible energy values** to examine. An energy value is called a *shell* in the paper and verification logs.
4. **Exclude the remaining candidates.** Mathematical bounds and exact finite computations show that no candidate above the lower bound can be the Gram matrix of an allowed H.

Together, these steps establish the maximum. The paper explains the mathematical reductions; the verification programs regenerate and check the finite candidate families.

## A matrix attaining the maximum

The matrix is stored as a 15 × 15 array of exponents in [data/benchmark_normalized.json](data/benchmark_normalized.json). Replace each stored exponent e by ωᵉ: **0 means 1, 1 means ω, and 2 means ω²**. Its first row and first column are all 1; this normalization is called *dephasing*.

The same matrix is displayed in [the paper's matrix appendix](paper/sections/A_matrix.tex). [data/benchmark_canonical.json](data/benchmark_canonical.json) gives an equivalent representative obtained by permuting rows and columns and multiplying them by cube roots of unity. These operations preserve the determinant magnitude.

This representative has a simple pattern of row inner products: if two rows are joined whenever their inner product is nonzero, the resulting graph is **seven triangles sharing one central vertex**. Its determinant, spectrum, orthogonal row subsets, and symmetries are checked by [tools/verify_structure.py](tools/verify_structure.py).

The displayed matrix has at most **seven** mutually orthogonal rows. This is a property of this particular matrix; the bound of **nine** above concerns all length-15 vectors over μ₃. See [the structural summary](docs/STRUCTURE.md) for the detailed Gram matrix and symmetry groups.

## Verify the result

Run these commands from the repository root. Use **Python 3.11**, as in CI, and have a **C++17 compiler** available (`g++` by default). The compiler is used for the complete symmetry enumeration. No GPU is required.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r proof/requirements.txt
python verify_all.py
```

The dependencies are pinned in [proof/requirements.txt](proof/requirements.txt). Keep Python assertions enabled: do not run with `python -O` or set `PYTHONOPTIMIZE`.

The unified verifier runs three stages:

| Stage | What it checks |
| --- | --- |
| [proof/run_final.py](proof/run_final.py) | The complete maximality certificate, covering all 38 energy values in the initial, coarser reduction |
| [tools/verify_sharp_reduction.py](tools/verify_sharp_reduction.py) | The stronger reduction to 23 energy values and the additional exact exclusion checks |
| [tools/verify_structure.py](tools/verify_structure.py) | The supplied matrix, its exact determinant, Gram structure, orthogonal subsets, and symmetries |

The 38-case replay is retained as an additional verification route. The paper uses the sharper 23-case reduction; these are two reductions of the same problem.

On success, the unified verifier ends with:

```text
SUCCESS: maximality, sharp reduction, and extremizer structure all verified.
```

The programs regenerate intermediate files and verification summaries locally. See [proof/README.md](proof/README.md) for the certificate layout and [data/README.md](data/README.md) for the generated data. The [GitHub Actions workflow](https://github.com/ddy314/order15-mu3-maxdet/actions/workflows/verify.yml) runs the same unified verifier.

## Read the paper

Start with [the introduction](paper/sections/01_introduction.tex) for the theorem and motivation. The full manuscript is assembled by [paper/main.tex](paper/main.tex); [paper/README.md](paper/README.md) maps its sections to their source files.

To build the PDF, install a LaTeX distribution providing `pdflatex`, then run:

```bash
make paper
```

This produces `paper/main.pdf`. The PDF is generated locally and is not tracked in the repository.

## Mathematical inputs and scope

This is a **computer-assisted proof**, with mathematical arguments in the paper and exact finite checks in the code. It is not a proof-assistant formalization.

Two published classification results are used as external inputs; their classifications are not rerun by this repository:

- **Lampio–Östergård (2011)**, *Classification of difference matrices over cyclic groups*, J. Statist. Plann. Inference 141, 1194–1207. Supplies the sharp bound of nine mutually orthogonal length-15 vectors over μ₃. [DOI](https://doi.org/10.1016/j.jspi.2010.09.023).
- **Todorov–Bogdanova (2020)**, *Ternary equidistant codes of length 11 ≤ n ≤ 15*, J. Math. Comput. Sci. 10, 2713–2721. Supplies the weaker code bound of twelve used in the coarser reduction.

The Gram-matrix search method follows earlier maximal-determinant work by Moyssiadis–Kounias and others. The known attaining matrix and the methodological background are credited in [the introduction](paper/sections/01_introduction.tex) and [bibliography](paper/sections/references.tex).

The result settles the maximum for order 15. Classification of all maximizing matrices and a uniform asymptotic gap from the Hadamard bound remain outside the proved claims.

## Repository guide

| Path | Contents |
| --- | --- |
| [paper/](paper/) | LaTeX manuscript and explicit matrix appendix |
| [proof/](proof/) | Complete replayable maximality certificate |
| [tools/](tools/) | Sharp-reduction and matrix-structure verifiers |
| [data/](data/) | Exact matrix inputs, stored as exponents |
| [docs/STRUCTURE.md](docs/STRUCTURE.md) | Detailed structure and symmetries of the supplied matrix |
| [verify_all.py](verify_all.py) | Entry point for all verification stages |

For the history of the manuscript changes, see [the Chinese revision notes](docs/REVISION_NOTES_ZH.md).
