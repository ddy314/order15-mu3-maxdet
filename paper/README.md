# Manuscript

The paper source is intentionally split so the mathematical argument can be reviewed section by section.

| File | Contents |
| --- | --- |
| [main.tex](main.tex) | Preamble, title, abstract, and section assembly |
| [sections/01_introduction.tex](sections/01_introduction.tex) | Statement of the theorem, literature background, Gram-search strategy |
| [sections/02_arithmetic.tex](sections/02_arithmetic.tex) | Eisenstein arithmetic, orthogonality, sharp 23-shell reduction |
| [sections/03_search_tools.tex](sections/03_search_tools.tex) | Schur, moment, rank, projection, and spectral-rectangle tools |
| [sections/04_finite_verification.tex](sections/04_finite_verification.tex) | Uniform treatment of all remaining energy shells |
| [sections/05_extremizer.tex](sections/05_extremizer.tex) | Gram structure, spectrum, orthogonal subsets, automorphism group |
| [sections/06_generalization.tex](sections/06_generalization.tex) | What extends to other orders and alphabets |
| [sections/07_reproducibility.tex](sections/07_reproducibility.tex) | Reproduction commands, external inputs, scope |
| [sections/A_matrix.tex](sections/A_matrix.tex) | Explicit dephased extremizer and switching to canonical form |
| [sections/B_q81.tex](sections/B_q81.tex) | Detailed boundary argument for Q=81 |
| [sections/C_replay.tex](sections/C_replay.tex) | Map from proof obligations to executable verifiers |
| [sections/references.tex](sections/references.tex) | Bibliography |
| [matrix_dephased.tex](matrix_dephased.tex) | Displayed 15x15 exponent matrix |

Compile from this directory with three LaTeX passes:

    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -interaction=nonstopmode -halt-on-error main.tex

The generated PDF is intentionally not tracked. The mathematical source, exact benchmark data, and replayable verification programs are the maintained artifacts.
