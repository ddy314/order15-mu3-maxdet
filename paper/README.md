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
| [sections/06_generalization.tex](sections/06_generalization.tex) | Pure and mixed energy gaps, first equality layer, and odd norm obstruction |
| [sections/09_weighing_construction.tex](sections/09_weighing_construction.tex) | General border formula and explicit order-51 construction |
| [sections/D_energy_recovery.tex](sections/D_energy_recovery.tex) | Minimum energy and inclusive equality audit |
| [sections/E_orthogonality_recovery.tex](sections/E_orthogonality_recovery.tex) | Independent capacity-nine certificate |
| [sections/F_factorizations_recovery.tex](sections/F_factorizations_recovery.tex) | Complete Fano factorization and eighteen classes |
| [sections/G_explicit_recovery.tex](sections/G_explicit_recovery.tex) | Order-25 and sixth-root matrices |
| [sections/H_sixth_improvement.tex](sections/H_sixth_improvement.tex) | New sixth-root energy bound Q >= 10 |
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

## arXiv submission files

From the repository root, run:

    make arxiv

This compiles the manuscript and creates:

- `output/pdf/order15_arxiv_source.zip`: the LaTeX files required to build the paper, with `main.tex` at the archive root. Upload this archive to arXiv and select PDFLaTeX.
- `output/pdf/order15_arxiv.pdf`: a separate reading copy for checking the compiled result.

The source archive excludes build products, repository documentation, and revision notes. All references are embedded in `sections/references.tex`, so no BibTeX step or external bibliography file is needed. The mathematical verification programs remain available in the companion repository linked in the paper.

The author affiliation is University of Electronic Science and Technology of China, Chengdu, China. A one-sentence AI use declaration discloses AI participation in the research and the author's responsibility for all content, following the principles in [arXiv's official policy](https://info.arxiv.org/help/moderation/index.html#policy-for-authors-use-of-generative-ai-language-tools).
