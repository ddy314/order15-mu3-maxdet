# Order-15 maximal determinant over the third roots of unity

This repository contains a computer-assisted proof of the exact maximal determinant for 15 x 15 matrices with entries in

\[
\mu_3=\{1,\omega,\omega^2\},\qquad \omega^2+\omega+1=0.
\]

The main result is

\[
\max_{H\in\mu_3^{15\times15}} |\det H|^2
=277868041444786176
=2^{22}3^{20}19,
\]

or equivalently

\[
\max |\det H|=120932352\sqrt{19}.
\]

The benchmark matrix attaining this value is included in `proof/order15_mu3_unified_audit/data/benchmark.json`. The proof closes all 38 admissible finite Gram-energy shells and the infinite tail. It uses exact Eisenstein-integer arithmetic, published ternary equidistant-code bounds, structural Gram estimates, and finite exhaustive certificates. It does **not** claim a classification of all equality cases and is not a proof-assistant formalization.

## Paper

The current article source is [`paper/main.tex`](paper/main.tex):

> **The Maximal Determinant of Order 15 over the Third Roots of Unity**

It is intended as the readable mathematical account of the result. Compile with a standard LaTeX installation:

```bash
cd paper
pdflatex main.tex
pdflatex main.tex
```

Author metadata is deliberately left blank in the draft and should be filled in before formal circulation or submission.

## Reproduce the proof

Python 3.11+ is recommended. From the repository root:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r proof/requirements.txt
python proof/run_final.py
```

`run_final.py` regenerates the low-energy audit, the intermediate component/spectral catalogues, and the final `Q=96,99,105` certificates. On the reference run the complete replay took about ten minutes. Assertions must remain enabled; do not use `python -O`.

The final delivered ledger and replay summary are retained in `proof/certificates/` as a human-readable snapshot. The replay does not trust generated positive outputs: intermediate JSON catalogues and logs are ignored by Git and rebuilt from source.

## Repository layout

- `paper/` -- the article source.
- `proof/run_final.py` -- single full-proof entry point.
- `proof/order15_mu3_unified_audit/` -- exact arithmetic and the independently rebuilt low-energy baseline.
- `proof/new/` -- component/Schur, projection, rank and spectral certificates that close the intermediate shells.
- `proof/last_three/` -- exact certificates for `Q=96,99,105` and the final global aggregator.
- `proof/certificates/` -- compact snapshots of the final 38/38 ledger and the successful full replay.

Historical migration notes, superseded audit reports, development logs, duplicate archives, generated catalogues, and old partial PDFs were removed from the current tree after the proof was completed. They remain recoverable from Git history.

## External dependency

The proof uses the published value `B_3(15,10)=12` from T. Todorov and G. Bogdanova, *Ternary equidistant codes of length 11 <= n <= 15*, J. Math. Comput. Sci. 10 (2020), 2713-2721, doi:10.28919/jmcs/4964. That finite code classification is cited rather than re-proved here.