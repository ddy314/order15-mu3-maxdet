# Structural data

The tracked JSON files in this directory are the exact exponent matrices used by the revised paper:

- `benchmark_original.json` -- the historical extremizer;
- `benchmark_normalized.json` -- the dephased representative;
- `benchmark_canonical.json` -- a monomially equivalent representative with the structured Gram matrix used in the paper.

The verifiers regenerate three larger derived files locally:

- `structure_verified.json` -- determinant, Gram, support statistics, orthogonal-subset counts, and automorphism-group summary;
- `automorphisms.json` -- all 336 projective monomial automorphisms checked entry-by-entry;
- `sharp_reduction_verified.json` -- the exact 38/29/23 energy sets and explicit witnesses for the sharper reduction.

These derived files are intentionally ignored by Git because the programs are the source of truth. Regenerate them with

```bash
python proof/run_final.py
python tools/verify_sharp_reduction.py
python tools/verify_structure.py
```
