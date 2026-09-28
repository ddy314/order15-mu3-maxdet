# Proof certificate

This directory contains the inherited replayable certificate for the exact order-15 third-root maximal determinant theorem. Its public entry point remains

```bash
python proof/run_final.py
```

The historical three-stage layout is preserved because the verified programs import one another through these paths:

- `order15_mu3_unified_audit/` reconstructs exact arithmetic, the benchmark, and the low-energy baseline;
- `new/` reconstructs the component, Schur, projection, rank, characteristic-polynomial, and pure-color catalogues;
- `last_three/` reconstructs the former terminal `Q=96,99,105` candidates, runs the two independent complete column-equation scans, and aggregates the global 38/38 result.

The September 2026 revision adds a stronger published orthogonality input, `M_3(15)=9`, which makes several of those terminal calculations redundant in the *new presentation*. They are deliberately retained here as an independent regression layer: the revised proof becomes shorter, while the older exact computations still replay and agree.

The new checks live outside this historical tree:

```bash
python tools/verify_sharp_reduction.py
python tools/verify_structure.py
```

or run everything with

```bash
python verify_all.py
```

Generated JSON files and replay logs under the historical proof directories are normally regenerated locally rather than treated as trusted evidence. Compact structural inputs for the revised presentation are kept under top-level `data/`.

The benchmark exponent matrix in the inherited audit has exact determinant

```text
604661760 + 241864704 * omega
```

with Eisenstein norm

```text
277868041444786176 = 2^22 * 3^20 * 19.
```

The finite proof uses published external classification results; see the root README and the paper for exact citations. No uniqueness classification of maximizing matrices is asserted.
