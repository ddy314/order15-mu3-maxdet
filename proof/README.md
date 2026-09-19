# Proof certificate

This directory is the replayable certificate for the order-15 third-root maximal determinant theorem. The only public entry point needed for a full verification is:

```bash
python run_final.py
```

Run it from this directory or as `python proof/run_final.py` from the repository root after installing `requirements.txt`.

The source is intentionally divided into three historical stages because the verified programs import one another through these paths. `order15_mu3_unified_audit/` reconstructs the exact arithmetic, benchmark and first 19 closed shells; `new/` reconstructs the component, Schur, projection, rank, characteristic-polynomial and pure-color catalogues used to reach 35/38; `last_three/` regenerates the final `Q=96,99,105` candidates, performs two independent complete column-equation scans, verifies exact annihilators, proves the residual `Q=105` obstruction, and joins all shells into the final theorem.

Generated JSON files and replay logs are deliberately absent from Git. The scripts write them locally during verification. The compact delivered snapshots in `certificates/` record the final theorem and the reference full replay, but `run_final.py` regenerates the proof rather than treating those snapshots as trusted input.

The benchmark exponent matrix is `order15_mu3_unified_audit/data/benchmark.json`. Its exact determinant is

```text
604661760 + 241864704 * omega
```

with Eisenstein norm

```text
277868041444786176 = 2^22 * 3^20 * 19.
```

The computer-assisted theorem depends externally on the published classification result `B_3(15,10)=12` of Todorov--Bogdanova (2020). No uniqueness classification of maximizing matrices is asserted.