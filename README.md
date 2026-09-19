# Order-15 maximal determinants over third roots of unity

Independent research on 15 × 15 matrices over {1, ω, ω²}, where ω² + ω + 1 = 0. Extracted from [order55-binary-circulant-maxdet](https://github.com/dongxuelian2/order55-binary-circulant-maxdet); the order-55 binary-circulant theorem is a different problem.

The complete 38/38 computer-assisted proof is now included in [`complete_proof/`](complete_proof/README_FINAL_zh.md). It proves global optimality:

```text
max |det H|² = 277868041444786176 = 2²² · 3²⁰ · 19
max |det H| = 120932352 · √19
```

The final package closes the three former frontier shells `Q=96, 99, 105`, combines them with the independently replayed 35 shells and the `Q≥171` tail bound, and records no remaining shells. The result is a computer-assisted mathematical proof, not a proof-assistant formalization; it also does not classify all equality matrices. The external dependency `B₃(15,10)=12` is cited in the final report and was not re-proved here.

## Reproduce

```bash
uv sync --extra dev --extra sat
uv run pytest -q
uv run python independent_audit/run_all.py --audit-only
```

The first command installs the research and solver dependencies. The audit requires assertions enabled. A passing historical test suite validates its assertions, not necessarily comparison against the record. For the completed result, install `complete_proof/requirements.txt` in an isolated environment and run:

```bash
cd complete_proof
python run_final.py
```

The final entry point reruns the inherited audit and all new finite certificates with two independent exhaustive column algorithms. It returns success only when all 38 shells are closed and the certificate has no remaining obligations; the delivered replay took about ten minutes in its recorded environment.

## Evidence and trust boundary

The complete package is now the global proof record. The independent audit remains the trusted replay basis for the inherited shells, while remote research sources are retained as historical material; passing historical tests alone does not certify their mathematical claims.

- [Independent audit](independent_audit/README_zh.md): rebuilt exact verifiers, shell status, logs and remaining obligations for the former 19/38 baseline.
- [Complete 38/38 proof](complete_proof/README_FINAL_zh.md): final report, source programs, exact certificates, replay logs and PDFs.
- [Audit report](independent_audit/report/order15_mu3_independent_audit.pdf).
- [Original research](historical_remote/order15_mu3/README.md), [historical handoff](historical_remote/docs/HANDOFF.md), and [historical frontier](historical_remote/docs/REMAINING_FRONTIER.md): preserved from the source main commit, including claims not yet independently accepted.
- [Remote proof audit](REMOTE_AUDIT.md): scope errors, conjugation issue, replay results and claims withheld.
- [Migration verification](MIGRATION.md) and [source hashes](provenance/migration_manifest.json).
- `provenance/research-reaggregate.patch` preserves the unmerged research branch changes, including CI definitions, without silently applying them.
- `artifacts/` retains the earlier delivered local files unchanged. The complete package keeps its original zip under `provenance/` and its unpacked, replayable sources under `complete_proof/`.

No unrestricted order-15 theorem follows from the order-55 result. The completed result instead combines explicit exact finite certificates with the documented structural and infinite-tail arguments in `complete_proof/report/FINAL_PROOF_zh.md`.
