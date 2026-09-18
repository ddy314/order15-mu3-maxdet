# Order-15 maximal determinants over third roots of unity

Independent research on 15 × 15 matrices over {1, ω, ω²}, where ω² + ω + 1 = 0. Extracted from [order55-binary-circulant-maxdet](https://github.com/dongxuelian2/order55-binary-circulant-maxdet); the order-55 binary-circulant theorem is a different problem.

**Global optimality remains unproved by the independent audit.** The exactly verified record has |det H|² = 277868041444786176 = 2²² · 3²⁰ · 19. The audit closes 19 of 38 relevant energy shells; 19 remain incompletely verified: 87, 90, 96, 99, 105, 108, 114, 117, 123, 126, 132, 135, 141, 144, 150, 153, 159, 162, 168. Historical closure claims and passing tests do not override this audit status.

## Reproduce

```bash
uv sync --extra dev --extra sat
uv run pytest -q
uv run python independent_audit/run_all.py --audit-only
```

The first command installs the research and solver dependencies. The audit requires assertions enabled. Without `--audit-only`, its runner exits 2 while global proof obligations remain; this is intentional. A passing historical test suite validates its assertions, not necessarily comparison against the record.

## Evidence and trust boundary

The delivered local audit is the trusted baseline. Remote research sources are retained as historical, unverified material; passing their tests alone does not certify their mathematical claims.

- [Independent audit](independent_audit/README_zh.md): rebuilt exact verifiers, shell status, logs and remaining obligations; authoritative for independently replayed scope.
- [Audit report](independent_audit/report/order15_mu3_independent_audit.pdf).
- [Original research](historical_remote/order15_mu3/README.md), [historical handoff](historical_remote/docs/HANDOFF.md), and [historical frontier](historical_remote/docs/REMAINING_FRONTIER.md): preserved from the source main commit, including claims not yet independently accepted.
- [Remote proof audit](REMOTE_AUDIT.md): scope errors, conjugation issue, replay results and claims withheld.
- [Migration verification](MIGRATION.md) and [source hashes](provenance/migration_manifest.json).
- `provenance/research-reaggregate.patch` preserves the unmerged research branch changes, including CI definitions, without silently applying them.
- `artifacts/` retains the three delivered local files unchanged. Replays update files under `independent_audit/results/`.

No unrestricted order-15 theorem follows from the order-55 result. Finite enumeration certifies only its explicit domain; incomplete shell coverage remains a mathematical gap.
