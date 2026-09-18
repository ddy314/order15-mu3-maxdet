# Migration and comparison record

The independent local audit is the trusted baseline. Historical remote sources have been reviewed against it and retained under `historical_remote/` with explicit restrictions in [REMOTE_AUDIT.md](REMOTE_AUDIT.md). They are not accepted globally verified proofs.

Source main: `827672a2c60dad7064c7955d9c460261d866d25e`. Research branch: `007af7d8f7f644396aee27dd87c68c8f4169efd5`. All 186 order15-specific source files are copied byte-for-byte and hashed in `provenance/migration_manifest.json`. The separate source branch patch preserves all six changed paths. Shared order55 `boxed_moment` code and tests stay in order55; no order55-specific native code, certificates, or paper is moved.

All three local artifacts are preserved unchanged. Fixed ZIP inputs match the extracted audit; the separately delivered PDF and CSV match their ZIP copies. The source and local benchmark matrices have exactly the same 225 exponents, and their exact determinant computations agree. The local eight-step replay passes and retains 19 incomplete shells.

```bash
uv run python scripts/verify_migration.py
uv run python scripts/compare_remote_bounds.py
uv run python scripts/check_remote_orbit_orientation.py
uv run python independent_audit/run_all.py --audit-only
```

The order55 cleanup removes only the 186 order15-specific files, removes order15-only dependency extras and its environment section, and replaces the README's order15 description with the new repository address. Its order55 theorem, source, native verifiers, certificates, paper, and shared mathematical utility are unchanged. Four remaining Python tests pass. Cleanup is published only after the new repository has received the preserved artifacts and audit report. Git history and the original research branch are not rewritten or deleted.
