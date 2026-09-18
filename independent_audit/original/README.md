# Order-15 μ3：指定 Q-shell 归档包

本包整理的 shell：

`87, 90, 96, 99, 105, 108, 114, 117, 123, 126, 132, 135, 141, 144, 150`

它们在 2026-09-15 的最终 record-level ledger 中全部为 **CLOSED**。

## 推荐阅读顺序

1. `proofs/Q87_Q150_proof_digest_zh.md` — 针对这 15 个 Q 的中文证明整理。
2. `proofs/final_record_level_ledger.md` — 上个对话最终 ledger 原件。
3. `process/TIMELINE.md` — 关键研究过程、失败估计与修正。
4. `programs/verify_q99_final_orbits_exact.py` — Q99 独立 exact 证书。
5. `results/q99_certificate_output.json` — 本次打包时实际重跑结果。
6. `inventory/repository_program_inventory.md` — 原仓库 verifier/模块清单及哪些源码未随聊天保存。
7. `NOTES_archive_scope.md` — 说明本 ZIP 与完整本机 Git 仓库之间的边界。

## 目录

- `proofs/`：最终 ledger + 本次整理的 15-shell 证明。
- `programs/`：可恢复的程序与归档完整性检查器。
- `process/`：原始过程日志与整理后的时间线。
- `context/`：旧 research handoff，用于理解证明架构的演变。
- `results/`：状态表与程序重跑输出。
- `inventory/`：模块清单与 SHA-256 manifest。

## 注意

`context/research_handoff_20260914.md` 是**中间 checkpoint**，其中仍写有 “project is not solved yet”。它不能覆盖 2026-09-15 的最终 ledger；保留它是为了过程审计。
