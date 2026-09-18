# 程序与证书清单

## 已实际包含在本 ZIP 中

### `programs/verify_q99_final_orbits_exact.py`
独立、纯 Python、精确 Eisenstein arithmetic 的 Q=99 最终 orbit obstruction。
本次打包时已经重新运行；输出见 `results/q99_certificate_output.json`。

### `programs/verify_bundle_integrity.py`
仅检查本归档自身是否完整：
- 15 个指定 Q 是否都出现在最终 ledger 且标记为 CLOSED；
- Q99 独立证书是否运行通过。
它不是其余 shell 的数学 verifier。

## 在过程日志/原仓库中被明确引用，但聊天文件库没有保存源码快照

- `order15_mu3.scripts.verify_aggregate_upper`
- `order15_mu3.scripts.z3_q90_top_commuting`
- `order15_mu3.scripts.verify_pre_q144_audit`
- Q141 certificate / Schur phase modules
- Q126 K4-forest / triangle / unicyclic / bicyclic exact catalogue modules
- energy-9 top/second obstruction modules
- Q135 support / spectral-cap certificate modules
- Q144 / Q150 state-pairing and componentwise closure modules

这些模块在当时位于用户本机仓库（日志中的典型路径为 `E:\maximal determinant`）。
本 ZIP 只包含聊天/Library 中能恢复出来的源码与过程材料，不能冒充整个本机 Git 仓库快照。

## 主要重放命令（来自 handoff）

```powershell
python -m order15_mu3.scripts.verify_aggregate_upper
```

历史日志中也出现过：

```powershell
& 'E:\maximal determinant\.venv\Scripts\python.exe' -m order15_mu3.scripts.verify_aggregate_upper
```

以及 Q90 的 solver 路线：

```powershell
python -m order15_mu3.scripts.z3_q90_top_commuting --timeout-seconds 60
```
