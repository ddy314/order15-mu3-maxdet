# 15 阶 μ₃ 最大行列式：2026-09-18 进展包

**本次新增闭合 16 个完整壳层，总计 35/38；尚余 Q=96,99,105。全局最大性仍未证明。** 纪录平方值为 B=277868041444786176=2²²·3²⁰·19。新增闭合的是 87,90,108,114,117,123,126,132,135,141,144,150,153,159,162,168。完整证明与枚举覆盖说明见 `report/PROOF_zh.md` 和 PDF；逐壳层结果见 `results/all_Q_progress.csv`。

这份工作保留内部孤立点，以分量 Schur 界、双侧核投影、低秩整数多项式、精确相位目录及新的单边缺陷引理替换遗漏分支的旧论证。Q108 扩展目录实际计算了 15,899,506 个相位状态。结构引理有书面证明，程序负责有限的算术检查与目录复演。原审计引用的已发表码界 B₃(15,10)=12 继续作为外部定理使用；本包没有重做原论文的码分类。

## 运行

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_progress.py --include-baseline --audit-only
```

默认重建全部新目录；`--include-baseline` 同时重跑原审计八项基础验证。`python run_progress.py --cached --audit-only` 会先核验目录哈希，再重算下游证明，不重新生成相位目录。去掉 `--audit-only` 后，只要仍有未闭合壳层就返回退出码 2。禁止 `python -O`。依赖版本是本次实际执行环境的版本。

## 证据与文件

`results/progress_ledger.json` 是本次状态汇总；`new/` 保存全部新源码、目录和精确结果；`report/` 说明每个手工引理及有限覆盖；`provenance/` 保存源码和目录哈希、来源与原 ZIP。`order15_mu3_unified_audit/` 中的 19/19 状态是原审计历史，不能当成本次总表。总表将继承证明与本次新增证明分开列出。

尚余三层的具体对象和候选值写在报告第 14 节及 `new/size13_energy18.json`。这些抽象 Gram 候选尚未证明能实现为 μ₃ 矩阵，也未被全部排除。没有修改远端仓库；可将本目录完整加入 `ddy314/order15-mu3-maxdet` 作为独立进展目录，避免覆盖历史证据。
