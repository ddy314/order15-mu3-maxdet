# 15 阶三次单位根最大行列式：38 / 38 完整计算机辅助证明

本包在原 35 层复核之上补齐 Q=96、99、105。合并原能量、颜色约化和 Q≥171 的全尾部界，得到

`max |det H|² = 277868041444786176 = 2^22 * 3^20 * 19`，其中 `H ∈ {1,ω,ω²}^(15×15)`。

相应最大模为 `120932352 * sqrt(19)`。原纪录矩阵达到等号；本结论不分类所有等号矩阵。数学证明依赖报告中列明的结构引理、精确有限计算和 Todorov–Bogdanova (2020) 的已发表码界 B3(15,10)=12；它尚非证明助手形式化。

## 从哪里开始

先读 `report/order15_mu3_complete_proof_38_of_38.pdf` 或其 Markdown 源文件 `report/FINAL_PROOF_zh.md`。最后一个 Q105 候选的排除已经简化为 Eisenstein 范数模 4 的矛盾，以及两对行各需 8 个不同位置却互不重叠的 16>15 矛盾。Q99 的最后纯色配对全部由精确谱迹不等式排除。

最终机器可读定理与每层状态在 `results/final_ledger.json`；完整执行记录在 `results/final_replay_summary.json` 和 `results/final_replay_logs/`。`last_three/` 保留全部新程序、规范 Gram 候选、完整允许列列表、193 个精确湮灭向量和全部 16 个 Q99 同谱配对证书。

## 重新执行

需要 Python 3.11 或更高版本，版本兼容性以四个依赖的安装要求为准。本次实际使用 Python 3.13.5。建议创建独立虚拟环境，安装本包固定依赖后执行：

```bash
python -m pip install -r requirements.txt
python run_final.py
```

已有 uv 的环境也可使用：

```bash
uv venv --python 3.13
uv pip install -r requirements.txt
uv run --no-project python run_final.py
```

默认入口从头重建原审计、全部前阶段目录、全部新增候选和两套完整的列方程枚举。每套新列算法检查 194 × 3^14 = 927895986 个赋值；两套算法逐候选的完整输出必须一致。程序还检查精确逆矩阵、行列式、能量、颜色同余和湮灭向量。所有阶段通过且无余项时，新的入口返回 0。不得使用 `python -O`。

`python run_final.py --cached-baseline` 只复用四个经哈希核验的旧目录；全部新增计算、独立直接列枚举和原审计仍会执行。如果之前完整重建改变了目录中的耗时字段，旧的缓存快照哈希可能不同；此时直接运行无缓存的默认入口，不要为了通过检查而随意改哈希。

## 历史文件与来源

`run_progress.py`、`results/progress_ledger.json`、`report/PROOF_zh.md` 及其 35 层 PDF 保持历史语义。它们仍报告前一阶段的三个缺口；新结论由 `run_final.py` 及 `final_ledger.json` 给出。原 19 层审计完整保存在 `order15_mu3_unified_audit/`，原始 ZIP 也保存在 `provenance/`。

`provenance/source_sha256.json` 固定旧证明源文件，`provenance/final_source_sha256.json` 固定本次新增证明程序。交付快照中的文件哈希另列于 `provenance/final_snapshot_sha256.txt`；重新执行会改写耗时和日志，因此交付快照哈希用于核验下载材料，不用于阻止正常复演。

所有孤立点在必要分支中均被保留。此次没有使用原仓库中失效的行列颜色同向轨道引理，也没有把旧的 CLOSED 标签当作证明。

## 引文勘误

所用外部码界的正确文献是 Todor Todorov and Galina Bogdanova, *Ternary equidistant codes of length 11≤n≤15*, Journal of Mathematical and Computational Science 10 (2020), No. 6, 2713–2721, DOI 10.28919/jmcs/4964。Proposition 14 与 Table 1 给出 B3(15,10)=12。旧 35 层报告中的期号和页码有误，本续篇明确更正；码界数值不变，本次没有重新执行该外部码分类。

本包在本地工作容器中生成；没有向远端 GitHub 仓库写入或推送任何修改。
