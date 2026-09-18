# 远端 order15 证明声明审查

审查对象：原仓库 main `827672a2c60dad7064c7955d9c460261d866d25e`；未合并研究分支 `007af7d8f7f644396aee27dd87c68c8f4169efd5` 的六文件差异另存于 `provenance/research-reaggregate.patch`。本地交付包为可信基线，未修改原 ZIP、PDF、CSV。

## 结论

**不能采信“全部做完”。** 原 main 的 `order15_mu3/reports/CURRENT_STATUS.md` 本身明确写着 “Exact maximality proved: no”。它声称完成的是 Q=153,159,162,168 尾部，不是全局最优性。其尾部推理还依赖未成立的同方向颜色轨道推断；因此本次不把尾部测试通过升级为完整的 Q≤150 定理。本地独立复核的 19 个未完整验证壳层保持原状态。

远端材料完整保留在 `historical_remote/`，作为研究来源与待修复资料，不覆盖 `independent_audit/` 的结论。没有发现最大行列式命题的反例；本次发现的是证明步骤与覆盖范围的问题。

## 1. 关键共轭方向问题：失效/需补证

定位：`historical_remote/order15_mu3/scripts/verify_q78_color_orbit_obstruction.py` 的 `certificate()`，关于两侧 singleton-minus-large 差相同、可把两侧交叉项归入同一个 δ 轨道的段落。其程序只检验人为给定的 δ 与共轭 δ 两轨道不相交，没有从 H 推导两侧实际轨道。

令 H_ij=ω^e_ij，行颜色 r_i=Σ_j e_ij mod3，列颜色 c_j=Σ_i e_ij mod3，δ=1−ω。由于 ω^t≡1+t(ω−1) mod3，在 Eisenstein 环模理想 (3) 下：

- G=HH* 满足 G_ik≡−(r_i−r_k)δ mod3；
- K=H*H 满足 K_jl≡−(c_l−c_j)δ mod3。

当两侧颜色型均为 (14,1)，总指数和给出相同的 singleton-minus-large 差 d∈{1,2}。但是两侧“large → singleton”交叉项分别为 dδ 和 −dδ mod3，**方向相反**。乘三次单位根不改变 δ 的模3余类，因为 δ(ω−1)=−δ²=3ω。

在范数3的情形，若行侧交叉项在 δ·μ₃，列侧交叉项应在 −δ·μ₃。代入 GH=HK 时列侧还要取共轭，实际比较为 δ·μ₃ 与 −conj(δ)·μ₃。两者相同，因为 −conj(δ)=δω²。因此原程序的 δ·μ₃ 与 conj(δ)·μ₃ 不相交，不能提供所声称的矛盾。

独立精确核对：

```bash
uv run python scripts/check_remote_orbit_orientation.py
```

该程序同时生成一个可逆 15×15 三次单位根矩阵：行列颜色型均为 (14,1)，有内部孤立行，交叉余类确实相反；完整指数矩阵与精确行列式写入 `provenance/remote_orbit_orientation.json`。其 Q=1347、|det H|²=4421773530027，**不在低能量范围，也不是纪录反例**；其作用仅是核对颜色方向并否定不带额外条件的孤立行排除。低能量情况下是否可由其他证明排除，必须另证。

受影响的引用至少包括 `verify_q123_boundary.py`、`verify_q132_boundary.py`、`verify_q141_boundary.py`、`verify_q150_boundary.py` 以及 `verify_tail_structural_reduction.py::_size14_rows_q159()`。后者删除有内部孤立点的分支，`verify_tail_final_closure.py::q159_size14_certificate()` 也只枚举无孤立点配置。尾部五个测试虽通过，却未检查这个前提；不能据此接受整个 Q=159 或 Q≤150 定理。

## 2. Q150 比较目标：有限证书，未达到纪录级关闭

`verify_q150_boundary.py` 的文档和不等式都比较旧 Q153 迹包络。对其六组 (internal,cross,lambda-cap)，本次用精确有理数重新计算 stationary_upper_dimension(14,internal)·(15−cross/cap)，确认五组仍高于 B=277868041444786176。

```bash
uv run python scripts/compare_remote_bounds.py
```

详细分子、分母和与 B 的比值见 `provenance/remote_bound_comparison.json`。这证明旧界不足，不证明存在超过纪录的矩阵；也不否认本地独立复核已用更强方法关闭无孤立点子分支。

## 3. 实际执行结果

| 核验 | 实际结果 | 可支持的范围 |
|---|---|---|
| 本地八项独立复核 `independent_audit/run_all.py --audit-only` | 全部通过 | 19 壳层关闭、19 壳层未完整验证；不改变原状态 |
| 远端 `verify_benchmark` | 精确通过，det=(604661760,241864704)，范数 B | 与本地纪录及225个指数完全一致 |
| 远端算术、搜索、迹界、颜色界、SAT模型构造等指定测试 | 15 passed | 这些测试的实际断言与有限检查 |
| 远端尾部结构及最终关闭测试 | 5 passed，18.75秒，1个 NetworkX 警告 | 已编码计算重放成功；不验证上述遗漏前提 |
| 直接执行远端尾部 CLI | 非零退出，TypeError: Fraction is not JSON serializable | certificate()计算及断言走完，但输出序列化失败；不能报告 CLI 成功 |
| 全量历史 pytest | 在汇总证书计算期间主动中断 | 未完成，不能称全量通过 |
| 移除 order15 后的 order55 pytest | 4 passed | 剩余 Python 测试；未重跑大型原生55证书 |
| 迁移哈希核验 | 186 个远端原文件、3个本地交付文件、固定 ZIP 内容全部一致 | 来源保全；不是数学证明 |

15个指定测试的命令为 `uv run pytest -q historical_remote/tests/test_mu3.py historical_remote/tests/test_mu3_search.py historical_remote/tests/test_mu3_circulant.py historical_remote/tests/test_mu3_trace.py historical_remote/tests/test_mu3_color_bounds.py historical_remote/tests/test_mu3_staged.py historical_remote/tests/test_mu3_support_sat.py historical_remote/tests/test_mu3_sat_q78_paw.py historical_remote/tests/test_mu3_sat_q90_top.py`。尾部命令为 `uv run pytest -q historical_remote/tests/test_mu3_tail_final_closure.py historical_remote/tests/test_mu3_tail_structural_reduction.py`。

未合并研究分支仅保存差异，未重跑其全量汇总；它仍使用同一尾部结论，未独立补足方向和孤立分支证明。本次未完成对每个远端细粒度引理的独立证明，不把未审查部分标为错误或已完成。

## 后续数学工作

首先修复或替代颜色轨道引理，并显式覆盖 Q159 等壳层的孤立内部支撑分支；再证明完整分支覆盖和纪录级界。不得通过修改返回的 closed_shells、跳过分支或固定常数来“关闭”证明。新仓库以本地独立复核包作为当前基线，历史关闭声明保留供审查。
