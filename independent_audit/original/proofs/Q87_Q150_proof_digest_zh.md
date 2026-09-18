# 指定 Q-shell 证明整理（Q=87…150）

## 1. 问题与基准

研究对象是
\[
H\in \mu_3^{15\times 15},\qquad \mu_3=\{1,\omega,\omega^2\},
\]
行 Gram 矩阵
\[
G=HH^*,
\]
以及非对角能量
\[
Q=\sum_{i<j}|G_{ij}|^2.
\]

已知显式基准矩阵满足
\[
D_0=|\det H_0|^2
=277868041444786176
=2^{22}3^{20}19.
\]

最终 ledger 使用的全局约化为：

1. 精确 trace/variance 优化把严格反例限制在 \(Q\le 168\)；\(Q\ge171\) 已低于 \(D_0\)。
2. 行颜色由指数和 mod 3 定义，并满足
   \[
   |G_{ij}|^2/3\equiv (s_i-s_j)^2\pmod 3.
   \]
   所以同色内积范数为 \(0\bmod 9\)，异色为 \(3\bmod 9\)，进而
   \[
   Q\equiv0\ \text{或}\ 6\pmod 9.
   \]
3. 严格反例只需考虑颜色分拆
   \[
   (15,0,0),\ (14,1,0),\ (13,2,0),\ (13,1,1).
   \]
4. 反复使用精确 ternary code bound
   \[
   B_3(15,10)=12,
   \]
   因而 13 个及以上同色行不可能完全两两正交。

下面只整理本次指定的 15 个 shell。

---

## 2. 两个反复使用的结构引理

### 2.1 12+1 defect 引理

若 12 个长度 15 的 \(\mu_3\) 行两两正交，第 13 行只与其中一个出现单个非零同色内积，则 rank-3 complement tight frame 的正交图满足 affine-ADE 型约束 \(Om=2m\)。

对 affine-\(A_2\) 基底使用 Parseval，并对其余 affine-ADE 局部类型使用 rank-3 障碍，可以排除这里需要的所有可达同色 defect norm，特别是 norm \(9\)。

因此它统一消掉
\[
Q=87,90,96,99,105,108
\]
中的 size-13、\(e=9\) 分支。

### 2.2 energy-18 near-equidistant 引理

13 行同色块若内部能量为 18，则支撑只能是 \(P_3\) 或 \(2K_2\)。

- \(P_3\)：删掉中间行后得到 12 个互相正交行；affine-ADE + \(3\times3\) contingency 约束只剩局部 \(C_4/C_5/E_6\) 或 \(D_4\)。前者违反 rank 3；\(D_4\) 违反精确 ETF projection identity
  \[
  3\sum M_i^2=(\sum M_i)^2.
  \]
- \(2K_2\)：rank-2 complement 立即排除三类符号；最后的 negative-negative 情形只剩两个精确 \(3\times3\) joint table。每个产生 6579 个候选 completion row；对应精确三元正交图有 2,436,912 条边，且不存在 9-clique（完全整数 difference-count adjacency，无浮点依赖）。

所以相关 size-13 的内部能量 18 分支全部不可能。

---

## 3. 各指定 shell

### Q=87 — CLOSED

最终结构闭合：
- \((13,2)\) 的三种能量分配分别由 orthogonal-class/code bound 或 12+1 defect 引理排除；
- \((14,1)\times(14,1)\) 分支由 isolated-orbit obstruction 排除。

研究过程中还有一条独立的有限计算路线：24 个 raw \(Q=87\) Gram 压缩为两个 monomial orbit，四个有序 orbit pair 的 intertwining 检查全部精确 UNSAT。这条路线保存在过程日志中，但最终 ledger 采用了更结构化的闭合表述。

### Q=90 — CLOSED

- size-13 的关键低能分支由 12+1 defect / affine-ADE 排除；
- all-same 最终只剩唯一
  \[
  C_5^-+5K_2^+
  \]
  orbit，由 eigenspace block / Cauchy contradiction 排除。

过程日志中另有独立计算证据：18 个 raw \(Q=90\) Gram 压缩为两个 monomial orbit，四个有序 orbit pair 全部 exact UNSAT。早期使用过 `z3_q90_top_commuting` 的线性 intertwining 形式；这比原先的 cardinality model 更强、更快。

### Q=96 — CLOSED

- size-13 \(e=9\)：12+1 defect；
- size-13 \(e=18\)：\(P_3/2K_2\) 结构不可能；
- 其余 \((14,1)\times(14,1)\) 支撑均带 isolate，因而被已有 isolated-support 障碍排除。

### Q=99 — CLOSED

size-13：
- \(e=9\) 由 defect 引理排除；
- \(e=18\) 由 \(P_3/2K_2\) 引理排除。

all-same 最终仅剩两个 cospectral orbit：
\[
A=4K_2^-\sqcup K_3^+\sqcup C_4^+,
\]
\[
B=3K_2^-\sqcup K_3^+\sqcup T_6,
\]
共同特征多项式
\[
(x-21)^2(x-18)^4(x-15)^2(x-12)^6(x-9).
\]

15-特征空间只支撑在 \(A\) 的 \(C_4\) 和 \(B\) 的 \(T_6\) leaf-difference 方向，因此局部块必须满足 intertwining/commutant 方程。包内独立程序 `verify_q99_final_orbits_exact.py` 精确枚举得到：

- A/B：144 个 \(\mu_3\)-valued local intertwiners，residual-compatible 为 0；
- A/A：225 个 local commutants，residual-compatible 为 0；
- B/B：729 个 local commutants，residual-compatible 为 0；
- B/A 由 A/B 的共轭转置排除。

本 ZIP 已实际重新运行该脚本，并保存输出 JSON。

### Q=105 — CLOSED

- size-13 \(e=9,e=18\)：结构不可能；
- \(e=27\)：使用精确上界；
- \((14,1)\times(14,1)\) 的唯一无-isolate \(7K_2\) 分支最大值恰好等于 \(D_0\)。

因此该 shell 不能产生严格超过已知纪录的矩阵。

### Q=108 — CLOSED

- all-same：仓库中的 exact theorem；
- size-13 \(e=9\)：defect；
- \(e=18\)：exact two-hub Schur；
- \(e=27\)：inverse-loss bound。

### Q=114 — CLOSED

由两部分合并：
- size-13：exact residual / inverse-loss bounds；
- size-14：componentwise Schur bounds。

所有允许颜色分配均低于严格反例阈值。

### Q=117 — CLOSED

- all-same：exact support theorem；
- size-13：exact inverse-loss / Schur bounds。

### Q=123 — CLOSED

- size-13：exact hub losses；
- size-14：componentwise bounds。

### Q=126 — CLOSED

这是中途最难的旧瓶颈之一。

最终 ledger：
- all-same 由 cubic log majorant + exact weighted Fourier 完整闭合；
- 全部 size-13 allocation 精确低于纪录。

过程中关键修正包括：

1. no-\(K_5\) 谱半径瓶颈用严格的 \(K_4\)-component tradeoff 解决：component energy 与 forced isolates 不能同时取各自最坏值，逐个精确 \((m,r)\) 分支均低于相应 envelope。
2. 一个最后的 12-edge、10-vertex 结构分支含一个 weight-3 edge 与 \(K_4\)。其余六条非-\(K_4\) 边组成 rooted forest 并新增六个顶点；Perron edge-moving 把极值归约到“六片叶全部接在同一 \(K_4\) 顶点”。精确 positive-definiteness 检查得到
   \[
   \lambda_{\max}<27.
   \]
3. 过程中还做了 pure/weighted \(K_4\) forests、\(K_5\)、\(K_5+\)leaf、isolated-\(K_4\)、\(K_6-e\)、4–10 顶点 triangle/unicyclic/bicyclic exact phase catalogues。它们用于把 coarse cap 逐步替换成可重放的有限上界。

### Q=132 — CLOSED

- size-13：exact refinements；
- size-14：no-isolate component enumeration。

### Q=135 — CLOSED

最终闭合依赖：
- cubic log majorant；
- complete high-moment triangle-core enumeration；
- exact integer Eisenstein Fourier；
- size-13 全部低于纪录。

最后一个 coarse spectral cap 在过程阶段通过 weighted Motzkin–Straus 收紧：
\[
\lambda_{\max}\le 15+\sqrt{202.5}<29.24.
\]
使用精确有理 cap 29.24 后，该分支降到下一 envelope 以下。此前还按 spectral component 是否含 \(K_5/K_6\) 分拆，并将 component energy 与 forced isolates 配对，避免不相容的“同时最坏”估计。

### Q=141 — CLOSED

最终 ledger：
- exact size-13 residual；
- size-14 no-isolate support bounds。

过程中的关键可重放节点：
- \(Q=141\) certificate 成功 replay；
- 关闭全部 size-13 allocation 与四类 no-isolate 99-energy size-14 support；
- 对紧的 \(e=36\) 情形，明确分离 determinant-maximizing support 与 spectral-radius-maximizing support，避免把不相容极值拼在一起；
- 后续 componentwise enumeration 对 \(K_2+K_3\) 和 paw 支撑分别做了 95,659,380 与 12,990,780 个压缩 phase assignment，得到约
  \[
  248.4274\times10^{15},\quad 246.9806\times10^{15},
  \]
  均显著低于后续 \(Q=153\) envelope。

### Q=144 — CLOSED

- all-same：exact \(K_6\) boundary / arithmetic；
- size-13：exact Schur / state pairing。

过程里曾出现一次失败的 generic energy-54 estimate；该估计没有被当作证明使用，而是改成 exact internal enumeration 后才闭合。这一点在过程日志中保留，便于审计“失败估计 → 精确替代”的链条。

### Q=150 — CLOSED

- size-13：exact state pairing；
- size-14：componentwise closure；
- 包括边界
  \[
  e=108,\qquad K_2+T_{12}.
  \]

在推进到 \(Q=150\) envelope 的途中，旧的 \(Q=90\) energy-9 phase maximum 一度成为全局瓶颈；对其完成精确 orbit/intertwining UNSAT 后，才允许全局上界继续下降。这个过程解释了为什么 ZIP 中保留 \(Q=87/90\) 的 solver/intertwining 日志，即使最终 ledger 已用更结构化表述收口。

---

## 4. 最终状态

上述 15 个指定 shell：
\[
87,90,96,99,105,108,114,117,123,126,132,135,141,144,150
\]
在最终 record-level ledger 中全部为 **CLOSED**。

完整项目的 ledger 进一步声称所有 \(Q\le168\) 的可容许 shell 均已 record-level closed，而 \(Q\ge171\) 由精确 trace/variance certificate 低于基准，因此得到
\[
\max_{H\in\mu_3^{15\times15}}|\det H|^2
=277868041444786176.
\]

本文件只负责整理指定 shell；完整最终断言请连同 `final_record_level_ledger.md`、仓库中的原始 verifier 以及 clean replay 一起审计。
