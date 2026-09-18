# 过程时间线（与指定 Q-shell 直接相关）

## 2026-09-13：Q141 / Q126 / Q144 / Q135 / Q90 / Q87

- Q141 verifier replay 成功，关闭全部 size-13 allocation 与四类 no-isolate 99-energy size-14 support。
- Q126 的 no-K5 瓶颈改用 K4-component tradeoff；随后把一个 12-edge / 10-vertex / weight-3 + K4 结构归约为六叶星附着到同一 K4 顶点，并以 exact PD 检查得到 λmax<27。
- Q144 的 generic energy-54 estimate 失败，随后改用 exact internal enumeration；失败估计没有进入最终证明。
- Q135 all-same 的粗 cap 通过 K5/K6 component split 与 forced-isolate coupling 收紧。
- Q90 顶层 phase class 先尝试较重 cardinality model，随后改成线性 GH=HK/intertwining 模型。18 个 raw Gram → 2 monomial orbits；4 个 ordered orbit pairs 全 UNSAT。
- 同一个顶值出现在 Q87 (13,2)；24 个 raw Gram → 2 monomial orbits；4 个 ordered orbit pairs全 UNSAT。
- energy-9 cap 因此下降到 289220156718750000，已经低于当时的 Q150 trace envelope。

## 2026-09-13 后续：Q141 深化

- 全局瓶颈曾转移到 Q141、size-13、internal-energy-36/minimal-cross。
- 对 K2+K3 与 paw 做 componentwise exact phase enumeration：
  - K2+K3：95,659,380 个压缩 assignment；
  - paw：12,990,780 个压缩 assignment。
- 对应 maxima 约 248.4274×10^15 和 246.9806×10^15，均远低于 Q153 envelope。
- 这关闭了当时最后一个可见的 pre-Q144 Q141 障碍。

## 2026-09-14：pre-Q144 精确清理

- 继续用 exact K4/K5/triangle catalogues 替代共享 coarse cap。
- 10-vertex unicyclic triangle 299-type exact layer最大值：
  246462947301769920。
- K5+leaf exact maximum：
  272942383546368000。
- 最后只剩 Q135 coarse spectral cap；weighted Motzkin–Straus 给出
  λmax ≤ 15+sqrt(202.5) < 29.24，
  使用 exact rational cap 后关闭。
- 得到 certified pre-Q144 maximum：
  287687522593613600，
  严格小于 Q153 envelope：
  287710230040327200。

## 2026-09-15：最终 ledger

- 指定的 15 个 Q-shell 全部记录为 CLOSED。
- 新增结构化 12+1 defect、energy-18 near-equidistant，以及 Q99 local eigenspace obstruction，减少对旧 solver-heavy 路线的依赖。
- Q99 独立 Python 证书保留并可直接执行。
