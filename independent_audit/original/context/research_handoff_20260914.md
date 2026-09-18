# Order-15 μ3 Maximal Determinant — Research Handoff

## 0. Purpose

This handoff is for continuing the exact maximal-determinant proof at order 15 over the cube roots of unity

\[
\mu_3=\{1,\omega,\omega^2\},\qquad 1+\omega+\omega^2=0.
\]

The project is **not solved yet**. The current task is to continue from the existing certified checkpoint and close the remaining energy shells, not to restart literature review, baseline search, or heuristic optimization.

---

## 1. Problem

Let \(H\in\mu_3^{15\times15}\). Define

\[
M_3(15)=\max_H |\det H|,\qquad D_3(15)=M_3(15)^2.
\]

The known explicit lower-bound matrix has

\[
D_0
=
2^{22}3^{20}\cdot 19
=
277868041444786176,
\]

hence

\[
M_0=2^{11}3^{10}\sqrt{19}.
\]

Its ratio to the Hadamard bound is about

\[
\frac{M_0}{15^{15/2}}\approx 0.7965900126.
\]

The research goal is to prove this record is globally maximal, or otherwise sharpen the rigorous interval as far as possible.

---

## 2. Current proof architecture

The proof works with the row Gram matrix

\[
G=HH^*,
\]

and the off-diagonal energy

\[
Q=\sum_{i<j}|G_{ij}|^2.
\]

Major ingredients already implemented and used:

- dephasing and exact Eisenstein-integer arithmetic;
- allowed order-15 Gram entry catalogue;
- trace / spectral stability bounds;
- row-color congruence from ternary exponent sums;
- color-partition reduction;
- support-graph decomposition;
- exact component determinant factorization;
- Schur-complement upper bounds;
- kernel / nullity / rank obstructions;
- two-sided row/column Gram compatibility;
- \(GH=HK\) / intertwining obstructions;
- finite phase enumeration modulo switching / monomial symmetry;
- exact SAT / Z3 UNSAT certificates where needed;
- Eisenstein-norm and \(3^{14}\)-divisibility arithmetic sieve;
- regression tests and replayable certificate scripts.

Do **not** discard this framework or rebuild from scratch.

---

## 3. Important proved structural facts

### 3.1 Color congruence

For row exponent sums \(s_i\in\mathbb Z/3\mathbb Z\),

\[
\frac{|G_{ij}|^2}{3}\equiv (s_i-s_j)^2\pmod 3.
\]

Consequences:

- orthogonal rows have the same color;
- norm-3 Gram pairs have different colors;
- global color partitions are heavily restricted.

### 3.2 Friendship-family result

The known record has friendship-graph support \(F_7\). The entire \(F_7\) support family was exhaustively certified, and the known record is maximal inside that family.

### 3.3 Low-energy elimination machinery

Many low-energy shells have been fully eliminated using combinations of:

- exact support enumeration;
- Eisenstein norm arithmetic;
- kernel/rank constraints;
- Schur bounds;
- exact intertwining / realizability obstructions.

The proof frontier has already been pushed far beyond the initial \(Q\le168\) necessary-condition regime.

---

## 4. Latest exact state from the last run

The latest run completed the analytic / finite cleanup needed to push every known **pre-\(Q=144\)** obstruction below the \(Q=153\) envelope.

The key comparison is:

\[
U_{\text{pre-}Q144}
=
287687522593613600
\]

versus

\[
U_{Q=153}
=
287710230040327200.
\]

Thus the former is strictly smaller by about

\[
2.27\times 10^{13}.
\]

Important: at the time of handoff, the **full unmocked aggregate replay was still running**. It had not reported an assertion failure, but it had not yet emitted the final PASS line. Therefore:

> First finish / rerun the full aggregate replay before promoting the \(Q=153\) frontier to fully certified status.

Primary aggregate command used in the project:

```powershell
& 'E:\maximal determinant\.venv\Scripts\python.exe' -m order15_mu3.scripts.verify_aggregate_upper
```

If the environment path differs, use the repository venv and the same module.

---

## 5. Latest exact finite certificates added

The latest cleanup includes, among others:

- pure \(K_4\) component-factorized bounds;
- weighted \(K_4\) forest cases;
- pure \(K_4\)+forest exact phase catalogues;
- weighted \(K_5\), pure \(K_5\), and \(K_5+\)leaf branches;
- isolated-\(K_4\) exact phases;
- \(K_6-e\) exact transform;
- active triangle components on 4–10 vertices;
- unicyclic / bicyclic triangle catalogues;
- exact \(Q=87\) / \(Q=90\) Gram-orbit intertwining UNSAT checks;
- refined \(Q=135\) weighted Motzkin–Straus spectral cap.

Representative exact maxima from the latest run include:

\[
K_6-e:\quad
286978140000000000,
\]

\[
K_5+\text{leaf}:\quad
272942383546368000,
\]

\[
\text{10-vertex unicyclic triangle layer}:\quad
246462947301769920,
\]

with many other component maxima already stored in certificate scripts and regression invariants.

Use the repository as the source of truth for the complete list.

---

## 6. Remaining frontier

Assuming the full replay passes, the proof should proceed from the \(Q=153\) boundary.

The remaining admissible energy-shell map is expected to be very short:

\[
Q=153,\ 159,\ 162,\ 168,
\]

with the congruence-forbidden intermediate shells skipped.

The working endgame is:

\[
Q=153
\rightarrow159
\rightarrow162
\rightarrow168
\rightarrow
\text{trace cutoff}.
\]

At sufficiently high \(Q\) (expected by \(Q\ge171\) under the current stationary trace machinery), the continuous trace upper bound itself drops below the known record, so the tail closes automatically.

Treat this shell map as the continuation target, but **recompute / verify it from the repository formulas before using it as a theorem statement**.

---

## 7. Recommended continuation strategy

### Step 1 — Replay first

Run:

```powershell
python -m order15_mu3.scripts.verify_aggregate_upper
```

and the full test suite.

Do not rely on copied numerical summaries if the replay fails.

### Step 2 — Identify the exact \(Q=153\) bottleneck

Recompute all color partitions and exact branch maxima at \(Q=153\). Rank them.

Do not broaden the search. Attack only the branch(es) that actually attain the global envelope.

### Step 3 — Prefer structural closure over brute-force SAT

Priority order:

1. exact component factorization;
2. Schur / Sylvester / Perron refinements;
3. moment constraints, especially
   \[
   \operatorname{tr}(E)=0,\quad
   \operatorname{tr}(E^2)=Q,\quad
   \operatorname{tr}(E^3)
   \]
   where support structure fixes or constrains the third moment;
4. kernel / nullity / two-sided rank obstruction;
5. exact \(GH=HK\) intertwining;
6. finite phase enumeration modulo switching / monomial symmetry;
7. SAT / Z3 only after the phase/orbit class is small.

Never treat solver timeout as UNSAT.

### Step 4 — After closing \(Q=153\)

Immediately promote the global bound to the next admissible shell and repeat. Do not stop at an intermediate milestone.

### Step 5 — If exact maximality is achieved

Continue with:

- independent replay;
- clean-room verification of the record determinant;
- minimal trusted computing base;
- proof/certificate compression;
- classification of maximizers up to monomial equivalence if feasible;
- automorphism group / Gram spectrum;
- near-maximizer gap;
- structural explanation of the \(\sim80\%\) phenomenon;
- possible generalization to globally norm-obstructed \(\mu_3\) orders.

---

## 8. Things that should NOT be redone

Do not spend time on:

- generic heuristic matrix search;
- another broad literature survey unless a precise new lemma needs citation;
- re-deriving the known benchmark;
- redoing the circulant scan;
- repeating old SAT formulations that already timed out;
- recomputing expensive exact catalogues if the certificate scripts already contain and replay them;
- replacing exact finite certificates by weaker floating-point bounds.

If a previous result is needed, replay the existing script rather than rediscovering it.

---

## 9. Provenance / rigor rules

Keep these distinctions explicit:

- **proved / replayable exact result**
- **finite exhaustive computation**
- **SAT/Z3 UNSAT certificate**
- **continuous analytic upper bound**
- **heuristic observation**
- **timeout / unknown**

Never upgrade the last two categories into proof.

Any new numerical upper must be passed through:

1. exact rational arithmetic where possible;
2. \(3^{14}\)-divisibility;
3. Eisenstein norm admissibility;
4. row/column realizability constraints when relevant.

---

## 10. Immediate research target

The next run should aim to answer:

> Can the complete \(Q=153\) shell be pushed below the current record, or at least below the next admissible-shell trace envelope, using the existing component/moment/intertwining machinery?

If yes, continue immediately to \(Q=159\), then \(162\), then \(168\).

Do not voluntarily stop while the execution budget remains and the exact problem is unresolved.

---

# Prompt to give the next Codex / research agent

```text
CONTINUE THE EXISTING ORDER-15 μ3 MAXIMAL-DETERMINANT PROOF FROM THE CURRENT REPOSITORY CHECKPOINT.

Do not restart the project, do not redo generic literature review, do not rerun heuristic matrix search, and do not repeat already-certified low-energy work except as replay/verification.

Problem:
For H ∈ μ3^(15×15), μ3={1,ω,ω²}, determine the exact maximum of |det H|.

Known explicit record:
D0 = |det H0|² = 2^22 * 3^20 * 19 = 277868041444786176,
so M0 = 2^11 * 3^10 * sqrt(19).

The repository already contains:
- exact Eisenstein arithmetic,
- trace/stability bounds,
- color congruence,
- support-graph enumeration,
- kernel/rank obstructions,
- Schur/component bounds,
- exact phase catalogues,
- GH=HK intertwining tests,
- SAT/Z3 exact realizability certificates,
- arithmetic norm/divisibility sieves,
- regression tests and aggregate verifiers.

CURRENT CHECKPOINT:
The latest analytic/finite cleanup gives

pre-Q144 maximum = 287687522593613600
Q153 envelope      = 287710230040327200

so all currently known pre-Q144 branches are below Q153.

HOWEVER, the last run ended while the FULL UNMOCKED aggregate replay was still running.
Therefore your FIRST TASK is:

1. Run the full exact aggregate replay:
   python -m order15_mu3.scripts.verify_aggregate_upper
2. Run the full test suite.
3. If either fails, repair the proof/certificate and do not promote the bound until replay passes.
4. If replay passes, treat Q=153 as the active frontier.

Then continue the proof.

EXPECTED REMAINING ENERGY-SHELL MAP:
Q = 153, 159, 162, 168,
with congruence-forbidden shells skipped, and the stationary trace tail expected to fall below the known record by Q >= 171.

VERIFY this shell map from the repository formulas before relying on it.

For each surviving shell:
A. enumerate the allowed color partitions;
B. rank exact branch upper bounds;
C. attack only the branch(es) attaining the global maximum;
D. prefer, in this order:
   - component factorization,
   - exact Schur/Sylvester/Perron bounds,
   - trace-moment constraints including tr(E^3),
   - kernel/nullity/two-sided rank,
   - GH=HK intertwining,
   - finite switching/monomial phase enumeration,
   - exact SAT/Z3 only after the class is small;
E. apply 3^14 divisibility and Eisenstein-norm arithmetic;
F. update replayable certificate scripts and regression tests;
G. immediately promote to the next admissible shell after closure.

Never interpret timeout as UNSAT.
Never replace a replayable exact certificate by a weaker numerical estimate.
Never stop merely because a new upper bound is obtained.

PRIMARY GOAL:
Prove the known record is globally maximal:
D3(15) = 277868041444786176.

If exact maximality is reached, continue with:
- independent clean replay,
- proof compression,
- minimal trusted computing base,
- maximizer classification up to monomial equivalence if feasible,
- automorphism/Gram-spectrum analysis,
- near-maximizer gap,
- structural explanation of the ~80% Hadamard-bound realization,
- generalization to globally norm-obstructed μ3 orders.

Maintain a concise CURRENT_STATUS / handoff file after every major promotion.
Do not voluntarily terminate while execution budget remains and a mathematically meaningful next step exists.
```
