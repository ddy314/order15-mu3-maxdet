# Exact continuation certificates

Run `python verify_all.py` from the repository root on a clean checkout. This first generates the original proof catalogues. Then `python verify_continuation.py` can replay just the continuation.

The source uses Python 3.11, the pinned original dependencies and a C++17 compiler. Compiled executables live in temporary directories; no binaries are committed. Assertions must stay enabled. Certificate JSON files are regenerated rather than used as caches, and every input hash is recorded in `data/revision/replay_summary.json`.

| Program | Coverage |
| --- | --- |
| `orthogonality/verify.cpp` | Five complete balanced-row graphs, exact maximum cliques, explicit nine-row witnesses, all simple graphs through six vertices as controls |
| `energy/verify_pure_energy.py` | All 676 labelled states and 373 component spectra through eleven pure energy units; integer norm and spectral-rectangle exclusions |
| `energy/verify_minimum_support.py` | All 24 majority supports at energy 105; all leaf-neighbor frame inequalities |
| `energy/audit_equality.py` | Closed determinant threshold, all 5570 pure Q108 assemblies, high pure and mixed bounds; records original input hashes |
| `energy/verify_minimum_phases.py` | Produces 34 exclusion certificates and independently checks full column coverage, colorings, and all high-color neighborhoods |
| `energy/check_phase_certificates.cpp` | Independent scaled-inverse multiplication and generic quadratic scan; rejects missing candidate columns |
| `factorization/verify_direct.py` | All 5,806,080 row transformations; all 81 direct class comparisons and exact stabilizer orders |
| `sixth_roots/verify_gap.py` | New Q9/Q10 exclusions: 264 leafless supports, 343158 labelled states; norm and rational frame/rectangle obstructions |
| `../tools/verify_frame_controls.py` | 72 one-entry repairs, 2520 two-entry changes, 81 pure tensor controls for the general energy theorems |
| `../tools/verify_fano_completions.py` | All 792 column count types, 675 potential patterns, 141 unique completions, and exhaustive quadrangle control |
| `../tools/verify_continuation.py` | Exact witnesses, two independent determinant representations, 168 Fano actions and nine orbits |

The analytic reductions in the paper establish why these finite families are complete. In particular, pure isolates are excluded by the sparse-core proof, not a determinant-above-threshold filter. The final energy/equality relation has separate forward (column certificates) and reverse (inclusive determinant audit) arguments.

The reconstructed coloring certificates differ from those mentioned in the interrupted chat: the new producer uses more colors and explicitly certifies every higher-color vertex's entire neighborhood. The checker proves the same absence of fifteen-cliques. No old certificate hash is claimed for these new files.
