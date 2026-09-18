# Source and reconstruction notes

The uploaded archive is retained file-for-file in `original/`. It explicitly describes itself as a recovery of conversation materials, not the complete original working repository. Historical `CLOSED` labels are not adopted as verification evidence.

The repository examined through the GitHub connector was `dongxuelian2/order55-binary-circulant-maxdet`. The available refs in the examined snapshot were main `827672a2c60dad7064c7955d9c460261d866d25e` and research/reaggregate-after-tail-closure `007af7d8f7f644396aee27dd87c68c8f4169efd5`. This delivery does not contain a complete checkout. Binary archive retrieval into the working container failed; no newer full local workspace was recovered.

The benchmark exponent matrix was recovered from `order15_mu3/data/literature_matrices/nunez_ponasso_m15.json`, at the latter commit; its repository blob SHA is `b140314ff2f25696ad7426536d41c33fa5ddb920`. The local JSON has added provenance metadata and therefore has its own, different file hash.

The inspected `order15_mu3/scripts/verify_q150_boundary.py` begins with the description “Exact color/support certificate placing Q=150 below Q=153.” Its comparison is `upper < q153`. The six size14 (internal,cross,lambda-cap) triples are (63,87,25), (72,78,127/5), (81,69,261/10), (90,60,27), (99,51,138/5), (108,42,563/20). The new `verify_component_bounds.py` recomputes these expressions, shows that five still exceed B, and independently replaces the entire no-isolate size14 branch with stronger component/tree bounds. The size13 and isolate branches remain separate obligations.

The inspected research-branch aggregate program combines earlier, Q144 and Q150 bounds and then applies arithmetic sieves. Its assertions do not establish that the returned global upper bound equals B. Historical tail-closure code depends on sparse support and spectral reductions which were not completely recovered and replayed in this delivery. A successful historical CI result, where encountered, has not been substituted for a full replay.

All files in `programs/` are new reconstructions written for this audit. They do not purport to be byte-for-byte copies of missing repository modules. The original Q99 program is preserved and is used only for comparing its explicitly enumerated local sets against independently derived complete local solution sets. No unavailable repository certificate is mocked.

External theorem used: Todor Todorov and Galina Bogdanova, “Ternary equidistant codes of length 11 through 15”, J. Math. Comput. Sci. 10 (2020), 2713–2721, DOI 10.28919/jmcs/4964. Proposition 14 on printed p.2717 and Table1 on p.2720 give B_3(15,10)=12. The PDF was read through the web tool, including the table image; its classification computation is not rerun here. Source: https://scik.org/index.php/jmcs/article/viewFile/4964/2448.
