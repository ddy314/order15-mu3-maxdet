# Final supplement: maximal determinant at order 15 over third roots of unity

Let ω²+ω+1=0, and let H have order 15 and all entries in {1,ω,ω²}. Combining the retained 35-shell proof with the three new shell exclusions in this package gives the computer-assisted theorem

**max |det H|² = 277868041444786176 = 2²² · 3²⁰ · 19.**

Thus max |det H| = 120932352√19. The previously recorded matrix attains this value; its 225 exponents and three exact determinant computations are retained. No classification or uniqueness statement about equality matrices is asserted.

The proof retains the published theorem B₃(15,10)=12 of Todorov and Bogdanova and the structural lemmas in the original audit and the 35-shell supplement. The new report is `FINAL_PROOF_zh.md` and its PDF. This is a computer-assisted mathematical proof, not a proof-assistant formalization.

## 1. Exact column feasibility

For any invertible factor H with HH*=G, necessarily H*G⁻¹H=I. Every column v therefore satisfies v*G⁻¹v=1. Fixing its first coordinate to 1 loses no columns up to a common third root, leaving exactly 3¹⁴ possible columns per G. Two distinct exact algorithms, incremental and direct, exhaust this whole set for every candidate. Both give identical complete lists.

The three shells share a size-13 residual with internal energy 18, cross energy 78, and outside squared modulus 0, 3, or 9. Its internal support is P₃ plus 10 isolates or two disjoint edges plus 9 isolates. Normalization uses only third-root row switches so that the allowed factor alphabet is preserved. Exhausting the remaining core phases, cross ratios, isolate multiplicities, and outside Gram entry produces 142 norm-admissible representatives above the record: 64 for Q96, 67 for Q99, and 11 for Q105. Of these, 141 have deficient allowed-column span.

## 2. The final Q105 candidate has an analytic contradiction

The remaining G has majority block A=diag(C,C,15I₉), where C has diagonal 15 and off-diagonal 3. Both hub cross vectors are z=(1−ω)1₁₃, and the hub mutual entry is 3. Its determinant is 281238577200000000.

Write a possible column as (a,b,c,d,p₁,…,p₉,u,v), put δ=1−ω, s=a+b+c+d, P=Σpⱼ and k=s/18+P/15. The inverse quadratic form is exactly

q = x*A⁻¹x + (15/392)|u+v−2δ̄k|² + |u−v|²/24.

If both a≠b and c≠d, the first term is 163/180. If u≠v, q≥371/360>1. If u=v, the equality q=1 instead forces

**N(30δu−5s−6P)=1666.**

This is impossible because Eisenstein norms have residues only 0,1,3 modulo 4, whereas 1666 has residue 2. Consequently every permissible column has a=b or c=d. Each of the two row pairs has Gram inner product 3 and hence exactly 8 disagreements across the 15 columns. Those two disagreement sets must be disjoint, requiring 16 distinct positions among 15. This proves non-realizability without relying on the column enumeration for the final column restriction. The enumeration of 5887 normalized allowed columns independently confirms it.

The inherited double-(14,1) bounds cover the other Q105 color possibility, including the equality branch at the record.

## 3. The final Q96 mixed supports

Once size-13 is excluded, both Gram color types must be (14,1). No-isolate supports would already cost energy at least 105. The inherited two-sided rank constraints leave exactly three rank signatures. A fresh support/phase reconstruction identifies them as balanced P₃+C₄ with 7 isolates, C₄ with 10 isolates, or K₂,₃ with 9 isolates. All excess cross-energy placements, including placements on isolated vertices, are retained. The three normalized enumerations have 7776, 4536 and 1296 states; only the first family leaves candidates above the record and norm sieve, namely 52 Grams. All 52 have allowed-column span of dimension at most 9.

## 4. The final Q99 pure-Gram pairs

Both Grams must now be all-same-color. With T=(G−15I)/3 and S=(K−15I)/3, TH=HS. Two-sided kernel projection eliminates all isolated-support candidates from the complete component catalogue. The no-isolate case has 632 assemblies, with 14 above-record norm-admissible component-spectral assemblies and 16 unordered equal-whole-spectrum pairs. All component decompositions are retained, rather than choosing a single representative per whole spectrum.

For unions I,J of complete components, X=H[I,J] satisfies T[I]X=XS[J] and XX*≤G[I]. If the greatest common divisor of the two component characteristic polynomials is

p(x)=xᵈ+c₁xᵈ⁻¹+…,

then a necessary condition is

**|I||J| = ‖X‖²_F ≤ 15d−3c₁.**

Indeed only common eigenspaces can be coupled, with rank bounded by the smaller eigenvalue multiplicity and each singular value squared bounded by 15+3λ. Summing these separate eigenvalue contributions gives the stated integer bound. Every one of the 16 pairs has a violating rectangle. The three pairs surviving a weaker bound give respectively 40≤30, 20≤18, and 36≤30.

## 5. Reproduction and scope

Run `python run_final.py` from the package root after installing `requirements.txt`. This reruns the original audit, all previous catalogues, every new candidate generation, both complete column algorithms on all 194 candidate Grams, and all exact certificates. There are 193 explicit Eisenstein annihilators for the deficient-span cases; the last Gram is excluded analytically as above. Both algorithms together execute 1855791972 column assignments. The new aggregator verifies all coverage joins before reporting 38/38 and a verified infinite tail beginning at Q171.

Correct bibliography for the external theorem: T. Todorov and G. Bogdanova, “Ternary equidistant codes of length 11≤n≤15,” *Journal of Mathematical and Computational Science* 10 (2020), No. 6, pp. 2713–2721, DOI 10.28919/jmcs/4964, Proposition 14 and Table 1. The issue and page range in the historical 35-shell report were miscited; the numerical theorem used is unchanged. The external code classification was not rerun here.
