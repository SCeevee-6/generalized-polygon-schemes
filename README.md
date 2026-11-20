# generalized-polygon-schemes

Association-scheme computations around generalized hexagons and octagons and the Tits group T = 2F4(2)'. Mostly fusions of the point-distance scheme of a generalized polygon, plus orbital algebras and graph searches for T acting on cosets of its subgroups.

Work in progress, not peer reviewed. `fusion/fusion_note.tex` is the draft of the paper.

## Results

- **Fusions.** Take a generalized 2d-gon of order (s,t), d = 3 or 4, with integers s,t >= 2. Its point-distance scheme has a proper fusion iff s = t, and then the fusion is unique: {1,...,d-1} | {d}. Exact proof from the intersection array (Groebner bases over Q, exact real-root isolation, integer scan) in `fusion/certify_fusions.py` and `fusion/certify_note.md`. The hexagon case is known. For octagons I haven't found it in the literature, but I only did keyword searches.
- **Hexagons, s = t = q.** The fusion is strongly regular with parameters ((q^6-1)/(q-1), q(q+1)(q^2+1), q^3+q^2+q-1, (q+1)(q^2+1)), same as the collinearity graph of Q(6,q).
- **Octagons.** s = t is ruled out by Feit-Higman (2st has to be a square), so the Ree-Tits octagon GO(2,4) has no proper fusion.
- **T on 1755 points, (s,t) = (2,4).** Intersection array {10,8,8,8;1,1,1,5}. The eigenvalues of A_1 are 10, 5, 1, -3, -5 with multiplicities 1, 351, 650, 675, 78. The only vanishing non-principal Krein parameter is q^4_44 = 0. Of the 15 fusions of {1,2,3,4} only the two trivial ones close. The scheme is not Q-polynomial in any ordering (De Bruyn-Vanhove 2015) and the graph is not tight (the Jurisic-Koolen-Terwilliger bound reads 0 >= -20, strict; `tits1755/tightness_check.py`). Full write-up in `docs/tits1755_orbital_algebra.md`.
- **Other actions of T.** On cosets of the 8 maximal subgroups the orbital configurations are commutative for degrees 1600 (twice), 1755, 2304 and 2925, and symmetric except for 2304 (L2(25)). They are non-commutative for A6.2^2 (degree 12480, rank 24, twice) and 5^2:4A4 (degree 14976, rank 25). See `tits1755/coherent_configurations_note.md`.
- **Fusions of the commutative actions** (`tits1755/coset_fusions_note.md`). Degree 1600: the one proper fusion is the known SRG(1600,351,94,72) from Crnkovic-Rukavina-Svob (2020). Degree 2304: three fusions (the T.2 orbital scheme of rank 6, its symmetrisation of rank 6, and the symmetrisation of the T.2 scheme of rank 5), none P- or Q-polynomial. Degree 2925: one fusion, the distance scheme of the dual octagon GO(4,2); P- but not Q-polynomial.
- **Non-commutative actions.** The orbital algebras are 11 C + M_2(C) + M_3(C) (degree 12480) and 9 C + 4 M_2(C) (degree 14976). The centre has dimension 13 in both cases, computed exactly from the structure constants and consistent with the character multiplicities.
- **SRG/DRG search for T beyond Crnkovic-Rukavina-Svob** (`search/search_note.md`). Exhaustive search over all unions of orbitals for every T-action of degree 7020 to 17550 (TomLib classes), including the primitive actions on 12480 and 14976 points. No SRG and no DRG of diameter >= 3 on 12480 or 14976 points, so together with CRS 2020 the primitive case is done. On 7020, 8775, 11700 and the 17550-point actions only complete multipartite SRGs show up. No new graph.
- **Brouwer's DRG tables** (`drgfusion/drgfusion_note.md`). Closed fusions of all tabled arrays of diameter <= 7, with the fused SRGs on at most 1300 vertices looked up in Brouwer's SRG table. Nothing new came out: of the six arrays that looked open, two fail the absolute bound, two are deleted in the BCN additions, and two follow from Koval's October 2026 preprint (conditional).
- **Fusions of the non-commutative T schemes** (`search/fusions/fusion_note.md`). The 12480-point configurations have no proper fusion of rank >= 3. Every fusion of rank >= 3 of the 14976-point configuration refines the rank-16 2F4(2) orbital configuration and is non-commutative.
- **Regular near hexagons and octagons** (`fusion/near_polygon_note.md`, `fusion/near_octagon_note.md`). Diameter 3 is BCN 4.2.17. For regular near octagons (s, t2, t3, t) with s >= 2 the integer parameter sets with a proper fusion are: (s, m, s(m+1), (m+1)(s^2+1)-1) with {1,2}|{3,4} (m = s is B_4(s)/C_4(s)), (s,0,0,s), (2, i, 2i, 2i+1) (only the thin octagon (2,1) and H(4,3) are feasible), and H(4,4). This is proved for ten of the thirteen partitions. For {1,4}|{2,3}, {1,3}|{2,4} and {1,2,4}|{3} one residual surface each is left, with no further integer point for s <= 200, t3 <= 2000. All examples are known.
- **GO(2,4) triple intersection numbers** (`triple/triple_note.md`). With q^4_44 = 0 the LP ranges only reproduce the octagon geometry for one triple type, so nothing on uniqueness.
- **SRGs on open parameter sets from primitive groups** (`search/primitive/primitive_note.md`). All 2164 primitive groups of rank >= 4 on the 430 open degrees v <= 1300 are settled, and none realises an open parameter set in Brouwer's table (as of 1 Nov 2026). All of it is exact now: the 14 groups that used to rest on floating-point MILPs were redone with integer character sums in CP-SAT (12 affine groups) and a fixed-parameter DFS (the two non-commutative ones).

- **Distance-regular graphs with a primitive group** (`search/drg/drg_note.md`). All 4244 primitive groups (not Sym/Alt) on the 156 degrees below 4096 that carry a primitive BCN array of diameter >= 3. 86 of the 117 arrays that survive the suborbit filter are realised, all by known graphs. For the other 31 no primitive group works: 8 are known not to exist, one depends on Koval's preprint, three exist with imprimitive groups (Shi-Krotov-Sole graph, Brouwer graph Br(3), a coset graph on 1024 vertices; the last two built and checked here), one is reported nonexistent but unconfirmed, and for 18 I found no status, so any graph with one of those arrays has an imprimitive automorphism group.
- **The 18 open arrays on abelian groups** (`search/cayley/cayley_note.md`). 13 of them fail the Krein integrality test for translation schemes. Four more ({52,42,16;1,6,28}, {49,40,22;1,5,28}, {54,42,11;1,3,44}, {63,48,24;1,6,28}) are ruled out on every abelian group by a quotient of order 125, checked with CP-SAT and SCIP. For {56,42,20;1,6,28} only Z_9 x Z_3^4 is left. No new graph.

## Layout

| Path | Content |
|---|---|
| `fusion/` | `fusion_note.tex` (draft); `verify_fusions.py` (exact sympy check for all integer (s,t) in [2,30]^2, hexagon cubic, Tits data, Krein and SRG parameters); `certify_fusions.py` (symbolic (s,t) proof); `near_polygon_fusions.py` (regular near polygons); `near_octagon_decompose.py`, `near_octagon_check.py`, `near_octagon_families.py` (regular near octagons); notes and outputs |
| `qpoly/` | Q-polynomial ordering check for point schemes of generalized 2d-gons, d = 2,3,4 |
| `thin/` | Fusions of thin octagon schemes of order (1,t) and (s,1); `verify_t10.py` |
| `tits1755/` | The 1755-point scheme (`scheme_from_array.py`, `compute_tits1755.py`, `tightness_check.py`), coset schemes of the maximal subgroups (`coset_schemes.py`, `cosets/`, `out_cosets/`), GAP scripts in `gap/` |
| `search/` | SRG/DRG search over unions of orbitals for T-actions (`srg_drg_search.py`, `dfs_numba.py`), fusions of the non-commutative schemes (`fusion_closure.py`), fixed-parameter search (`dfs_fixed.py`), SRG search over primitive groups (`primitive/`), DRG search over primitive groups (`drg/`) |
| `drgfusion/` | Fusions of every array in Brouwer's DRG tables, SRG cross-check, feasibility sweep |
| `triple/` | Triple intersection number LPs for GO(2,4) |
| `code/`, `data/` | ATLAS orbital table, Bannai-Muzychuk cross-check, and the CSVs they and the other scripts write |
| `docs/` | Report on the 1755-point scheme, literature notes |

Most folders have a note next to the scripts saying how to run them and what the output means.

## Reproducing

    pip install -r requirements.txt
    python fusion/verify_fusions.py
    python qpoly/qpoly_check.py
    python tits1755/scheme_from_array.py --atlas
    python code/atlas_orbital_table.py
    python code/bannai_muzychuk_crosscheck.py     # about 7 minutes

The GAP side needs GAP 4.16 with ctbllib and atlasrep (`sh tits1755/gap/run_gap.sh`). `drgfusion/fetch_and_crosscheck.py`, `code/atlas_orbital_table.py` and `triple/go24_triple_lp.py` download data, so they need network access.

What has been checked on the current commit: `scheme_from_array.py --atlas` exits 0, the 125 p^k_ij from the intersection array equal the counts from the ATLAS permutation action, and the distance in R_1 equals the orbital label at all 1755 points. Both intersection tables in `data/` and `tits1755/out/` agree with it cell by cell. `bannai_muzychuk_crosscheck.py` compiles, but I haven't re-run it since the last change; the CSV is from the previous run.

## Known limitations

- `fusion_note.tex`: the first real compile (pdfTeX, TeX Live 2025 under WSL) stopped on a missing `booktabs.sty`. The preamble now loads booktabs, array, longtable and hyperref only if present. Not recompiled since. Cites and refs resolve and braces balance. The hand proof in Section 3 still uses the sign analysis of eigenvalue differences; the exact proof is cited in Section 7.
- Statements about what is new come from keyword searches (`docs/literature.md`), not a proper review.
- In `coset_schemes.py` the eigenmatrices and Krein parameters are floating point. Fusion tests, intersection numbers and centre dimensions are exact.
- The permutations in `tits1755/cosets/` come from GAP's `FactorCosetAction`, so point numbering depends on GAP and the class labels in `out_cosets/` refer to that numbering.
- Character table checks (GAP 4.16.1, ctbllib 1.3.11, `tits1755/gap/tits_chars.g`): the library table of 2F4(2)' has 22 classes and its degrees square-sum to |T| = 17971200. The permutation character on 1755 points is 1 + 78 + 351 + 650 + 675, multiplicity-free with <pi,pi> = 5, and each eigenspace character computed exactly from Q and the fixed-point counts equals the matching library irreducible. `tits1755/verify_characters.py` repeats this in Python against `gap/tits_chartable.csv`, and a Dixon-Schneider table computed from the generators agrees with the library one (`gap/tits_ctG_output.txt`).
- The two intersection tables in `tits1755/out/` and `data/` have different column names (`p_k_ij` and `p^k_ij`) and the same values, since two scripts write them.

## To do

- Compile `fusion/fusion_note.tex` and rewrite Section 3 around the exact proof.
- Settle the three residual surfaces in `fusion/near_octagon_note.md`.
- Check the literature for the 18 arrays in `search/drg/drg_note.md` with no known status, and find the paper on {20,16,5;1,1,16}.
