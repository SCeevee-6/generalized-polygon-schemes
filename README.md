# generalized-polygon-schemes

Association-scheme computations for point graphs of generalized polygons: fusions of the point-distance scheme, Q-polynomial orderings, and the orbital algebra of the Tits group on 1755 points.

Status: project working notes, not peer reviewed. `fusion/fusion_note.tex` is a draft.

## Results

- **Fusions.** For a generalized 2d-gon of order (s,t), d in {3,4}, with integers s,t >= 2, the point-distance scheme has a proper fusion if and only if s = t. The fusion is then unique: {1,...,d-1} | {d}. This now has an exact proof from the intersection array (Groebner bases over QQ, exact real-root isolation, plus an integer scan): `fusion/certify_fusions.py`, `fusion/certify_note.md`. The hexagon case s = t is known; for the octagon statement I found no source, but my search was not exhaustive.
- **Hexagons, s = t = q.** The fusion is strongly regular with parameters ((q^6-1)/(q-1), q(q+1)(q^2+1), q^3+q^2+q-1, (q+1)(q^2+1)), the parameters of the collinearity graph of Q(6,q).
- **Octagons.** The case s = t is excluded by the Feit-Higman condition (2st must be a square), so the Ree-Tits octagon GO(2,4) has no proper fusion.
- **Tits group on 1755 points, (s,t) = (2,4).** Intersection array {10,8,8,8;1,1,1,5}. Eigenvalues of A_1 are 10, 5, 1, -3, -5 with multiplicities 1, 351, 650, 675, 78. The only vanishing non-principal Krein parameter is q^4_44 = 0. Of the 15 fusions of {1,2,3,4}, only the two trivial ones close. The scheme is not Q-polynomial in any ordering (known, De Bruyn–Vanhove 2015) and the graph is not tight (Jurišić–Koolen–Terwilliger bound 0 >= -20 is strict; `tits1755/tightness_check.py`).
- **Other actions of T.** On cosets of the 8 maximal subgroups, the orbital configurations are commutative for degrees 1600 (twice), 1755, 2304 and 2925. All of these are symmetric except 2304 (L2(25)). They are non-commutative for A6.2^2 (degree 12480, rank 24, twice) and 5^2:4A4 (degree 14976, rank 25). See `tits1755/coherent_configurations_note.md`.
- **Fusions of the other commutative actions** (`tits1755/coset_fusions_note.md`, computed from the explicit coset actions). Degree 1600: the one proper fusion is the known SRG(1600,351,94,72) (Crnković–Rukavina–Švob 2020). Degree 2304: three fusions, the T.2 orbital scheme (rank 6, non-symmetric), the symmetrisation (rank 6), and the symmetrisation of the T.2 scheme (rank 5); none is P- or Q-polynomial. Degree 2925: one fusion, the distance scheme of the dual octagon GO(4,2), which is the T.2 orbital scheme; P- but not Q-polynomial.
- **Non-commutative actions.** The orbital algebras are 11 C + M_2(C) + M_3(C) (degree 12480, both classes) and 9 C + 4 M_2(C) (degree 14976). Centre dimension 13 in both cases, computed exactly from the structure constants, which agrees with the character multiplicities.
- **SRGs and DRGs with a transitive T-action, beyond Crnković–Rukavina–Švob** (`search/search_note.md`). Exhaustive search over all unions of orbitals for every T-action of degree 7020 to 17550 (TomLib classes), including the primitive actions on 12480 and 14976 points. No SRG and no DRG of diameter >= 3 on 12480 or 14976 points, so with CRS 2020 the primitive case is complete. On 7020, 8775, 11700 and two of the three 17550-point actions only complete multipartite SRGs occur. The third 17550-point action (rank 54) is still running. New as a search range; no new graph.
- **Forced fusions in Brouwer's DRG tables** (`drgfusion/drgfusion_note.md`). All closed fusions of all tabled arrays of diameter <= 7, with every fused SRG on at most 1300 vertices looked up in Brouwer's SRG table. Besides the three arrays already ruled out this way on Brouwer's page, six more come out: {14,10,3;1,5,12} (fused SRG violates the absolute bound), {19,12,5;1,4,15} (Degraer 2007), and, conditional on Kovaľ's Oct 2026 preprints, {39,30,4;1,5,36}, {112,77,16;1,16,88}, {21,16,10,1;1,2,16,21}, {25,24,2,1;1,2,24,25}. The method is known; I did not find these six arrays ruled out in the sources I checked.
- **Fusions of the non-commutative T schemes** (`search/fusions/fusion_note.md`). The 12480-point configurations (both classes) have no proper fusion of rank >= 3 at all. Every fusion of rank >= 3 of the 14976-point configuration refines the rank-16 2F4(2) orbital configuration and is non-commutative. New as far as I found, and small.
- **Regular near hexagons and octagons** (`fusion/near_polygon_note.md`). Diameter 3 is BCN 4.2.17 (known). For near octagons a box scan (s <= 8, t <= 24) finds closing thick fusions only for generalized octagons with s = t (excluded by Feit–Higman), Hamming graphs and the B_4/C_4 dual polar graphs, all known. No general theorem.
- **GO(2,4) triple intersection numbers** (`triple/triple_note.md`). With q^4_44 = 0 the LP ranges of the triple intersection numbers only reproduce octagon geometry for one triple type. Nothing new for uniqueness.

## Layout

| Path | Content | How to run |
|---|---|---|
| `fusion/` | Draft note `fusion_note.tex`; exact sympy verifier `verify_fusions.py` (all integer (s,t) in [2,30]^2, hexagon cubic, Tits data, Krein parameters, SRG parameters); expected stdout `verify_fusions_output.txt` | `python fusion/verify_fusions.py` (sympy) |
| `qpoly/` | Q-polynomial ordering check for point schemes of generalized 2d-gons, d = 2,3,4 (`qpoly_check.py`, output, `qpoly_note.md`). Not re-run in the latest commit. | `python qpoly/qpoly_check.py` (sympy) |
| `tits1755/scheme_from_array.py` | Derivation of the 1755-point algebra from the intersection array alone (stdlib, exact Fractions): p^k_ij, P, Q, multiplicities, Krein parameters, Q-polynomial test over 24 orderings, closure-violation counts for all 15 fusions. `--atlas` checks the result against the ATLAS generators. | `python tits1755/scheme_from_array.py --atlas` (about 10 s) |
| `tits1755/out_array/` | Outputs of `scheme_from_array.py`: `orbital_valencies.csv`, `intersection_numbers.csv`, `P_matrix.csv`, `Q_matrix.csv`, `krein_parameters.csv`, `fusion_results.csv`, `run_output.txt` | — |
| `tits1755/figures/` | `tits_1755_scheme.png`: intersection diagram of the 5 orbitals and closure violations for the 15 fusions | — |
| `tits1755/compute_tits1755.py`, `tits1755/octagon_fusion_symbolic.py`, `tits1755/out/` | Earlier computation from the raw ATLAS generators, with its outputs (`intersection_numbers.csv`, `krein_parameters.csv`, `fusion_test.csv`, `summary.json`), and a symbolic octagon fusion test | `python tits1755/compute_tits1755.py` (needs the ATLAS generators) |
| `code/` | `atlas_orbital_table.py`: orbital table from the ATLAS v3 permutation generators. `bannai_muzychuk_crosscheck.py`: Groebner-basis column criterion over the 13 proper partitions (about 7 minutes). | see Reproducing |
| `data/` | `intersection_numbers_independent_verification.csv`, `krein_parameters_independent.csv`, `bannai_muzychuk_crosscheck.csv`, `symbolic_fusion_test.csv` | outputs of the scripts above |
| `thin/` | Fusions of thin octagon schemes, orders (1,t) for t = 2, 8, 18, 32, 50 and (s,1) for s = 2, 8. For (1,t) the only proper fusions are the bipartition {1,3}\|{2,4} and the halved fusion {2}\|{1,3}\|{4}, whose class 2 is two copies of the GQ(t,t) collinearity graph. Order (2,1) has the primitive fusion {3}\|{1,2,4}, which is SRG(45,16,8,4), the parameters of T(10). Order (8,1) has none. None is Q-polynomial. `thin_fusions.py` writes `out_thin/thin_fusions.json`; `verify_thin_explicit.py` checks W(2), W(8); `thin_note.md` | `python thin/thin_fusions.py` (under 1 s) |
| `tits1755/gap/` | GAP checks: `tits_chars.g` (library table, permutation character, exact eigenspace characters), `tits_ctG.g` (Dixon–Schneider table versus the library), `tits_maxes.g` (coherent configurations on cosets of the 8 maximal subgroups), with their outputs and `maxes_coherent_configurations.csv`. Notes in `tits1755/coherent_configurations_note.md` | `sh tits1755/gap/run_gap.sh` (GAP 4 with ctbllib and atlasrep; fetches the ATLAS files) |
| `tits1755/coset_schemes.py`, `tits1755/cosets/`, `tits1755/out_cosets/` | Orbital schemes of T on cosets of maxes 1, 2, 4, 5, 6, 7, 8 from the two generator images (written by `tits1755/gap/coset_actions.g`): exact stabiliser orbits from Schreier generators, integer p^k_ij, centre dimension, eigenmatrix, P/Q-polynomial tests, exhaustive fusion test. Note: `tits1755/coset_fusions_note.md` | `python tits1755/coset_schemes.py tits1755/cosets/max5.txt tits1755/out_cosets/max5` (numpy, scipy, sympy; about 50 s for degree 14976) |
| `fusion/certify_fusions.py` | Exact fusion test for generalized hexagons and octagons with symbolic (s,t); output `certify_fusions_output.txt`, note `certify_note.md` | `python fusion/certify_fusions.py` (run from `fusion/`, about 15 s) |
| `search/` | `srg_drg_search.py` (exhaustive SRG/DRG search over unions of orbitals: Gray code, or a CP-SAT model with no-good cuts), `build_scheme.py` (orbital scheme from two generator images), `validate_cuts.py`; schemes in `search/schemes/`, results in `search/out/`; note `search_note.md`. Coset actions of the TomLib classes in `tits1755/cosets_tom/` (GAP: `tits1755/gap/tom_coset_actions.g`) | `python search/srg_drg_search.py search/schemes/tom422_intersection_numbers.csv out.json gray` (numpy; `cuts` needs ortools) |
| `drgfusion/` | `drg_table_fusions.py` (closed fusions of every array in Brouwer's DRG tables), `fetch_and_crosscheck.py` (downloads the tables, looks up fused SRGs in Brouwer's SRG table); outputs `table_fusions.csv`, `srg_crosscheck.csv`, `nonexistent_fusions.csv`; note `drgfusion_note.md` | `cd drgfusion; python fetch_and_crosscheck.py` (about 3 minutes, needs network) |
| `search/fusion_closure.py`, `search/fusions/` | Fusions generated by one symmetric relation (coherent closure inside the orbital algebra), sharded runs for 12480 and 14976 points; note `fusion_note.md` | `python search/fusion_closure.py tits1755/out_cosets/max8_intersection_numbers.csv out.json 0 4` (four shards, about 30 min each), then `merge` |
| `fusion/near_polygon_fusions.py` | Symbolic p^k_ij of regular near hexagons/octagons (s, t_2, ..., t), Groebner bases of all closure ideals, integer box scan; output `near_polygon_fusions_output.txt`, note `near_polygon_note.md` | `python fusion/near_polygon_fusions.py 8 24` (about 40 min) |
| `triple/` | `go24_triple_lp.py`: LP ranges of triple intersection numbers of {10,8,8,8;1,1,1,5} with q^4_44 = 0, checked on the ATLAS octagon; output and `triple_note.md` | `python triple/go24_triple_lp.py` (needs network for the ATLAS generators) |
| `docs/handoff_next_chat.md` | Prompt for the next working session: angles aimed at a new result | — |
| `docs/` | `tits_1755_orbital_algebra_report.md`: report on the 1755-point scheme. `novelty_audit.md`: literature audit from an earlier session; its octagon-only framing is superseded by the unified theorem, and its exceptional-parameter table by `fusion/certify_note.md`. `literature_check_2026-10-31.md`: Crossref/OpenAlex search (raw hits in `data/lit_search_2026-10-31.csv`). | — |

The earlier copies of `fusion_note.tex`, `verify_fusions.py` and `verify_fusions_output.txt` at the repository root were removed. They were byte-identical to the copies in `fusion/`.

## Reproducing

    pip install -r requirements.txt
    python fusion/verify_fusions.py
    python qpoly/qpoly_check.py
    python tits1755/scheme_from_array.py --atlas
    python code/atlas_orbital_table.py
    python code/bannai_muzychuk_crosscheck.py

Checks run on the contents of this commit:

- `tits1755/scheme_from_array.py --atlas` exits 0. The 125 p^k_ij derived from the intersection array equal the counts taken from the ATLAS permutation action. Graph distance in R_1 equals the orbital label at all 1755 points.
- The intersection tables in `data/intersection_numbers_independent_verification.csv` and `tits1755/out/intersection_numbers.csv` both agree with the table derived by `scheme_from_array.py`, cell by cell (125 entries each).
- `code/bannai_muzychuk_crosscheck.py` compiles. It was not re-run, because the run takes about 7 minutes. Its CSV comes from the producing run.

## Known limitations

- `code/bannai_muzychuk_crosscheck.py` now computes the last column, `feasible_thick_octagon`. Isolated solutions are tested for integer s, t >= 2 with 2st a perfect square. For curves, the curve is solved for m and every integer s in [2, 5000] is tested. All 13 partitions come out `none`, the same answer as the earlier hand check. Sanity check: the point (s, m) = (8, 8) is correctly flagged as (s, t) = (8, 4).
- `fusion/fusion_note.tex`: the first real compile (pdfTeX, TeX Live 2025 under WSL) stopped because `booktabs.sty` is not installed there. The preamble now loads booktabs, array, longtable and hyperref only if they exist and falls back to plain LaTeX otherwise. It has not been recompiled since that change. Static checks pass: every cite and ref resolves, and braces and environments balance.
- Literature novelty claims come from a search, not an exhaustive review (`docs/literature_check_2026-10-31.md`).
- `fusion/fusion_note.tex` keeps its hand proof (sign analysis of eigenvalue differences). Section 7 now also cites the exact Groebner check `fusion/certify_fusions.py`, and the hexagon cubic is factored explicitly as (s^2+st+t^2)(st-(t+1)^2).
- In `coset_schemes.py` the eigenmatrices and Krein parameters are floating point. The fusion test, intersection numbers and centre dimensions are exact.
- The coset-action permutations in `tits1755/cosets/` come from GAP's `FactorCosetAction`; the point numbering depends on GAP, so class labels in `out_cosets/` refer to that numbering.
- Character table, checked with GAP 4.16.1 and ctbllib 1.3.11 (`tits1755/gap/tits_chars.g`, output in `tits1755/gap/tits_chars_output.txt`). The library table of 2F4(2)' has 22 classes, and its sorted degrees equal the list 1, 26², 27², 78, 300, 325, 351³, 624², 650, 675, 702², 1300², 1728, 2048², whose squares sum to |T| = 17971200. The ATLAS class words give representatives whose orders and centralizer orders match the library, and no two of them are conjugate. The permutation character on 1755 points, computed both from fixed points and as 1_H^T from the library max 2.2^8.f20 (order 10240), is 1 + 78 + 351 + 650 + 675 (library X.1, X.6, X.9, X.14, X.15). It is multiplicity-free, with <pi,pi> = 5. Each eigenspace character chi_e(g) = (1/n) sum_j Q_je N_j(g), evaluated exactly on all 22 classes, equals that library irreducible. GAP does this in `tits_chars.g`, and `tits1755/verify_characters.py` repeats it in Python against `gap/tits_chartable.csv`. A table computed from scratch from the generators by Dixon–Schneider is equivalent to the library table (`gap/tits_ctG_output.txt`). An earlier version of this README said the 4000-sample check compared against the character table. That was wrong, and the sampling script has been replaced.
- The two intersection tables in `tits1755/out/` and `data/` have different column names (`p_k_ij` and `p^k_ij`) and the same values. They are kept because two different scripts write them.

## Open directions

- Q-polynomial question: settled. Thick octagons are covered by De Bruyn–Vanhove (2015). `qpoly/` reproduces it.
- Fusions of the commutative actions (degrees 1600, 2304, 2925): done. Fusions of the non-commutative schemes (23 and 24 non-trivial classes) were not searched.
- SRG/DRG search for T: done up to 17550 vertices except one action (TomLib class 417, rank 54), which is running.
- Forced fusions in the DRG tables: done for diameter <= 7 and fused SRGs on <= 1300 vertices. The six arrays in `drgfusion/drgfusion_note.md` should go to Brouwer, after a check of BCN's printed tables and Vidali's thesis.
- Near polygons, non-commutative fusions and triple intersection numbers: done to the extent described in the notes; none of them gave a publishable new result.
- Compile `fusion/fusion_note.tex` (needs a TeX installation outside the sandbox) and rewrite its Section 3 around the exact proof.
