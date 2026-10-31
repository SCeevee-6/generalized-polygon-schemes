# generalized-polygon-schemes

Association-scheme computations for point graphs of generalized polygons: fusions of the point-distance scheme, Q-polynomial orderings, and the orbital algebra of the Tits group on 1755 points.

Status: project working notes, not peer reviewed. `fusion/fusion_note.tex` is a draft.

## Results

- **Fusions.** For a generalized 2d-gon of order (s,t), d in {3,4}, with integers s,t >= 2, the point-distance scheme has a proper fusion if and only if s = t. The fusion is then unique: {1,...,d-1} | {d}.
- **Hexagons, s = t = q.** The fusion is strongly regular with parameters ((q^6-1)/(q-1), q(q+1)(q^2+1), q^3+q^2+q-1, (q+1)(q^2+1)), the parameters of the collinearity graph of Q(6,q).
- **Octagons.** The case s = t is excluded by the Feit-Higman condition (2st must be a square), so the Ree-Tits octagon GO(2,4) has no proper fusion.
- **Tits group on 1755 points, (s,t) = (2,4).** Intersection array {10,8,8,8;1,1,1,5}. Eigenvalues of A_1 are 10, 5, 1, -3, -5 with multiplicities 1, 351, 650, 675, 78. The only vanishing non-principal Krein parameter is q^4_44 = 0. Of the 15 fusions of {1,2,3,4}, only the two trivial ones close. The scheme is not Q-polynomial in any ordering (known, De Bruyn–Vanhove 2015) and the graph is not tight (Jurišić–Koolen–Terwilliger bound 0 >= -20 is strict; `tits1755/tightness_check.py`).
- **Other actions of T.** On cosets of the 8 maximal subgroups, the orbital configurations are commutative for degrees 1600 (twice), 1755, 2304 and 2925. All of these are symmetric except 2304 (L2(25)). They are non-commutative for A6.2^2 (degree 12480, rank 24, twice) and 5^2:4A4 (degree 14976, rank 25). See `tits1755/coherent_configurations_note.md`.

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
| `docs/` | `tits_1755_orbital_algebra_report.md`: report on the 1755-point scheme. `novelty_audit.md`: literature audit from an earlier session; its octagon-only framing is superseded by the unified theorem. | — |

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
- `fusion/fusion_note.tex` has not been compiled yet. tectonic and TinyTeX pdflatex both fail inside the sandbox on this Windows host, so the file was checked only statically: every cite and ref resolves, and braces and environments balance. Bibliography: Crnković–Rukavina–Švob (pages corrected to 641–652), Cameron–Goethals–Seidel, Feit–Higman, Van Maldeghem, Brouwer–Cohen–Neumaier (Springer 1989, doi:10.1007/978-3-642-74341-2), De Bruyn–Vanhove (Combinatorica 35 (2015) 181–208) and Jurišić–Koolen–Terwilliger (J. Algebraic Combin. 12 (2000) 163–197) match Crossref. For Tits (1983), Crossref gives the DOI 10.2307/2374268 and start page 539; the range 539–594 matches several independent reference lists. Bannai–Ito is not in Crossref. It was checked against the Stanford SearchWorks and Trove (NLA) catalogue records: Benjamin/Cummings, Menlo Park, 1984, Mathematics Lecture Note Series 58, xxiv+425 pp., ISBN 0-8053-0490-8 (hardcover) and 0-8053-0491-6 (paperback).
- Literature novelty claims come from a search, not an exhaustive review.
- Character table, checked with GAP 4.16.1 and ctbllib 1.3.11 (`tits1755/gap/tits_chars.g`, output in `tits1755/gap/tits_chars_output.txt`). The library table of 2F4(2)' has 22 classes, and its sorted degrees equal the list 1, 26², 27², 78, 300, 325, 351³, 624², 650, 675, 702², 1300², 1728, 2048², whose squares sum to |T| = 17971200. The ATLAS class words give representatives whose orders and centralizer orders match the library, and no two of them are conjugate. The permutation character on 1755 points, computed both from fixed points and as 1_H^T from the library max 2.2^8.f20 (order 10240), is 1 + 78 + 351 + 650 + 675 (library X.1, X.6, X.9, X.14, X.15). It is multiplicity-free, with <pi,pi> = 5. Each eigenspace character chi_e(g) = (1/n) sum_j Q_je N_j(g), evaluated exactly on all 22 classes, equals that library irreducible. GAP does this in `tits_chars.g`, and `tits1755/verify_characters.py` repeats it in Python against `gap/tits_chartable.csv`. A table computed from scratch from the generators by Dixon–Schneider is equivalent to the library table (`gap/tits_ctG_output.txt`). An earlier version of this README said the 4000-sample check compared against the character table. That was wrong, and the sampling script has been replaced.
- The two intersection tables in `tits1755/out/` and `data/` have different column names (`p_k_ij` and `p^k_ij`) and the same values. They are kept because two different scripts write them.

## Open directions

- Q-polynomial question: settled. Thick octagons are covered by De Bruyn–Vanhove (2015). `qpoly/` reproduces it.
- Fusions of the other commutative actions (degrees 1600, 2304, 2925) and of the three non-commutative ones were not computed.
