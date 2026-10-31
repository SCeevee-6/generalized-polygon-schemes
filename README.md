# generalized-polygon-schemes

Association-scheme computations for point graphs of generalized polygons: fusions of the point-distance scheme, Q-polynomial orderings, and the orbital algebra of the Tits group on 1755 points.

Status: project working notes, not peer reviewed. `fusion/fusion_note.tex` is a draft.

## Results

- **Fusions.** For a generalized 2d-gon of order (s,t), d in {3,4}, with integers s,t >= 2, the point-distance scheme has a proper fusion if and only if s = t. The fusion is then unique: {1,...,d-1} | {d}.
- **Hexagons, s = t = q.** The fusion is strongly regular with parameters ((q^6-1)/(q-1), q(q+1)(q^2+1), q^3+q^2+q-1, (q+1)(q^2+1)), the parameters of the collinearity graph of Q(6,q).
- **Octagons.** The case s = t is excluded by the Feit-Higman condition (2st must be a square), so the Ree-Tits octagon GO(2,4) has no proper fusion.
- **Tits group on 1755 points, (s,t) = (2,4).** Intersection array {10,8,8,8;1,1,1,5}. Eigenvalues of A_1 are 10, 5, 1, -3, -5 with multiplicities 1, 351, 650, 675, 78. The only vanishing non-principal Krein parameter is q^4_44 = 0. Of the 15 fusions of {1,2,3,4}, only the two trivial ones close.

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
| `thin/` | Thin-octagon GO(1,t) fusion check (t = 2, 8, 18, 32, 50): only trivial, bipartite-parity and halved-GQ(t,t) fusions close; no primitive fusion; 0/24 Q-polynomial orderings. `thin_fusions.py`, `thin_fusions.json`, `verify_thin_explicit.py`, `thin_note.md` | `python thin/thin_fusions.py` |
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
- `fusion/fusion_note.tex` was only checked statically; no TeX engine was run on it. Bibliography: Crnković–Rukavina–Švob (pages corrected to 641–652), Cameron–Goethals–Seidel, Feit–Higman, Van Maldeghem and Brouwer–Cohen–Neumaier (Springer 1989, doi:10.1007/978-3-642-74341-2) match Crossref. For Tits (1983), Crossref gives the DOI 10.2307/2374268 and start page 539; the range 539–594 matches several independent reference lists. Bannai–Ito (1984) is not in Crossref and is still unchecked.
- Literature novelty claims come from a search, not an exhaustive review.
- Verified (`tits1755/verify_group_order.py`, output in `tits1755/verify_group_order_output.txt`): |T| = 17971200 by Schreier–Sims on the ATLAS generators, matching the ATLAS page order; the 5 eigenspaces of A_1 (dimensions 1, 351, 650, 675, 78) are irreducible, because the rank is 5 and the permutation character is therefore multiplicity-free with 5 constituents. Checked against the character table by `tits1755/verify_characters.py` (output in `tits1755/verify_characters_output.txt` and `tits1755/character_values_by_order.csv`): the eigenspace characters chi_i(g) = (1/n) sum_j Q_ji N_j(g) were computed for 4000 random elements, covering every element order 2, 3, 4, 5, 6, 8, 10, 12, 13, 16. All values are rational integers. The mean of chi^2 is 1.0, 0.92, 0.90, 0.83, 0.95, consistent with 1 within sampling error (not quantified). The characters of degree 351, 650 and 78 vanish on elements of order 13, and those of degree 650 and 675 vanish on elements of order 5, as 13- and 5-defect-zero characters must. The degree list 1, 26², 27², 78, 300, 325, 351³, 624², 650, 675, 702², 1300², 1728, 2048² has sum of squares 17971200 = |T| and contains all five degrees. That list was not taken from a downloaded copy of the ATLAS table.
- The two intersection tables in `tits1755/out/` and `data/` have different column names (`p_k_ij` and `p^k_ij`) and the same values. They are kept because two different scripts write them.

## Open directions

- Q-polynomial classification of hexagon and octagon point schemes. `qpoly/` is the start of this.
- Whether the fusion result extends to other coherent configurations of T on Omega.
