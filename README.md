# generalized-polygon-schemes

Association-scheme computations for point graphs of generalized polygons: fusions of the point-distance scheme and Q-polynomial orderings.

Status: project working notes, not peer reviewed. `fusion_note.tex` is a draft.

## Results

- **Fusions.** For a generalized 2d-gon of order (s,t), d in {3,4}, with integers s,t >= 2, the point-distance scheme has a proper fusion if and only if s = t. The fusion is then unique: {1,...,d-1} | {d}.
- **Hexagons, s = t = q.** The fusion is strongly regular with parameters ((q^6-1)/(q-1), q(q+1)(q^2+1), q^3+q^2+q-1, (q+1)(q^2+1)), the parameters of the collinearity graph of Q(6,q).
- **Octagons.** The case s = t is excluded by the Feit-Higman condition (2st must be a square), so the Ree-Tits octagon GO(2,4) has no proper fusion.
- **Tits group on 1755 points, (s,t) = (2,4).** Eigenvalues of A_1 are 10, 5, 1, -3, -5. Multiplicities are 1, 351, 650, 675, 78. The only vanishing non-principal Krein parameter is q^4_44 = 0.

## Files

| Path | Content | How to reproduce |
|---|---|---|
| `verify_fusions.py` | Exact sympy checks: eigenmatrix recursion, direct fusion test on all integer (s,t) in [2,30]^2, hexagon cubic, Tits-group data, Krein parameters, SRG parameters | `python verify_fusions.py`; expected stdout in `verify_fusions_output.txt` |
| `code/atlas_orbital_table.py` | Orbital algebra of T on 1755 points from the ATLAS v3 permutation generators (Schreier generators, union-find suborbits, transporter counts) | Needs internet access to brauer.maths.qmul.ac.uk. Writes `intersection_numbers_independent_verification.csv` to the working directory |
| `code/bannai_muzychuk_crosscheck.py` | Groebner-basis column criterion over the 13 proper partitions, listing real-positive (s,m) solutions | About 7 minutes. Writes `bannai_muzychuk_crosscheck.csv` to the working directory |
| `data/intersection_numbers_independent_verification.csv` | 125 intersection numbers p^k_ij of the 1755-point scheme | Output of `code/atlas_orbital_table.py` |
| `data/bannai_muzychuk_crosscheck.csv` | Real-positive solutions per partition | Output of `code/bannai_muzychuk_crosscheck.py` |
| `data/krein_parameters_independent.csv` | Krein parameters of the 1755-point scheme | From an earlier session. `verify_fusions.py` re-derives only the vanishing parameter |
| `fusion_note.tex` | Draft note: theorems, proofs, Tits-group tables, Appendix with the octagon sign table | Compile with pdflatex or Overleaf (standard packages) |
| `docs/novelty_audit.md` | Literature audit from an earlier session | Predates the unified theorem; its octagon-only framing is superseded |
| `qpoly/` | Q-polynomial ordering check from an earlier session (`qpoly_check.py`, output, `qpoly_note.md`) | Not re-run in this commit |

## Repository layout from earlier commits

Earlier commits on `main` (same account, 22:44) added folders that overlap with the files above:

- `fusion/`: `verify_fusions.py`, `verify_fusions_output.txt` and `fusion_note.tex`. These are byte-identical to the root copies.
- `tits1755/`: `compute_tits1755.py` (earlier ATLAS orbital script) and `out/` (intersection numbers, Krein parameters, fusion test, summary). `code/atlas_orbital_table.py` is the version reproduced in this commit.

Both folders were left as they were. The root copies are the ones this README describes.

## Reproducing

    pip install -r requirements.txt
    python verify_fusions.py
    python code/atlas_orbital_table.py
    python code/bannai_muzychuk_crosscheck.py

Checks run on the contents of this commit:

- `verify_fusions.py` stdout equals `verify_fusions_output.txt`.
- `code/atlas_orbital_table.py` run against the live ATLAS generators reproduces `data/intersection_numbers_independent_verification.csv` exactly (125 rows).
- `code/bannai_muzychuk_crosscheck.py` compiles. It was not re-run, because the run takes about 7 minutes. Its CSV comes from the producing run.

## Known limitations

- `code/bannai_muzychuk_crosscheck.py` writes the last column, `feasible_thick_octagon`, as the constant `none`. The column is not computed. Checked by hand against the solution list: each real-positive entry has s = 1, t = 1, a non-integer s or t, or lies on the diagonal s = t (the curve written `M**2 - 2*S**2`), which the Feit-Higman condition excludes. The computed feasibility test is the grid test in `verify_fusions.py`.
- `fusion_note.tex` was checked statically only; no TeX engine ran in the environment that produced it. Bibliography entries other than Crnkovic-Rukavina-Svob (2020) and Cameron-Goethals-Seidel (1978) were written from memory and need verification.
- Literature novelty claims come from a search, not an exhaustive review.
- Not verified: |T| by Schreier-Sims, and the multiplicities against the ATLAS character table.

## Open directions

- Q-polynomial classification of hexagon and octagon point schemes. The `qpoly/` folder is the start of this.
- Whether the fusion result extends to other coherent configurations of T on Omega.
