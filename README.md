# generalized-polygon-schemes

Association-scheme computations for point graphs of generalized polygons: fusion closure,
Q-polynomial orderings, and the orbital algebra of the Tits group on 1755 points.

## Layout

| Path | Contents |
|---|---|
| `tits1755/compute_tits1755.py` | From the ATLAS v3 generators of T = 2F4(2)' (TF42G1-p1755B0): suborbits, 125 intersection numbers, P, Q, Krein parameters, Q-polynomial test, 15-partition fusion test |
| `tits1755/octagon_fusion_symbolic.py` | Symbolic Bose-Mesner algebra of a generalized octagon GO(s,t); fusion closure differences per partition; cross-check against the numeric table |
| `tits1755/out/` | Outputs of the two scripts (`summary.json`, CSV tables, `octagon_fusion_symbolic.txt`) |
| `fusion/` | Earlier fusion verification for generalized hexagons and octagons: `verify_fusions.py`, `verify_fusions_output.txt`, `fusion_note.tex` |
| `qpoly/` | Q-polynomial check for generalized polygon schemes (`qpoly_check.py`, output, note) |
| `docs/novelty_audit.md` | Literature audit of the 1755-point project and the claimed new result |
| `docs/tits_1755_orbital_algebra_report.md` | Orbital-algebra report for the 1755-point scheme (version 2) |

## Reproduce

```
pip install -r requirements.txt        # sympy
python tits1755/compute_tits1755.py    # ~10 s; fetches ATLAS generators into tits1755/cache/
python tits1755/octagon_fusion_symbolic.py          # ~5 s; add --solve for sympy.solve (slow)
```

`cache/` is git-ignored because it holds third-party ATLAS files.

## Results reproduced by the scripts

- Suborbits of the point stabiliser: 1, 10, 80, 640, 1024.
- Eigenvalues 10, 5, 1, -3, -5 with multiplicities 1, 351, 650, 675, 78.
- BFS distance from 0 on the R_1 graph equals the suborbit label at all 1755 points.
- All Krein parameters are non-negative. The only vanishing Krein parameter with all indices >= 1 is q^4_44 = 0 (eigen-indices).
- No ordering of the eigenspaces passes the triangle-zero Q-polynomial test.
- Of the 15 partitions of {1,2,3,4}, only the trivial and identity partitions have zero closure violations.
- The symbolic octagon matrices reproduce all 125 numeric intersection numbers at (s,t) = (2,4).

## Verification status and caveats

- The graph is published in Crnkovic, Rukavina and Svob, J. Algebraic Combin. 51 (2020), DOI 10.1007/s10801-019-00888-5 (arXiv:1809.10197). Its Table 13 lists the same intersection array and automorphism group. The Krein, fusion and Q-polynomial results are not in that paper.
- Not checked against a primary source in this work: the Cameron-Goethals-Seidel (1978) extremal-polygon result on vanishing Krein parameters, and the De Bruyn (2022) statement of its octagon case. The publisher page was blocked by the network allowlist during the audit.
- The tightness claim in the audit uses an identity that was not checked against a primary source. It is not computed by the scripts.
- The character-degree attribution of 1, 351, 650, 675, 78 to the ATLAS character table was not re-checked.
- The full group order |T| = 17,971,200 was not verified by Schreier-Sims. The code checks transitivity and the suborbit structure only.
- The octagon fusion theorem in `docs/novelty_audit.md` concerns thick generalized octagons. For formal s = t the fusion {1,2,3}|{4} closes, but s = t is excluded for real octagons by the Feit-Higman condition that 2st be a perfect square.
- The Feit-Higman condition and the non-existence claim it supports are cited from memory and not re-checked against a primary source in this work.
