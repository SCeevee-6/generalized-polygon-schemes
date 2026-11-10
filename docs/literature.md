# Literature search, 31 Oct 2026

Twelve queries, each run against Crossref (`query.bibliographic`, 8 hits) and OpenAlex (`search`, 8 hits). Raw hits are in `data/lit_search.csv` (186 rows). Queries:

fusion schemes of association schemes generalized polygon; Bannai Muzychuk criterion fusion association scheme; fusions of distance-regular graph association scheme; Q-polynomial generalized octagon; Q-polynomial generalized hexagon distance-regular; coherent configuration Tits group 2F4(2)'; multiplicity-free permutation characters Tits group; strongly regular graphs rank four Tits group; association scheme generalized octagon Ree-Tits; amorphic association schemes fusion distance-regular; Tits simple group distance-regular graph permutation representation; non-commutative association scheme sporadic group centralizer algebra.

## Found

- De Bruyn and Vanhove, On Q-polynomial regular near 2d-gons, Combinatorica (2015), doi 10.1007/s00493-014-3039-x. Covers the Q-polynomial question for thick octagons.
- De Bruyn, A combinatorial characterization of extremal generalized hexagons, doi 10.37236/9245; and Extremal generalized hexagons and triple intersection numbers, Discrete Math. (2021), doi 10.1016/j.disc.2020.112280.
- Breuer and Lux, Comm. Algebra 24 (1996), doi 10.1080/00927879608825701. The database shipped with ctbllib 1.3.11 has 26 sporadic groups and 12 automorphism groups and does not include the Tits group (`tits1755/gap/extra_checks_output.txt`).
- Crnkovic, Rukavina and Svob, J. Algebraic Combin. 51 (2020), doi 10.1007/s10801-019-00888-5: the 1755-point graph, and the SRG(1600,351,94,72) with automorphism group T (Theorem 8, Table 12).
- arXiv:2510.25410 (Strongly regular graphs of rank four) lists the same SRG(1600,351,94,72) for 2F4(2)' with stabiliser L3(3).2 as an exceptional rank-four example.
- General fusion literature, none of it on generalized polygons: Kharaghani et al., Non-commutative association schemes and their fusion association schemes (FFA 2018, doi 10.1016/j.ffa.2018.03.007); van Dam and coauthors on amorphic schemes and fusions of pairs (JCTA 2025, doi 10.1016/j.jcta.2025.106045); Herman et al. on isolating fusions (arXiv 2205.01810); Sankey, Weighted association schemes, fusions, and minimal coherent closures (doi 10.1007/s10801-014-0553-2).
- Ree-Tits octagon geometry: Coolsaet (embeddings), Van Maldeghem (Proc. LMS 1998, doi 10.1112/s0024611598000070), Belousov (automorphisms of GO(2,4)), Ihringer (partial ovoids).

## Not found

- Any paper on fusions of the distance scheme of a generalized octagon or hexagon.
- Any paper on the coherent configurations of T on cosets of A6.2^2 or 5^2:4A4.
- Any use of the Bannai-Muzychuk criterion on polygon schemes.

These are keyword searches. BCN (1989) and Van Maldeghem's book (1998) were not read page by page, so a miss here doesn't show the results are new.

## Later additions

- Crnkovic, Rukavina, Svob, arXiv:1809.10197 (= JACo 51 (2020)), read in full for Section 3.4: primitive T-actions up to 12000 points, imprimitive up to 7000; SRGs (1600,351,94,72) and (3510,693,180,126); DRGs GO(2,4), GO(4,2). OpenAlex lists three citing works (Borsuk counterexamples 2020; two Mathieu-group DRG papers 2021/2022); none extends the range for T.
- Fusions as a nonexistence tool are already used on Brouwer's DRG page (Lang; Galazidis; Huang via halved graphs) and implemented in Vidali's `sage-drg` (`all_fusions`, recursive `check_feasible`). So the method is known.
- Koval, SRG nonexistence preprints (Oct 2026), cited in Brouwer's SRG table as Koval26h; not found on arXiv.
- Diameter 4 (added 10 Nov 2026, web keyword searches). Fiol, A spectral characterization of strongly distance-regular graphs with diameter four, arXiv:1407.1392: Gamma_4 is strongly regular iff b_3 = a_4 + 1 and b_1 = b_3 c_3 (quoted from the van Dam-Koolen-Tanaka survey, Section 13.7; the paper itself was not read). arXiv:1608.00373 treats graphs whose distance-2 or distance-1-or-2 graph is strongly regular, with H(4,3) as an example. So the two-class fusions in `fusion/near_octagon_note.md` are special cases of known criteria. A classification of all fusions of the dual polar schemes was not found (searched: subschemes / fusion schemes of dual polar association schemes); the B_d/C_d merging of distances 1 and 2 is Ivanov-Muzychuk-Ustimenko.
- T(n) is determined by its parameters for n other than 8, so the SRG(45,16,8,4) from the thin octagon (2,1) is T(10); also checked directly (`thin/verify_t10.py`).
