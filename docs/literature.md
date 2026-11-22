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
- Vertex-primitive distance-regular graphs (added 20 Nov 2026, two web keyword searches). Distance-transitive graphs with a primitive group are nearly classified (introduction of 'Imprimitive distance-transitive graphs with primitive core of diameter at least 3', Michigan Math. J. 58 (2009)). GRAPE's `VertexTransitiveDRGs` handles one group at a time. A run over the whole primitive groups library against the BCN arrays was not found. Statuses of arrays: Vidali's sage-drg `drg/nonex.py` (read as text) and the DistanceRegular.org catalogue (www.math.mun.ca/distanceregular).
- 22 Nov 2026 (arXiv API, Crossref, OpenAlex and web keyword searches; every item below was checked to exist, and the reading level is given).
  - Gavrilyuk and Vidali, On local spectra of classical distance-regular graphs with b = 2, arXiv:2611.02986 (2 Nov 2026). Thm 2.6(a): no DRG with array {49,40,22;1,5,28}; statement and proof read. One of the 18 census arrays.
  - Chen, Makhnev and Klimin, On distance regular graphs with diameter 3 and degree 44, Proc. Inst. Math. NAS Belarus 32:1 (2024) 57-63. The abstract says {44,36,5;1,9,40}, {44,36,12;1,3,33} and {44,42,5;1,7,40} do not exist; full text not read. One of the 18 census arrays.
  - Cesarz and Williford, Four class association scheme double covers of strongly regular graphs, arXiv:2610.05420 (Oct 2026). Rules out the antipodal arrays {20,18,3,1;1,3,18,20}, {22,21,3,1;1,3,21,22}, {54,50,5,1;1,5,50,54}, {170,162,9,1;1,9,162,170} and a parity family. Antipodal, so outside the primitive census.
  - Fujisaki, Koolen and Tagami, Some properties of the twisted Grassmann graphs, Innov. Incidence Geom. 3 (2006) 81-87, doi 10.2140/iig.2006.3.81. Cited for twisted Grassmann graphs not being vertex-transitive; not read. Used only to say that a DRG need not be vertex-transitive.
  - Higman, Invariant relations, coherent configurations and generalized polygons, in Combinatorics (1975) 347-363, doi 10.1007/978-94-010-1826-5_18. Flag coherent configurations of generalized polygons; not read.
  - Leonard, Using Groebner bases to investigate flag algebras and association scheme fusion, in Algorithmic Algebraic Combinatorics and Groebner Bases (Springer 2009) 113-135, doi 10.1007/978-3-642-01960-9_3. Flag algebras of generalized polygons and their fusions; not read (described in the next item).
  - Colangelo, Monzillo and Siciliano, The association scheme on the set of flags of a finite generalized quadrangle, Discrete Math. 347 (2024) 114054, doi 10.1016/j.disc.2024.114054 (arXiv:2406.03942). Intersection numbers and all fusions for GQ flags; abstract and Remark 4.1 read.
  - Bamberg and Lansdown, Implications of vanishing Krein parameters on Delsarte designs, with applications in finite geometry, Algebr. Comb. 6 (2023), doi 10.5802/alco.246 (arXiv:2107.05207). Section 8 read: no nontrivial m-ovoids in a GO(s,s^2), via the same Krein zero (their q^2_22 = q^4_44 here, `triple/krein_convention.py`).
  - Cohen, O'Brien and Shpectorov, On the uniqueness of the generalized octagon of order (2,4), J. Algebra 421 (2015) 369-393, doi 10.1016/j.jalgebra.2014.08.034; De Bruyn, The uniqueness of a certain generalized octagon of order (2,4), Discrete Math. 338 (2015) 2125-2142, doi 10.1016/j.disc.2015.04.022. Conditional uniqueness results; not read beyond abstracts.
  - Petit and Van Maldeghem, The natural embedding of the Ree-Tits octagon, arXiv:2607.18560 (July 2026). Introduction read: pseudo-isometric embeddings match distances 1-4 with collinear, symplectic, special and opposite in the F4 metasymplectic space.
  - Buekenhout and Van Maldeghem, Remarks on finite generalized hexagons and octagons with a point-transitive automorphism group, in Finite Geometries and Combinatorics (1993) 89-102, doi 10.1017/cbo9780511526336.011. Not read.
  - Coolsaet and Van Maldeghem, J. Algebraic Combin. 12 (2000) 107-113, doi 10.1023/a:1026579608364: partial ovoids of GO(2,4) have at most 27 points.
  - Not found: a source for 'the distance-1,2,3 graph of GO(2,4) has exactly 2925 cliques of size 27' (now checked in `triple/line_cliques.py` and `triple/clique_certificate.py`, and proved for every finite GO of order (s,s^2) in `triple/max_nonopposite_proof.md`), and any explicit treatment of the 2F4(2)'-orbital configuration on the 8775 flags.
