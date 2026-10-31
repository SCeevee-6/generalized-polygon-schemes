# Forced fusions in Brouwer's tables of feasible DRG parameters

Labels: **known** (in the literature or on Brouwer's pages), **not found** (not ruled out in the sources I
checked; the check is a search, not a full review), **conditional** (depends on an unrefereed preprint).

## What the script does

`drg_table_fusions.py` takes every intersection array in Brouwer's machine-readable tables
(https://aeb.win.tue.nl/drg/drgtables.html: IA3P, IA4A/B/C/P, IA5*, IA6*, IA7*, and the extended IA4a, IA5a;
arrays of diameter above 7 are skipped), computes the p^k_ij exactly from the array, and tests every set
partition of {1,...,d} for closure. Closure depends only on the array, so a closed fusion holds in every graph
with that array. 18508 closed proper fusions come out, most of them the obvious ones of bipartite or antipodal
arrays. `fetch_and_crosscheck.py` then looks up every two-class fusion (a strongly regular graph) with v <= 1300
in Brouwer's SRG table and keeps those marked as not existing. Output: `table_fusions.csv`,
`srg_crosscheck.csv` (3614 rows), `nonexistent_fusions.csv` (26 rows, 9 arrays). The fused SRG parameters
of the nine arrays were checked a second time from the eigenvalues of the array (the fused graph has exactly
two non-principal eigenvalues).

## The method

**known.** If a distance-regular graph has a strongly regular distance-i or distance-{i,j} graph whose
parameters are known not to exist, the array is ruled out. Brouwer's DRG page already lists results of
exactly this kind (Lang: {44,35,3;1,5,42}; Galazidis: {56,45,12,1;1,12,45,56} and {36,35,27,6;1,9,30,36}),
BCN Section 4.2 treats strongly regular distance graphs, and Vidali's `sage-drg` package (jaanos/sage-drg,
`all_fusions`, `check_feasible(derived=...)`) checks the parameters of fusion schemes recursively. Krein
conditions and the absolute bound of a fused scheme follow from those of the original scheme, so the method
only gains something through results specific to the fused object (SRG nonexistence results, or the
2-distance-set bound v <= f(f+3)/2, which is stronger than the scheme absolute bound).

## Result of the sweep

| array | v | closed fusion | fused SRG | why the SRG does not exist | status of the array |
|---|---|---|---|---|---|
| {44,35,3; 1,5,42} | 375 | {1,2}\|{3} | (375,22,5,1) | mu = 1 | **known** (Lang, on Brouwer's page) |
| {56,45,12,1; 1,12,45,56} | 324 | {1,2}\|{3,4} | (324,57,0,12) | Gavrilyuk–Makhnev, Kaski–Östergård | **known** (Galazidis) |
| {36,35,27,6; 1,9,30,36} | 324 | {1,4}\|{2,3} | (324,57,0,12) | same | **known** (Galazidis) |
| {14,10,3; 1,5,12} | 50 | {1,3}\|{2} | (50,21,4,12) | absolute bound: eigenvalue -9 has multiplicity 7, and 50 > 7·10/2 | **not found** |
| {19,12,5; 1,4,15} | 96 | {1,3}\|{2} | (96,38,10,18) | Degraer (2007) | **not found** |
| {39,30,4; 1,5,36} | 300 | {1,2}\|{3} | (300,26,4,2) | Kovaľ, Oct 2026 preprint | **not found**, **conditional** |
| {112,77,16; 1,16,88} | 750 | {1,3}\|{2} | (750,210,55,60) | Kovaľ, Oct 2026 preprint | **not found**, **conditional** |
| {21,16,10,1; 1,2,16,21} | 300 | {1,4}\|{2,3} | (300,26,4,2) | Kovaľ, Oct 2026 preprint | **not found**, **conditional** |
| {25,24,2,1; 1,2,24,25} | 352 | {1,4}\|{2,3} and {1,2}\|{3,4} | (352,26,0,2) | Kovaľ, Oct 2026 preprint | **not found**, **conditional** |

"Not found" means: the array is still in Brouwer's current table, it is not in the list of ruled-out arrays on
his DRG page, and it is not in the `sporadic` nonexistence list of `sage-drg` (`drg/nonex.py`). I have not
checked BCN's printed tables, Vidali's thesis, or the papers that cite Degraer (2007).

The Kovaľ references (Brouwer's SRG table, anchor Koval26h: "The Bannai–Munemasa–Venkov lattice method for
strongly regular graphs: a uniform discriminant-form condition and new nonexistence results", preprint,
Oct 2026) are not on arXiv as far as I can tell, so I could not read them. The four arrays that depend on them
stand or fall with those preprints.

The case {14,10,3; 1,5,12} is worth a second look by hand. The eigenvalue -6 of the DRG has multiplicity 7, and
the fusion says that in the 7-dimensional spherical representation for this eigenvalue the cosines of
distances 1 and 3 are equal, so the 50 vertices form a 2-distance set in R^7, which has at most 35 points. If
this is right, the array should already fail a standard check, which makes me suspect it is ruled out somewhere
I have not looked.

## Assessment

Six arrays ruled out by a known method, four of them only conditionally. This is the kind of item Brouwer adds
to his page ("X notes that the array ... is ruled out since ..."), not a paper on its own. The sweep only covers
v <= 1300 for the SRG lookup; fused SRGs on more vertices were not checked against anything.

## Reproducing

    cd drgfusion
    python fetch_and_crosscheck.py      # downloads the tables into cache/ (not committed), about 3 minutes
