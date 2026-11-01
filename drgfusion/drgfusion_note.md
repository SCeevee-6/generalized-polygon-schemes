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
| {14,10,3; 1,5,12} | 50 | {1,3}\|{2} | (50,21,4,12) | absolute bound | **known**: the array itself fails the absolute bound (BCN Thm 2.3.3), q^h_33 = (7,1,1,0), 43 > 28 |
| {19,12,5; 1,4,15} | 96 | {1,3}\|{2} | (96,38,10,18) | Degraer (2007) | **known**: Neumaier, by counting 4-cliques (BCN additions, Ch. 4) |
| {39,30,4; 1,5,36} | 300 | {1,2}\|{3} | (300,26,4,2) | Kovaľ, Oct 2026 preprint | **conditional**, and an immediate corollary of that preprint |
| {112,77,16; 1,16,88} | 750 | {1,3}\|{2} | (750,210,55,60) | Kovaľ, Oct 2026 preprint | **known**: deleted in the BCN additions (Ch. 14, by BCN Thm 4.4.4) |
| {21,16,10,1; 1,2,16,21} | 300 | {1,4}\|{2,3} | (300,26,4,2) | Kovaľ, Oct 2026 preprint | **known**: the array fails the absolute bound, q^h_44 = (7,0,1,0,0), 43 > 28 |
| {25,24,2,1; 1,2,24,25} | 352 | {1,4}\|{2,3}, {1,2}\|{3,4} | (352,26,0,2) | Kovaľ, Oct 2026 preprint | **conditional**; the folded graph is a srg(176,25,0,4), which the same preprint also rules out, so no fusion is needed |

Exact Krein and absolute-bound numbers: `exact_check.py`, output in `exact_check_output.txt`.

## Literature check of 1 Nov 2026 (outcome: nothing new here)

The earlier "not found" labels were wrong for four of the six arrays. The reason is that Brouwer's machine-readable
files are the 1989 BCN tables of arrays that pass the elementary conditions, without the remarks column of the
printed book and without the later deletions. Being in those files says nothing about being open.

- `feasibility_sweep.py` runs the Krein conditions, the absolute bound and the Delsarte–Goethals–Seidel
  s-distance bound over all 5526 arrays in the files (including the extended IA4a, IA5a, IA6b–IA17b).
  44 arrays fail Krein or the absolute bound (`feasibility_sweep.csv`), among them {14,10,3;1,5,12} and
  {21,16,10,1;1,2,16,21}. These are presumably marked in the printed BCN tables; in any case they are ruled out
  by standard conditions, which `sage-drg` also checks.
- The s-distance bound (v <= C(m+s-1,s) + C(m+s-2,s-1) for the s distinct cosines of an eigenvalue of
  multiplicity m) rules out exactly one array, {14,10,3;1,5,12}, which the absolute bound already kills. It adds
  nothing on these tables.
- The BCN additions (BCN-ac.ps.gz, Brouwer's DRG page) delete {19,12,5;1,4,15} (Neumaier) and
  {112,77,16;1,16,88} (Thm 4.4.4).
- Fused SRGs on more than 1300 vertices: 718 primitive fused SRGs (240 with v > 1300, from 81 arrays) coming from
  arrays that pass Krein and the absolute bound. Krein conditions and the absolute bound of a fusion are
  inherited from the original scheme, so I tested the conditions that are not: Neumaier's claw bound (checked
  against Brouwer's SRG table: no existing graph violates my implementation, 24 nonexistent ones do) and the
  conference-graph condition. The only hit is Lang's {44,35,3;1,5,42}. The large triangle-free cases such as
  srg(3404,82,0,2) from {81,80,2,1;1,2,80,81} are open SRG parameter sets, so they give nothing either.
- Arrays of diameter 8 and more (30 arrays in all tables): every 2-class fusion found belongs to a classical
  family (hypercubes, Hamming graphs, doubled Odd graphs and their relatives), all of which exist.
- Not done: three-class fusions against van Dam's tables. Their Krein and absolute-bound conditions are inherited
  too, so only van Dam's case-specific results could help, and I stopped here.

So the item "forced fusions as a paper" is dropped: the method is known, two of the six arrays fail the absolute
bound, two are deleted in the BCN additions, and the last two follow at once from Kovaľ's preprint (one of them
without any fusion).

## Reproducing

    cd drgfusion
    python fetch_and_crosscheck.py      # downloads the tables into cache/ (not committed), about 3 minutes
    # IA6b ... IA17b are not fetched by the script; put them into cache/drg by hand as IA6b_ext.txt etc.
    python feasibility_sweep.py cache/drg feasibility_sweep.csv
    python exact_check.py "14,10,3;1,5,12" "21,16,10,1;1,2,16,21"
