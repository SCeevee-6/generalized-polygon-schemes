# Handoff prompt for the next chat (written 2026-10-31, after commit on top of 9961f7d)

Paste everything below the line into a new chat.

---

Continue the project SCeevee-6/generalized-polygon-schemes (GitHub, branch main). Read README.md,
tits1755/coset_fusions_note.md, fusion/certify_note.md and docs/literature_check_2026-10-31.md first.
Push with the Git Data API (blobs, tree with base_tree, commit, PATCH refs/heads/main with force false;
re-read the head right before the PATCH; never force-push). Public text must read as plain, direct prose
written by me. Label every result as known, new (not found in the sources checked), or unverified.

GOAL: find one result that is actually new and big enough for a short paper (EJC, Discrete Math, JACo,
Des. Codes Cryptogr., or an arXiv math.CO note). So far everything solid is a reproduction or an
elementary fact. Do a literature check before computing anything for each angle, and stop an angle as
soon as it turns out to be known.

WHAT EXISTS (all in the repo, all checked)
- 1755-point scheme of T = 2F4(2)' = GO(2,4) distance scheme; array {10,8,8,8;1,1,1,5}; known (Crnkovic-
  Rukavina-Svob, JACo 51 (2020)). Not Q-polynomial (De Bruyn-Vanhove 2015), not tight (JKT 2000).
- Exact fusion theorem (fusion/certify_fusions.py, Groebner over QQ from the intersection array): for
  generalized hexagons/octagons with s,t >= 2, a proper fusion exists iff s = t; octagons never (Feit-Higman).
  Hexagon part known; octagon part not found in the literature, but elementary.
- Coset schemes of T on all 8 maximal subgroups (tits1755/coset_schemes.py; GAP coset actions in
  tits1755/cosets/). 1600: SRG(1600,351,94,72), known. 2304: three fusions, none P/Q-polynomial. 2925: GO(4,2).
  12480 (A6.2^2, rank 24) and 14976 (5^2:4A4, rank 25): non-commutative, centre dim 13. Their fusions were
  NOT searched.
- Breuer-Lux database (ctbllib 1.3.11) does not contain the Tits group.

TOOLS
- GAP 4.16.1 through host.mcp("gap","run_gap",code=...) in the repl tool: about 60 s per call, no state
  between calls, files in /cygdrive/c/Users/shree/Downloads/gps_gap/ (ATLAS generators and max SLPs there).
  Pre-bind globals used inside functions; end long lines with ;;. Read results from stdout.
- Python kernel: numpy, scipy, sympy. coset_schemes.py handles degree 14976 in about 50 s (dense int8
  orbital matrix; memory ~1 GB). About 4 GB RAM free on this machine; degree above ~25000 needs a sparse
  rewrite.
- TeX compiles only on my side: wsl sh /mnt/c/Users/shree/Downloads/gps_gap/compile_tex.sh.

ANGLES, in the order I think they should be tried (justify any change of order)

1. Strongly regular and distance-regular graphs from T and 2F4(2) beyond the Crnkovic-Rukavina-Svob range.
   Their Theorem 8 covers primitive T-actions up to 12000 vertices and imprimitive ones up to 7000. Extend:
   (a) primitive: degrees 12480 (two classes) and 14976 for T, and 14976 (Fi22N5, rank 16) for 2F4(2);
   (b) imprimitive: transitive T-actions on cosets of subgroups of index 7000 to ~20000 (list subgroup
   classes with GAP's table of marks if TomLib has 2F4(2)', otherwise via maximal subgroups of maxes).
   Method: structure constants from coset_schemes.py; enumerate all unions S of orbitals closed under
   transposition and test A_S^2 in span{I, A_S, J} (SRG) and distance-regularity of each symmetric union
   (layers of the distance partition must be unions of orbitals with constant intersection numbers).
   About 2^20 subsets per action, cheap with the structure constants. Then compare parameters with Brouwer's
   tables and decide isomorphism or automorphism group in GAP (GRAPE/nauty if available, else by rank of the
   automorphism group on the orbital graph). Output: a theorem of the form "up to N vertices, the only
   SRGs/DRGs with a transitive T-action are ...". New if the search range is new, even if the answer is empty.

2. Fusions as a nonexistence tool for distance-regular graphs. If the parameters of a feasible array force a
   fusion (closure is a property of the intersection numbers alone), then a graph with that array contains an
   SRG or a smaller scheme with determined parameters. If those parameters are known not to exist (Brouwer's
   SRG table, absolute bound on the fused scheme, Krein conditions on the fused scheme, known nonexistence
   results), the array is ruled out. Run the exact fusion test over every feasible array in Brouwer's DRG
   parameter tables for diameter 3 and 4 (download from his site if the domain is allowed; otherwise ask me).
   First check the literature: Brouwer-Van Maldeghem 'Strongly regular graphs' (2022) ch. on DRG with SRG
   distance graphs, and BCN 4.2.17-4.2.18. New nonexistence results for open arrays would be publishable.

3. Parametric fusion classification for other families with arrays polynomial in parameters: regular near
   hexagons and near octagons with parameters (s, t_2, t), dual polar graphs, the distance-2 graphs of
   bipartite doubles, antipodal covers. Use the certify_fusions.py approach (Groebner over QQ, exact real roots,
   integer scan). Known cases to reproduce first: Johnson and Hamming (Muzychuk 1988 and later), dual polar
   graphs (check). The near-polygon case with general (s, t_2, t) may be new and would contain our octagon
   theorem as a special case.

4. Non-commutative schemes of T (12480, 14976): find all commutative fusions (fusion to a commutative scheme,
   test closure and commutativity of the block sums). Report the lattice of fusions. Probably low impact but
   unexplored.

5. Stretch: triple intersection numbers with respect to the vanishing Krein parameter q^4_44 = 0 for
   extremal octagons GO(q,q^2) (Coolsaet-Jurisic method), aiming at structural restrictions or a computer
   proof of uniqueness of GO(2,4). Uniqueness of GO(2,4) is open; only attempt after a literature check of
   Coolsaet, Jurisic, Vidali, De Bruyn.

HOUSEKEEPING
- fusion/fusion_note.tex: compile result from my WSL run; fix errors if it still fails.
- Do not re-claim anything listed as known. Do not overstate. Keep the chat output short and push each
  milestone with its code and output files.
