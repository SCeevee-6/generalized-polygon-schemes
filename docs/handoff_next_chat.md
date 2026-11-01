# Handoff for the next session (written 1 Nov 2026)

Continue the project SCeevee-6/generalized-polygon-schemes (GitHub, branch main). First read README.md,
drgfusion/drgfusion_note.md, search/search_note.md, search/fusions/fusion_note.md, fusion/near_polygon_note.md,
triple/triple_note.md and docs/literature_check_2026-10-31.md. Push with the Git Data API (blobs, tree with
base_tree, commit, PATCH refs/heads/main with force false; re-read the head right before the PATCH; never
force-push). Pushes can be stopped by the platform safety check; then ask me and do not try another route. Public text
must read as plain, direct prose written by me. Label every result as known, new (not found in the sources checked),
conditional, or unverified.

GOAL: get one result out that can be published, even a short one (arXiv math.CO note, then EJC / Discrete Math /
Australas. J. Combin. / Art Discrete Appl. Math.). Do a literature check before computing on each item and stop an
item as soon as it turns out to be known. Keep the chat output short.

## Where things stand

- Solid but known or elementary: the 1755-point scheme and GO(2,4) facts, the octagon fusion theorem, the commutative
  T-actions, near-polygon fusions (BCN 4.2.17 for d = 3; box scan for d = 4), triple intersection numbers for GO(2,4).
- New but small: no SRG/DRG for T on 12480 and 14976 points (primitive case of T complete with CRS 2020); only
  complete multipartite SRGs on 7020, 8775, 11700 and two 17550-point actions; the 12480 T-configuration has no
  proper fusion and every fusion of the 14976 one refines the rank-16 2F4(2) configuration.
- Best candidate for a paper: `drgfusion/`. Nine tabled DRG arrays have a forced strongly regular fusion whose
  parameters are known not to exist. Three are already on Brouwer's page. Six are not ruled out anywhere I looked:
  {14,10,3;1,5,12} (fused srg(50,21,4,12), absolute bound), {19,12,5;1,4,15} (srg(96,38,10,18), Degraer 2007), and,
  conditional on Kovaľ's Oct 2026 preprints (cited in Brouwer's SRG table as Koval26h, not on arXiv),
  {39,30,4;1,5,36}, {112,77,16;1,16,88}, {21,16,10,1;1,2,16,21}, {25,24,2,1;1,2,24,25}.
- Open gap: TomLib class 417 of 2F4(2)' (17550 points, rank 54, 37 units). CP-SAT ran 3 h with the complete
  multipartite unions excluded: no solution, no proof (UNKNOWN).
- `fusion/fusion_note.tex` has still not been compiled after the booktabs fix (I run it: `wsl sh
  /mnt/c/Users/shree/Downloads/gps_gap/compile_tex.sh`; read the .log from Downloads/gps_gap).

## Plan, in order

1. **Make drgfusion a paper ("Forced fusions and nonexistence of some distance-regular graphs").**
   a. Literature: check whether the six arrays are already ruled out: BCN's printed tables and additions
      (BCN-ac.ps.gz on Brouwer's DRG page), Vidali's thesis and sage-drg (`drg/nonex.py`, and whether
      `check_feasible(derived=2)` would catch them), papers citing Degraer (2007) and the 2-distance-set bound. Drop
      any that are known.
   b. Widen the net, all cheap with the existing code:
      - Fused SRGs with v > 1300: test them against general SRG conditions not implied by the DRG itself (absolute
        bound v <= f(f+3)/2, conference conditions, claw bound, mu = 1 results, Neumaier/Brouwer–Haemers bounds).
      - Three-class fusions of d >= 4 arrays: look them up in van Dam's tables of three-class schemes (check the
        allowlist; ask me for the domain) and apply the absolute bound in its s-distance-set form.
      - Spherical s-distance bound for every eigenvalue of every tabled array: if the cosines u_1(theta), ..., u_d(theta)
        take s distinct values other than 1, then v <= C(m+s-1, s) + C(m+s-2, s-1) (Delsarte–Goethals–Seidel). This is
        what kills {14,10,3;1,5,12}; a sweep may find more. Check first whether BCN or Brouwer's feasibility checks
        already apply exactly this.
      - The extended tables (IA4a, IA5a, d >= 6) and arrays of diameter 8 and above, using closure of single unions
        rather than all set partitions.
   c. For the Kovaľ-dependent cases, write to me so I can ask Kovaľ or Brouwer for the preprints. Do not present them
      as unconditional.
   d. Write `drgfusion/nonexistence_note.tex` (3–6 pages: method, table of ruled-out arrays with reasons, code
      availability). Keep it honest: method known, these arrays new.

2. **Close class 417 (17550 points).** In order of cost: (i) look for a reduction through the block systems (check in the
   table of marks whether H, of order 1024, lies in the flag stabiliser, class 423 of order 2048; if so the 17550
   points sit over the 8775 flags in blocks of size 2; the closed subsets found by `maximal_closed_subsets` show the
   block systems); the block partition of a
   T-invariant graph need not be equitable, so work out what it actually forces before relying on it;
   (ii) split the CP-SAT model by valency k (k is a sum of unit sizes; each k gives a small model) and use 8 workers;
   (iii) symmetry breaking with N_T(H)/H acting on the orbitals (compute in GAP). If it closes, the statement "up to
   17550 vertices the only non-trivial SRGs and DRGs with a transitive T-action are the four from CRS 2020" becomes
   complete and can be a section of the note in item 1 or a short separate note.

3. **A real discovery attempt: new SRGs on open parameter sets.** The SRG/DRG search code works for any transitive
   action. Literature check first: which groups have already had a full Crnković-style search (Crnković, Švob,
   Rukavina, Mikulić Crnković, Maksimović and others: G2(4), G2(5), O(7,3), T, Mathieu groups, U(3,3), L(2,q),
   O±(6,2), A8, ...) and whether primitive groups of small degree have been searched systematically. Then pick groups
   not yet covered (candidates: G2(3), 3D4(2), Sz(8), U4(3), U5(2), S6(2), O8±(2), L3(4), J2 imprimitive actions) and
   target the open (yellow) parameter sets of Brouwer's SRG table with v <= 1300: only run actions whose degree equals
   an open v. A graph on an open parameter set is publishable on its own.

4. Housekeeping: compile fusion_note.tex after my WSL run and fix errors; push each milestone with code and outputs.

## Practical notes (this host)

- Windows, PowerShell only; `python` is not on PATH in powershell. Run scripts from the python kernel with
  subprocess, or with `background: true` powershell calls using `$w\.venv\python\Scripts\python.exe`. Processes
  started with Start-Process are killed when the call ends. multiprocessing pipes are blocked; use shards in
  separate subprocesses (see `search/fusion_closure.py`).
- OR-Tools CP-SAT fails to load inside the python kernel after numpy-heavy imports (WinError 127); run it in a
  subprocess that imports ortools first.
- The file system is case-insensitive: Brouwer's IA4A and IA4a collide (we use the suffix _ext).
- GAP 4.16.1 via host.mcp("gap","run_gap",code=...) in repl: about 60 s per call, no state, files in
  /cygdrive/c/Users/shree/Downloads/gps_gap/; results can also be read from Downloads/gps_gap/queue/*.result if
  the connector times out. Pre-bind globals used inside functions; end long lines with ;;. Report the exact GAP code
  with its output.
- aeb.win.tue.nl and www.win.tue.nl (Brouwer's tables), brauer.maths.qmul.ac.uk (ATLAS), arXiv, Crossref and
  OpenAlex are reachable. `git clone` fails (schannel); use the GitHub API.
- About 5 GB RAM free; coset_schemes.py handles degree 17550 in about 20 s and 1.5 GB.
