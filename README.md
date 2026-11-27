# GO(2,4) uniqueness

This branch collects the repository's work toward deciding whether the
generalized octagon of order $(2,4)$ is unique. The general uniqueness
problem remains open.

## Status

A conditional theorem is proved by computer-assisted finite case analysis.
Every nine-tip grid in a generalized octagon of order $(2,4)$ has a centered
transversal. That configuration gives a full thin suboctagon of order
$(2,1)$, and De Bruyn's theorem identifies the ambient octagon as Ree–Tits.

The remaining step on this route is to force a grid in an arbitrary
generalized octagon. The repository also studies opposite components,
automorphism groups, triples, and incidence constraints that bear on the
same uniqueness problem. None closes the unrestricted question.

## Main proof files

- At collinear roots, outer components of sizes $256i,256j$ intersect in exactly
  $32ij$ points. All three component labels on a random opposite line are
  independent. [Collinear component counts](octagon/collinear_components.md)

- `thin/opposite_seed.md` derives the 21-point configuration from two
  opposite lines and explains the centered-transversal criterion.
- `thin/grid_path_closure.md` proves the conditional theorem and records
  the finite classification of grid directions.
- `thin/grid_path_closure.py` and `thin/grid_path_certificates.json`
  replay the direction-class computation.
- `thin/failed_pair_types.md` and `thin/failed_pair_types.py` enumerate
  the five local metrics when an opposite pair fails to give a grid.
- `thin/fixed_involution.md` constrains fixed structures of involutions in
  a hypothetical grid-free octagon.
- In a grid-free octagon every nonidentity involution fixes exactly the
  radius-four ball about a point, with 91 points and 45 lines.
  The earlier parity and quotient counts apply to this remaining case.
  [Central involutions](thin/central_involution.md)
- A repeat among eight leaf-map outputs forces a grid without an
  automorphism hypothesis. [Repeat test](thin/opposite_leaf_maps.md)
- Outside grid components the opposite-pair graph has girth at least eight.
  The nine distinguished points of a nine-cycle have no centered triple.
  [Short cycles and their line labels](thin/opposite_cycles.md),
  [nine-cycle simplex](thin/nine_cycle_simplex.md)
- In a grid-free octagon every opposite-pair component contains an
  indexed cycle of length at least twelve. At least 5,850 distinct
  indexed systems must contain such a cycle. Every assignment on the
  24-cell test tree has at least three reference residues mapping into
  large rooted pieces. [Three marked residues](thin/three_residue_marks.md),
  [two-mark count](thin/all_residue_marks.md),
  [cube-diagonal bound](thin/rooted_cube_diagonals.md)
- Spectral design rigidity, positive semidefiniteness on the twelve-cycle
  witnesses and global GF(2) ranks do not force a grid. The cycle space of
  the incidence graph gives a rank bound on the pair–apartment incidence.
  [Spectral and rank checks](thin/spectral_rank_checks.md)
- Every six-step non-backtracking walk with a shared endpoint line keeps
  that line at every even step. The sixth walk distribution and bounds at
  steps seven and eight still admit an all-18 projected moment witness.
  [Six-step endpoint returns](thin/nb_m6_cross.md),
  [non-backtracking moments](thin/nb_moment_sdp.md)
- Every indexed cycle of length at least twelve has at least one fifth
  of its edges labelled 3. At length twelve there must be at least four.
  [Cyclic comparisons and edge labels](thin/indexed_cycle_paths.md)
- A global complex couples the indexed systems at all root lines. Its
  boundary gives a binary rank bound and its real incidence matrices give
  exact projected moments. [Global compatibility](thin/opposite_boundary.md)
- In a grid-free octagon every line has at least 45 opposite partners
  whose 132-vertex auxiliary graph has two disjoint squares.
  [Opposite-root square count](thin/opposite_root_squares.md)
- A failed opposite pair of 16-cycle shape forces a rooted surface with
  at least fourteen faces. The attachment types through fourteen faces
  are constrained by exact local replays.
  [Small surfaces and the global count](thin/opposite_link_surfaces.md)
- `thin/rooted_grid.md` gives a girth-10 example showing why the nine
  rooted path conditions alone do not force a centered transversal.

Other grid constraints, exact marginal systems, and computational checks
are documented in the remaining files under `thin/`.

## Related uniqueness work

- `octagon/status.md` summarizes recognition and restrictions on opposite
  components.
- `octagon/` contains the associated proofs and verification programs.
- `triple/` studies triple intersection data and non-opposite sets.
- `tits1755/` and `data/` contain the Ree–Tits model, generators, and
  related association-scheme computations.

## Reproducing the central check

From the repository root, run

```sh
python thin/grid_path_closure.py
```

This enumerates the direction classes and replays their certificates.
The scripts and notes describe any additional dependencies and model data
required for their checks.

## Branch contents

The branch was cut from `main` at
`8adffaa3d12d046492327011f4b9a789967f2bd6`. It retains the repository files so
that cited computations and data are available together.
