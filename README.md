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
  [Short cycles and their line labels](thin/opposite_cycles.md)
- A global complex couples the indexed systems at all root lines. Its
  boundary gives a binary rank bound and its real incidence matrices give
  exact projected moments. [Global compatibility](thin/opposite_boundary.md)
- In a grid-free octagon every line has at least 32 opposite partners
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
