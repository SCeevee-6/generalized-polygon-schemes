# Where GO(2,4) uniqueness stands, and the sibling obstruction

The unrestricted uniqueness question remains open. The new exclusions cover seven explicit bare graphs on 512 vertices and every
256-point trace component with 2-group monodromy, including nonregular covers.
A second end-to-end implementation also reproduces the regular 256 exclusion. A symmetry error
reopened the earlier cross-edge search; the corrected search is now complete and
leaves nothing, so among elementary abelian regular 512 covers only the Ree–Tits
component survives.

## An obstruction on the bare outer graph

Let C be any 512-point opposite component. Every shell-three point carries
exactly two C-lines by the boundary-count proof in the README. A C-line
contains two points of C, so it corresponds to one edge of the bare graph.
Two shell-three points on the same inward line are siblings. Their four
C-lines form a block (a,b;c,d), with a,b on one sibling and c,d on the other.
The 1,280 edges of C must be partitioned into 320 such blocks.

Subdivide each edge of C once and call its subdivision vertex a line vertex.
Write D for half the distance between line vertices in this graph. For distinct
edges D is one plus the minimum outer distance between endpoints. Inserting
the two sibling points and their common inward line adds a tree whose path has
half-length 1 within either pair and 2 between the pairs. Girth 16 requires

* D(a,b), D(c,d) ≥ 7,
* D(a,c), D(a,d), D(b,c), D(b,d) ≥ 6.

Enumerate every block satisfying these conditions, over **all** edges. Fix an
outer line e. For each possible block Q containing e, insert its tree.
Every line f outside Q still needs a disjoint block R. The same distance lower
bounds apply to R in the enlarged graph because the interior of its own tree
is absent. For any two lines of R, travelling through a pair of lines of Q
and its new tree supplies an upper bound on their distance. One short route
is enough to forbid R. Overlap of the two outer pieces of the route causes
no difficulty: a shortest path in the enlarged graph is at most that long,
and its interior is disjoint from R's new inward tree.

If some f has no possible R, discard Q. If this happens for every Q at e,
then the bare graph cannot occur. Different choices of Q may have different
witnesses f. The deck group, the trace map and the colouring play no role
after the input edge list has been built.

`python octagon/sibling_obstruction.py --verify` checks all 46 saved labelled
inputs. The last switch rechecks every reported witness by breadth-first
distances in the literal incidence graph with Q inserted. The archived
[sibling obstruction output](sibling_obstruction_output.txt) gives the anchor
and one uncovered line for each candidate, in the verifier's enumeration order.
Only zero-based row 6 survives. BLISS identifies this row with both actual
components. In the known model all true sibling blocks can of course coexist.

| First row in stage4_fp0.json | Labelled rows | Bare automorphism order | Diameter | Outcome |
|---:|---:|---:|---:|---|
| 0 | 19 | 128 | 7 | excluded |
| 1 | 5 | 3072 | 8 | excluded |
| 3 | 4 | 2048 | 7 | excluded |
| 5 | 3 | 2048 | 7 | excluded |
| 6 | 1 | 10240 | 6 | Ree–Tits component |
| 10 | 9 | 512 | 8 | excluded |
| 20 | 3 | 6144 | 7 | excluded |
| 26 | 2 | 4096 | 7 | excluded |

All exclusions were run on individual rows. The isomorphism calculation only
groups the results for presentation.

## What happened to the elementary cover enumeration

I rebuilt the distance tests on all 1,064 saved stage-three inputs. Every
category's count of compatible sibling pairs agreed with the original data.

| Marked orbit | Saved inputs | Rejected at stage 3 | Left | Left after bare obstruction |
|---|---:|---:|---:|---:|
| Same base edge | 912 | 866 | 46 | 1 |
| Distinct base edges, historical run | 152 | 152 | 0 | 0 |

In a tree gauge, a marked translate is measured between the smaller endpoints
of the two base edges. Under pullback by h, transport the two root-to-endpoint
tree paths by h and reverse any image edge whose reference endpoint changed.
The XOR of the chord voltages along these paths and edge reversals must be
added to the old translate. For the same edge the contributions cancel. For a
cross pair, only four of the 32 stabilizers in the old program preserved the
forbidden-walk masks with the translate held fixed. Thus lexicographic pruning
could compare a valid marked object to an inadmissible one. This invalidates
the old inference of exhaustive exclusion, irrespective of its leaf counts.

`python octagon/check_mate_action.py` compares the full forbidden-walk sets,
as identities in the 25 independent binary chord coordinates, and checks the
actual column-support arrays sent to the search kernel.

    case 0: stabilizers 48, old walk-set symmetries 48, corrected 48, nonzero shifts 0
    case 1: stabilizers 32, old walk-set symmetries 4, corrected 32, nonzero shifts 28

With the corrected action there are 15,086 depth-eight cross prefixes (old:
9,914). Prefix 0 completed with zero leaves. The same-edge run has 5,925
prefixes; prefix 0 again rejected 128 leaves at stage two. The original
10,779,630,523-node same-edge generation has not been repeated independently. Its effective symmetries are unchanged, and all 25 unmarked chord
voltages must be nonzero because their fundamental cycles have length 4 or 5.
The corrected cross sweep has been run over all 15,086 prefixes (split into four
half-open ranges, no prefix truncated). It visited 1,442,239,465 nodes. Of the
leaves, 3,889 failed the stage-two H-sibling test, 18 gave disconnected covers and
96 were left. Stage three rejects 95 of them. The bare-graph sibling obstruction,
run with the witness recheck on all 96 regardless of the stage-three verdict,
excludes every one, including the stage-three survivor (row 22). Voltages and
per-category data are in `stage3_fp1_switched.json`, the obstruction witnesses in
`sibling_stage3_fp1_switched.json`, totals and survivor sequences in
`elem_fp1_switched_summary.json`. Rerun with

    python octagon/elem_run.py 1 LO HI
    python octagon/stage3.py 1 elem_fp1_switched_*.json   # writes stage3_fp1.json, saved here as stage3_fp1_switched.json
    python octagon/sibling_check.py octagon/stage3_fp1_switched.json

Old cross shards cannot fill gaps in the new indexing.

A further exact pruning rule applies to the residual case in which all
geometrically used mates cross between distinct base edges. In each category
let nᵢⱼ be the number of pairs connecting slots i and j among its four base
edges. Each slot has 32 lifts, so every row sum is 32 and nᵢᵢ=0. Subtracting
pairs of row sums gives n₁₂=n₃₄, n₁₃=n₂₄ and n₁₄=n₂₃. A used cross pair
therefore forces a mate between the complementary slots. This does not
require the deck group to preserve the pairing in a completion.

## Recognition and the boundary count

Both components of the ATLAS model give the recorded outputs below.
`verify_recognition.py` now asserts the truth comparisons, the six level sizes,
and the thin witness. Each successful seed generates 45 points on 30 lines
in the forced half, with point-degree 2, line-degree 3, incidence girth 16
and diameter 8. The seeds use the ATLAS T labelling (zero-based).

    component 0 mates forced True pairing correct True six merges correct True GO(2,1) seed [0, 146, 314, 1256]
    component 1 mates forced True pairing correct True six merges correct True GO(2,1) seed [0, 491, 1579, 1663]

The merge sizes are 640, 320, 80, 40, 10, 5, 1. In any completion siblings
at each step must pass the corresponding girth test. The admissibility
graph is checked to consist of disjoint cliques of the required sizes.
Thus an isomorphic copy of the bare component in an arbitrary octagon
forces a half containing the same subgeometry. De Bruyn's Main Theorem
then identifies the octagon. [De Bruyn, Discrete Mathematics 338 (2015),
2125–2142](https://doi.org/10.1016/j.disc.2015.04.022)

The component-size restriction is a proof from the intersection array.
For a proper component, the centred indicator on the outer layer attains
the upper nonprincipal Rayleigh quotient 5 because all five outer neighbours
stay in the component. At every shell-three point the corresponding
eigen-equation gives exactly |C|/256 incident C-lines. For the single
component all four lines belong to it directly. Both model components
have four outer neighbours at each of the 640 boundary points, i.e. two
C-lines. Sizes 256, 512, 768 and 1024 follow without an automorphism assumption.

The five-cycle computation was rerun: 2,167,056 invariant covers, 251 of
girth at least eight; the mate-count histogram is 217 of type (0,0), 28 of
type (0,1), five of type (1,1) and one of type (5,5). Four of the five
one-mate graphs fail the first inward merge at every boundary point. For
the five-mate graph the ten categories each have 320 mate pairs and zero
compatible sibling pairs. The known graph is the unique survivor in this
symmetry class.

## Regular 256-point components

The exact original run from the abstract base produced

    19132 245 48 27979
    [((0, 'bipartite'), 40), ((0, 'eight'), 40), ((0, 'four'), 16), ((0, 'sixteen'), 40), ((1, 'eight'), 16), ((1, 'four'), 152), ((1, 'gap'), 136)]

The separate implementation in [regular256](regular256/README.md) starts
from the 16 base vertices, enumerates the double and order-four quotients, and
checks the final two lifts. Outputs agree record by record at the double and
four-sheet stages. The second route omits the cyclic spectral-gap filter
altogether: all 152 cyclic quotients lead to zero final graphs. All 40 final
records in the elementary branch are bipartite. The proof of the least
eigenvalue bound excludes each of them. The regular 256 result has two
computational implementations with the same mathematical reduction. Neither
implementation assumes that later double lifts are regular over the base.

## Nonregular 256 covers with 2-group monodromy

Removing regularity of the intermediate four-sheet quotient still permits a
complete calculation. For each of the 48 admissible first quotients we tried
all switching classes of double lifts through its 49 chord variables, subject
to edge-tree capacity. There are 632 capacity records. A rational residual
certificate for a signed adjacency eigenvalue rejects 566; the other 66 yield
106 eight-sheet and 64 sixteen-sheet records. All 64 are bipartite. The
[proof and scope](regular256/PROOF.md) show that this excludes every 256
cover with 2-group monodromy. Any surviving 256 component has an odd prime
divisor of its monodromy order. The degree alone imposes no contradiction.

## Other constraints

| Result | Check | Scope |
|---|---|---|
| 13-overlap graph | 24,960 copies, two components of 12,480; shells 1/45/1170/9300/1964; BLISS order 17,971,200. Transport of the two local rules from one edge under a checked transitive edge action gave 3,510 classes of 80, then 1,755 centres and exactly the 8,775 true collinear pairs. | Recognizes the known overlap graph. An arbitrary octagon must first supply thin suboctagons. |
| Maximum flag families | Numerical pullback dimensions 351, 675, 675, 351; rank 2052. Both canonical families have 135 flags, zero opposition pairs and occupancies with constants 5 and 3. An exact incidence-operator argument now proves the potential description for any GO(2,4). | Existence of a mixed maximum family remains undecided. See [flag note](../flags/max_family_note.md). |
| Dual packing | 133 distinct indices. Among unordered pairs, 3,187 have dual distance 3 and 5,591 have distance 4. The disjoint balls cover 1,729 vertices. | Lower bound in the known dual. The upper bound 223 uses published nonexistence of a perfect packing and the defect-one argument. |
| Smith | Point opposition (binary rank, multiplicity of valuation one) = (26,272); line opposition = (246,350). | Exact integer-valued sums in the ker-to-coker calculation for the known model. |
| Cayley trace | HiGHS returned INFEASIBLE for order-625 group IDs 3,4,6,7,8,9,10,12,13,14. The trace and multiplicity equations and the eight table-based arithmetic exclusions were rechecked. | These are tests on the separate diameter-three array {52,42,16;1,6,28}. Groups 12 and 13 in the trace note still use solver infeasibility. No octagon component is excluded by that computation. |

The class sizes and the packing histogram can be reproduced with `python octagon/check_model.py`. The ATLAS checks used files TF42G1-p1755B0.g1 and .g2 (SHA-256
8428995642391eb0e360492598c1e7fbcd17e01f335676ce88bb467661d5413f and
ff46d1bb2d9a6affeb9faa4ce4eaff6fb723745cc19ddecd601f2921d8d4f813).
The sorted-line checks use `data/generators.npy` (the ATLAS generators c, d of 2F4(2) as permutations of the 1755 points).

These runs used Python 3.12, NumPy 2.3.5, SciPy 1.17.0, python-sat 1.9.dev24,
igraph 1.0.0, NetworkX 3.7 and Numba 0.68.1. Commands from the repository root (the recognition check can also use
the stored point permutations):

    python octagon/verify_recognition.py
    python octagon/check_model.py
    python octagon/check_mate_action.py
    python octagon/sibling_obstruction.py --verify
    OPENBLAS_NUM_THREADS=1 python octagon/regular256/check_quotients.py
    python octagon/regular256/check_terminal.py
    OPENBLAS_NUM_THREADS=1 python octagon/regular256/double_towers.py
    OPENBLAS_NUM_THREADS=1 python octagon/regular256/double_towers.py --sat-terminal
    python flags/flag_structure.py data/generators.npy
    python smith/smith_two.py data/generators.npy
    python search/cayley/cayley_trace.py search/cayley/out_normal/tables_625.txt '[52,42,16]' '[1,6,28]' /tmp/trace625.json

## Logical route to full uniqueness

The shortest available sufficient statement is that **every** GO(2,4)
contains one GO(2,1). Existence at one root suffices. The component route
would establish this by forcing one recognized bare graph at some root.
The current universal proof allows the following patterns at a fixed root.

| Sizes | Remaining obstacle to recognizing a counterexample |
|---|---|
| 1024 | A single connected 64-sheeted trace cover. |
| 768 + 256 | A 256 cover whose monodromy has an odd prime divisor must coexist with a 768 cover. |
| 512 + 512 | A counterexample must escape recognition on both components. |
| 512 + 256 + 256 | Both 256 covers must escape the double-tower exclusion and the 512 cover must escape recognition. |
| 256 + 256 + 256 + 256 | Each of the four monodromy groups must have an odd prime divisor. |

The elementary classification is now complete: an elementary abelian regular
512 component is the Ree–Tits one. The next finite task is regular 512 covers
with the other 50 groups of order 32. Even a complete classification gives no reason for an
arbitrary octagon to have an elementary component, a 512 component, or
this component pattern at every root. A root-independent forcing theorem
for a thin suboctagon, or a component argument covering all these
possibilities and forcing a recognized root, remains necessary.

A direct alternative is to force one four-cycle in the cubic graph on
opposite-line pairs and one centered transversal in the resulting grid.
Both conditions together force a thin suboctagon without a component or
monodromy hypothesis. If no transversal is centered, 25 direction classes
survive exact tests on 81 points, and 8 survive the completed-star tests. The [opposite-line argument](../thin/opposite_seed.md)
gives the proof, model counts, integer witnesses and the remaining existence
question. [A girth-ten construction](../thin/rooted_grid.md) realizes all nine
marked roots in a single incidence structure of the full size.

A fresh check of the 96 stored corrected cross-orbit inputs reproduced the
archived sibling witnesses byte for byte. The final line was

    covers 96 left after stage 3 1 left after sibling obstruction []

This replay starts at the output of the generation sweep.
