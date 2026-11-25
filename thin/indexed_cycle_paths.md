# Attachment restrictions along indexed paths

For two lines at line distance two, the
[128-vertex indexed system](opposite_cycles.md) is a union of
cycles. A cycle of length four is a grid. The following restrictions
apply along every longer cycle, regardless of its length.

The vertices of an indexed system alternate between its two
sides, and the branch on either side alternates between its two
values. Give a vertex the color consisting of its side and its
branch. The colors occur in a fixed order of period four. For
two distinct vertices four steps apart, compare their leaves
in that color. Use 0 for different intermediate lines, 1 for
the same intermediate line and different ancestor points, and
2 for the same ancestor point.

**Consecutive visits to a color always have distinct ancestor
points.** In particular, each permutation on 32 leaves obtained
by taking four steps has a zero $4\times4$ diagonal block for
each ancestor, except for fixed points corresponding to grids.
More explicitly, each of these blocks can contain only its
diagonal entries.

**The cyclic comparison word has no occurrence of 101.** This
includes occurrences crossing the choice of starting vertex.

## Replay on open paths

The [replay](indexed_cycle_paths.py) builds open paths of four
and six edges. It uses the same full lines and ancestor blocks
as the [eight-cycle attachment test](indexed_cycle_geometry.md).
Only one match is supplied at each end of the path. A supplied
match excludes other matches in its target branch. Unknown
matches in the other branch are retained.

A path of four edges has two distinct leaves in just one of the
four colors. There are three ways to place them in the ancestor
tree up to relabelling. The case with the same ancestor is
excluded by three-line propagation. A path of six edges gives
three consecutive comparisons in the traversal order, so the
replay tests all $3^3$ possibilities. The cases with any 2 and
the case 101 are excluded. All other binary words pass this
local propagation. Applying these short path checks at every
starting vertex gives both restrictions.

## Indexed twelve-cycles

There are three leaves in each of the four colors on a
twelve-cycle. Their ancestor points are pairwise distinct,
since the three cyclically successive comparisons account for
every pair in one color. An intermediate line has only two
available ancestor points. Thus either the three intermediate
lines are distinct, or exactly one pair of leaves shares an
intermediate line. Each color has four choices, giving $4^4$
labelled combinations before the comparison-word condition.
Exactly 136 have no cyclic occurrence of 101. This reduces
the indexed twelve-cycle cases from 256 to 136. The count
is of necessary attachment patterns.

Run the replay from the repository root with
`python thin/indexed_cycle_paths.py`.
