# Full lines in the local cycle test

The [rank-48 matching table](indexed_cycle_blocks.md) satisfies
the distance-block counts and the 179-element projector tests, but
it has no completion to a generalized octagon of order $(2,4)$.
This exclusion depends only on the distance-two matching entries.
The unspecified distance-three and distance-four positions can be
left free. Completing the lines also rules out 71 of the 81 labelled
attachment cases for an eight-cycle in the opposite-pair graph.

## What the full-line test uses

Fix roots $A,C$ at line distance two. The two rooted families are as
in the block test. On each side, name the root and the intermediate
lines of the groups being checked. For each such intermediate line,
name all eight leaf lines through its two outward points. Name all
three points of every line. Also name the common intermediate line
between $A$ and $C$.

Whenever a distance-two match between two leaves is given, name
its intermediate line and all three points on it. A matched leaf
uses its two outward points for its two matches, in the order of
the two branches on the other side. These attachments are distinct
because the root and the two matches are the three middle lines
of an opposite pair. The inward point is already used by the
path to the other root. Independently naming the unused point
of each added line permits coincident points until ruled out.

Use line distance domains for the root paths and the cross blocks.
On either side, the paths through the root give all distances
because they have incidence length at most eight. An intermediate
line on one side is opposite every intermediate line on the other
side by the eight-edge path through both roots. A leaf with its
two matches known has distance three to the opposite-side
intermediate lines over those matches, and distance four to the
others. Indeed a distance-three path continues through that
intermediate line to the opposite root in eight edges, so it
belongs to one of the known geodesics. The blocks containing a
match have the exact distances from the block test. If a leaf's
two matches are given, any other cross entry in its row is
restricted to distances three and four. The propagation applies
the existing [three-line constraints](grid_path_closure.md) to
these domains and the named point incidences.

## An impossible complete star

For the rank-48 construction, give the 128 matches and leave all
distance-three placements in the unmarked blocks unspecified.
The propagation keeps 423 point slots and 275 line slots. Let
$x$ be the row-side ancestor with branch, intermediate line and
child indices all zero. Take the match from the row with indices
$(1,0,0,0)$ to the column with indices $(1,0,1,0)$, and let $z$
be the third point of its intermediate line.

The replay forces $d(x,z)=3$ in point distance. It also forces
the ten distinct neighbors of $x$ on its five named lines all
to have distance four from $z$. A geodesic of length three
from $x$ to $z$ would start with a neighbor at distance two.
This contradiction eliminates the matching table itself.

There is a corresponding exclusion for a family of matching
tables. Keep the first three branch-to-branch bijections in the
abstract example and replace the fixed child flip in the fourth
by an arbitrary bit depending on $(u,p)$. Each zero bit gives a
four-cycle on its $(u,p)$ fiber since the return map fixes the
child there. If there are no zero bits, the result is the table
excluded by the complete star. Thus every octagon whose matching
table has this form contains a grid.

## Attachment patterns of an eight-cycle

By [the short-cycle result](opposite_cycles.md), an eight-cycle
in a non-grid component belongs to a system indexed by two lines.
It uses two leaf lines in each of the two branch classes on
each side. For each such pair of leaves, let its relation be
zero for different intermediate lines, one for the same
intermediate line and different ancestor points, and two for
the same ancestor point. These are all possibilities since the
two leaves are distinct.

There are $3^4=81$ labelled cases. Each specifies the positions
of the eight leaves in the rooted families up to relabelling
within the tree groups. Given those positions, name complete
groups for every intermediate line used and supply both known
matches of each cycle leaf. Leaves outside the cycle retain
unknown matches wherever allowed by the prescribed distances.
Two prescribed matches in one ancestor block violate the
one-entry rule. In 46 further cases the three-line propagation
finds an empty domain.

The ten surviving cases have only zero and one entries. They
consist of the all-zero case, the four choices with one entry
equal to one, the four choices with one entry equal to one on
each side, and the all-one case. In particular **the two leaves
in any one branch on an eight-cycle have different ancestor
points**. If both branch classes on one side use a common
intermediate line within their leaf pair, both classes on the
other side must do the same. The surviving cases are necessary
possibilities and do not assert that a cycle is realized.

Run [the exact replay](indexed_cycle_geometry.py) using
`python thin/indexed_cycle_geometry.py`. The exclusions leave
ten eight-cycle attachment cases and do not force a grid in an
arbitrary octagon.
