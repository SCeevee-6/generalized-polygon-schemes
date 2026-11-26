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

## Reading comparisons from the cubic edge labels

There is a geometric interpretation of the entire comparison word.
Number the vertices on an indexed cycle cyclically. Let $e_i$ be the
label of the edge $i,i+1$ and $s_i$ the comparison of the leaves at
$i,i+4$. In every grid-free component,
\[
 s_i=1\quad\Longleftrightarrow\quad e_{i+1}=e_{i+2}=4.
\]

To see this, write the middle lines at vertex $i+2$ as $C,D,E$, with
$C$ retained by both cycle neighbors. In the middle triples of those
two neighbors, the third lines beyond the endpoints of vertex $i+2$
are the leaves at $i$ and $i+4$. They are the middle lines on the
geodesics from the unused tip of $C$ to the unused tips of $D,E$,
respectively. Both lie three incidence steps from the tip of $C$.
If the two geodesics have different first lines, the route between
these middle lines through that tip has length six and is geodesic.
If they have the same first line, the 16-cycle shape shares exactly
that first edge. The route between the middle lines then has length
four. The latter case is exactly the case in which both retained
edges have label 4. [Failed edge labels](opposite_edge_labels.md)

The two leaves are in the same rooted branch. There their line
distance is three for different intermediate lines, two for different
ancestors on the same intermediate line, and one for the same ancestor.
This proves the formula, excludes comparison 2, and gives a geometric
proof that 101 is impossible. The [metric verifier](failed_pair_types.py)
also checks the two line distances in the failed metric orbits.

Each 16-shape vertex has one pair of incident label-4 edges. Their
common middle line together with either endpoint gives an index
containing that turn. An 18-shape vertex supplies no such turn. Thus
the total number of comparisons equal to one, counting each position
of each unoriented indexed cycle once, is **exactly $2N_{16}$**. In
particular, a component of 18-shape vertices permits only the
all-distinct-intermediate-line attachment pattern at lengths eight and
twelve. Any shared comparison forces 16-shape vertices at the three
successive positions $i+1,i+2,i+3$.

## A restriction on every longer cycle

On a cycle longer than eight, the four comparisons
$s_i,s_{i+1},s_{i+4},s_{i+5}$ cannot all be one. Otherwise the leaf at
$i+8$ uses the same ancestor as the leaf at $i$: two consecutive
changes between the two available ancestors on an intermediate line
return to the starting one. The same holds at $i+1,i+9$. The edges at
$i$ and $i+8$ would give two matches in the same ancestor block, which
has at most one. They would be the same edge. [Ancestor blocks](indexed_cycle_blocks.md)

Consequently the seven-edge pattern $444\,{*}\,444$ is impossible on
any indexed cycle of length at least twelve, for either value of the
middle label. **At least one fifth of its edges have label 3.**
Indeed every run of 4s has length at most six. If a single 3 separates
two runs, one of those runs has length at most two, so their lengths
sum to at most eight. If $g\ge2$ consecutive 3s separate the runs, their
lengths sum to at most $12\le8g$. Summing over the gaps bounds the
number of 4s by four times the number of 3s. The all-3 case also
satisfies the bound.

For length twelve, each color has at most one comparison equal to one.
If $b$ of the twelve edges have label 4, at least $2b-12$ adjacent
pairs of edges both have label 4. There are at most four such pairs
by the comparison formula. Thus **each indexed twelve-cycle contains
at least four edges labelled 3**. These necessary conditions leave
mixed-label cycles possible.

## Joint twelve-cycle metric replay

The 136 labelled comparison words at length twelve give 13 classes
under shift and reversal. For each class the verifier supplies one
assignment of the 32 labelled failed metrics to **all twelve vertices
at once**. Each vertex already names its three middle lines. Two of
the tip-to-tip geodesics use leaves two positions away on the cycle.
For the third, add two internal points and the three line spares,
allowing coincidences with other named slots. Impose each selected
36-slot metric and apply three-line propagation to the whole diagram.

All 13 assignments pass this check. Their scope is a single indexed
cycle and its named lines. Existence of an ambient completion and
compatibility with the other indexed systems remain open.

Run `python thin/indexed_twelve_metrics.py`.
[Assignments and verifier](indexed_twelve_metrics.py)

## Four systems through a 4-edge

Let a 4-edge join opposite pairs $U=\{A,B\}$ and $V=\{D,E\}$.
The four indices containing it are $\{A,D\},\{A,E\},
\{B,D\},\{B,E\}$. At $U$, the two indices using the middle line
retained by the other 4-edge make a 44 turn; the other two turn onto
the 3-edge. At $V$, the same statement holds with the roles of the
endpoint pairs interchanged. Since these are independent choices, the
four indices have the following exact distribution:

| Turn at $U$ | Turn at $V$ | Number of indices |
|---|---|---:|
| 44 | 44 | 1 |
| 44 | 43 | 1 |
| 43 | 44 | 1 |
| 43 | 43 | 1 |

Thus exactly one indexed cycle has the edge flanked by 3-edges, two
have a 44 turn at exactly one end, and one has a 44 turn at both
ends. In particular, an all-4 indexed cycle through this edge can only
use the last index. An indexed 8-cycle is not required to be all-4,
so this observation alone does not force the other indices to have long
cycles. Counting turns over edges also counts each 44 turn twice; it is
consistent with the existing total of $2N_{16}$ turns.
