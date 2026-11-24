# A full-size incidence structure with nine octagon roots

There is a connected incidence structure with 1,755 points, 2,925 lines,
three points per line and five lines per point in which a specified
three-by-three grid has no centered transversal. From **each of its nine
tips**, every line is reached by exactly one nonbacktracking incidence
path of length at most seven. The incidence girth is ten. Thus the nine
rooted path conditions, even imposed simultaneously in a single structure,
do not force a centered transversal.

The certificate lists the three zero-based point labels on every line.
The nine grid tips are points 0 through 8 in row order. The first 81
points and the 51 prescribed lines agree with the direction word
`000000000` in [the grid construction](opposite_seed.md). All six row
and column centers are at line distance two across the two families and
four within either family. Distances from the nine tips to all 2,925
lines give exactly the six row and column triples with a center. All
six transversals fail. [Line triples](rooted_grid_lines.txt).

## Verification

[check_rooted_grid.py](check_rooted_grid.py) takes the line triples as
input. It checks degrees, linearity, connectedness, all 51 prescribed
incidences, the tip profiles, and the centered triples by breadth-first
distances. It checks the full intersection table of each pair of distinct
tips and all cubic cosine sums on the nine tips, including repeated
indices, using their distances to every point.

With $B$ the point-by-line incidence matrix, nonbacktracking paths starting
at points have count matrices $W_0=I$, $W_1=B$ and $W_2=BB^T-5I$. At
length $j\ge3$, extend $W_{j-1}$ by the appropriate incidence matrix and
subtract $2W_{j-2}$ for odd $j$ or $4W_{j-2}$ for even $j$. These are the
numbers of ways to extend the old path one step before immediately reversing.
The verifier checks entries of $W_0+W_2+W_4$ and $W_1+W_3$ are at most one.
Bipartiteness then gives girth at least ten. At this girth two different
length-five paths with the same endpoints form a ten-cycle. Each ten-cycle
has five opposite point-line pairs, giving the exact count below.

Every row of $W_1+W_3+W_5+W_7$ totals $5+40+320+2560=2925$. Thus a row
has all entries one exactly when every line is reached within seven.
Checking all rows identifies the nine marked roots and measures the
failure elsewhere.

```text
$ python thin/check_rooted_grid.py
points 1755 lines 2925 degrees 5 3
marked shells [1, 10, 80, 640, 1024] required seed lines 51
centered rows and columns 6 centered transversals 0
point roots with a unique path of length at most seven to every line 9
incidence girth 10 cycles of length ten 2518 point diameter 5
incidence diameter 11 point-line pairs at distance above seven 1702775
```

## How the certificate was found

For a point $p$, use the nine distances to the marked tips as its profile.
For every prospective line, each coordinate of its three profiles must
have the form $(h,h+1,h+1)$ in some order. A point at tip distance one,
two or three must have exactly one incident line on which that coordinate
decreases. These conditions allow the search to assign a multiset of five
line incidences to each point before choosing individual incidences.

An integer program chose profile and line multiplicities invariant under
permutations of rows and columns. It included the 51 required lines, the
tip-pair marginals, the cubic cosine equations, and the partitions of
inward directions among the five lines through each point. Matching slots
of the same point profile and the same set of inward coordinates gave
the incidence structure. Swaps between such slots preserved all the
root conditions while removing cycles of lengths four, six and eight.
The published triples give a direct exact verification independent of
those search choices. The example shows why an argument using only the
nine marked roots needs an additional constraint elsewhere in the incidence
graph. A generalized octagon would require girth sixteen and incidence
diameter eight throughout.
