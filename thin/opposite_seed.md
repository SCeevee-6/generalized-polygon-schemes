# Forcing a thin suboctagon from two opposite lines

Let the ambient geometry be any generalized octagon of order $(2,4)$. Start
with opposite lines $A,B$. There is a configuration on 21 specified points
without making any assumption on automorphisms or on opposite components.
It gives the following sufficient condition for uniqueness.

1. Six specified points of the configuration have mutual point distance three.
   If they have even one common point at distance three, they have exactly
   three, and the resulting nine points have a unique arrangement over the
   prescribed lines in a three-by-three array.
2. Take one point from each row and each column of that array. If some such
   triple has a line at incidence distance three from all its points, the
   octagon contains a full suboctagon of order $(2,1)$.

The second condition has six possible witnesses. Either all six work or none
does. A single opposite pair satisfying both conditions suffices for De
Bruyn's Main Theorem. The outstanding existence question concerns **one
pair in an arbitrary octagon**. The exact counts below hold for the Ree–Tits
model. [De Bruyn, Discrete Mathematics 338 (2015), 2125–2142](https://doi.org/10.1016/j.disc.2015.04.022).

**Subsequent result.** The [short-path refutations](grid_path_closure.md) prove that condition 2 always holds once the grid of condition 1 exists. It now suffices to find one opposite-line pair satisfying condition 1 in an arbitrary octagon. The case counts below record successive earlier tests.

## Distances and two elementary facts

Write $d$ for point-graph distance and use incidence distance for a point
and a line. Two lines at line-graph distance $j$ have incidence distance
$2j$. The incidence graph has girth 16 and diameter eight. A nonbacktracking
path of length at most eight is geodesic, uniquely so below eight.
A point sees a three-point line at point distances $(h,h+1,h+1)$, since
its incidence distance to the line is odd and at most seven.

There are real unit vectors $e_x$ with
$\langle e_x,e_y\rangle=\phi(d(x,y))$, where $\phi(j)=(-1/2)^j$.
Indeed the point graph has intersection array
$\{10,8,8,8;1,1,1,5\}$. Its distance matrices are polynomials in the
adjacency matrix $H$. Substitution in the three-term recurrence gives
$HF=-5F$ for $F_{xy}=\phi(d(x,y))$. As $F$ is a polynomial in $H$, it acts
by one scalar on the $-5$ eigenspace and vanishes on the other eigenspaces.
Its positive trace makes it positive semidefinite. The three vectors on a
line sum to zero by taking the squared norm. For $m$ points at mutual
distance three the squared norm of their sum is $m(9-m)/8$. In particular,
nine is the maximum permitted by this bound, and at nine the sum is zero.

For points $y,z$ at distance three, there is exactly one point $p$ at
distance two from both. In the incidence graph, the three paths between
$y,z,p$ have lengths $6,4,4$. Their union is a tree by girth. Its branch
vertex is the middle line of the unique $y,z$ geodesic, and $p$ is the
third point on that line. This also proves that at most two points of a
mutually distance-three set can be at distance two from a given $p$.
For any two, their geodesics from $p$ must take the same first line and
then different points of that line. Only two such points are available.

## The six points and their possible extensions

There are three geodesics from $A$ to $B$, one starting at each point of
$A$, obtained by its unique projection to $B$. Their interiors are disjoint
by girth. Write their line sequences as $A,U_j,C_j,V_j,B$ with $j=1,2,3$.
Along each path call the four successive intersection points $a_j,l_j,r_j,b_j$.
Write $u_j,t_j,v_j$ for the unused points of $U_j,C_j,V_j$, respectively.
Put $S=\{u_1,u_2,u_3,v_1,v_2,v_3\}$.

All points of $S$ are mutually at distance three, as are the three $t_j$.
The distances from $t_j$ to $u_i,v_i$ are two if $i=j$ and four otherwise.
For a direct way to check these values, distances among the twelve path
points are read along a shortest path in the three arms, of incidence
length at most eight. Replace each unused point vector by minus the sum of
its two neighbours on its line and take the inner product. For instance,
for $i\ne j$, the four distances between $(a_i,l_i)$ and $(r_j,b_j)$
are $3,4,4,3$. Their cosines sum to $-1/8$, so $d(u_i,v_j)=3$.
The same calculation gives $d(l_i,t_j)=d(r_i,t_j)=4$ for $i\ne j$.
The full 21-point calculation is in [theta_profiles.py](theta_profiles.py).

The line sums on the arms and on $A,B$ give
$\sum_{s\in S}e_s=\sum_j e_{t_j}$. A point $z$ at distance three from all
of $S$ therefore satisfies $\sum_j\phi(d(z,t_j))=-3/4$.
The only possible distance multisets are $(1,3,3)$ and $(1,1,2)$.
The latter contradicts the mutual distances of the $t_j$.
Suppose $z$ is adjacent to $t_i$. Each of the other two $t_j$ is at distance
three from both $z$ and $t_i$, hence at distance two from the third point
$w$ of their line. The middle-line fact above puts $w$ on the middle line
$M$ of the geodesic between those two $t_j$. All three $t_j$ now lie at
incidence distance three from $M$ with distinct projections. In particular
$M$ is uniquely determined by any two of them.

Conversely, suppose such a center $M$ exists. The line joining each $t_j$
to its projection on $M$ differs from $C_j$, because a geodesic from $t_j$
to $t_i$ cannot start through $l_j$ or $r_j$. The three row lines $A,B,M$
and three column lines $C_1,C_2,C_3$ have line distance two across families.
Each row and column uses three different attachment points. Lines in the
same family are opposite by the resulting nonbacktracking paths of length
eight. Let $D_{ij}$ be the unique intermediate line between row $i$ and
column $j$, and $x_{ij}$ its unused point. Apply the preceding six-point
calculation to each pair of rows. All nine $x_{ij}$ have mutual distance
three. The three unused points in the new row are exactly the possible
extensions of $S$, since every extension determines $M$, its projection
from $t_i$, and the third point on that joining line. This proves part 1.

## A global count for the first step

Form a graph on unordered opposite-line pairs. Join two vertices if each
of the two lines of one pair is at line distance two from both lines of the
other pair. A vertex $U=\{A,B\}$ has exactly three neighbours, the pairs
from its three middle lines $C_1,C_2,C_3$. Those middle lines are mutually
opposite by the length-eight paths through $A$. Thus the graph is cubic.

A four-cycle through $U$ is equivalent to the grid in part 1. Indeed two
neighbours of $U$ can be written $\{C_1,C_2\}$ and $\{C_1,C_3\}$. Each has
$A,B$ as two of its three middle lines. To share a further neighbour they
must have the same third middle line $M$. The paths from $C_i$ to $C_j$ via
$A$ and $B$ already use $l_i,r_i$, so the path via $M$ uses $t_i$ and $t_j$.
This gives the center above. Conversely the pairs of rows and pairs of columns
of a grid exhaust each other's neighbours, making an entire connected
component isomorphic to $K_{3,3}$. In particular failure of the first
condition propagates to all three neighbouring pairs.

For the adjacency matrix $J$ of this graph and its $N=2925\cdot2048/2$
vertices, the first condition somewhere is equivalent to
$\operatorname{tr}(J^4)>15N$. In any simple cubic graph the closed
four-walks that retrace an edge contribute 15 per vertex. In a $K_{3,3}$
component the total is 27 per vertex. All four-cycles occur in such
components, so $(\operatorname{tr}(J^4)-15N)/72$ counts completed grids.
Evaluating this trace in an arbitrary octagon is a six-line incidence
count. A model computation gives its value only for the Ree–Tits octagon.

There is also a simple obstruction to forcing a grid by projection counts
at one fixed $A$. Its 96 lines at line distance two split into 32 for each
of the three points of $A$. Each opposite partner selects a triple, one in
each group. Fix two coordinates. These lines are opposite through $A$, so
their other two middle lines are precisely the two possible partners of
$A$ with those coordinates. Both are opposite $A$ by the paths of length
eight. Thus each two-coordinate choice occurs twice. A repeated *whole*
triple is equivalent to the grid. To see that the pair margins permit zero
repeats, on three copies of $\mathbb Z/32$ take the triples with third
coordinate $x+y$ or $x+y+1$. This example concerns the array of counts
only.

## One centered transversal forces thirty lines

Call a triple of the nine tips *centered* if some line is at incidence
distance three from all three. Its projections are distinct, and each
pair's unique geodesic has this center line as its midpoint. At most three
tips can project in this way to a line, since tips with a common projection
are at distance at most two.

Fix a tip $x$ and a line $L=\{x,p,q\}$. Let $a,b$ count the other tips whose
unique geodesic from $x$ starts via $p,q$, respectively. At $p$ the other
tips have distances $2,3,4$ with multiplicities $a,b,8-a-b$, by projection
onto $L$. Taking the inner product of $e_p$ with the zero sum of all nine
tip vectors gives $(-8+4a-2b+8-a-b)/16=0$. Thus $a=b\le2$ by the
middle-line fact. A used line accounts for either two or four other tips.
Two other tips using the same first **point** $p$ give a centered triple
with $x$ by that same fact. The converse follows from the path to the center.

The grid line $D_{ij}$ through $x=x_{ij}$ accounts for four tips, its row
and column companions. The four others therefore use either one more line
with two through each of its other points, or two more lines with one through
each point. In the first case $x$ belongs to two further centered triples.
In the second it belongs to none. Each further centered triple must be a
transversal. Two tips in a common row or column already have that row or
column as the midpoint of their geodesic, with capacity for exactly its
three tips. There are precisely two transversals through $x$. If one is
centered, so is the other. Among the six permutations of three columns,
the graph of transversals meeting in a cell is $K_{3,3}$ by permutation
parity. This propagates one centered transversal to all six.

There are now twelve different centers (three rows, three columns and six
transversals) and eighteen lines through the nine tips, two per tip. The
latter are all different since tips have mutual distance three. No center
contains a tip, since it already has three other tips at incidence distance
three. At each tip its two lines supply four neighbours, the projections
to its four different centers. These four are distinct by the direction
counts, and neighbours of different tips are distinct by their distance.
The resulting 45 points lie on the 30 chosen full lines. Every tip belongs
to two chosen lines. Each neighbour belongs to its tip line and its center.
All line places have been accounted for, so the incidence graph has the
claimed degrees. It is connected through the initial grid and the lines
attached to its tips.

Suppress the degree-two points. This gives a connected simple cubic graph on
30 line vertices with girth at least eight. The two depth-three trees growing
away from any edge are disjoint and exhaust $2(1+2+4+8)=30$ vertices.
In particular every vertex is closer to one endpoint of each edge than to
the other. Distance parity from any vertex bipartitions the graph. The same
trees give diameter at most four, and a breadth-first count gives diameter
at least four and girth exactly eight. Subdividing every edge restores
incidence girth 16 and diameter eight. For instance a vertex is within
three of one endpoint of every edge, and opposite-parity endpoints of two
edges are at distance at most three. This is a full GO$(2,1)$.

## What the computations establish

[check_opposite_seed.py](check_opposite_seed.py) uses the stored ATLAS point
permutations, verifies transitivity on lines, checks all partners opposite
line 0, and checks the cubic intersection array of every distinct resulting
copy through that line. Any full thin copy containing the opposite pair must
contain the three forced geodesics, the center between any two $t_j$, and
the midlines of the transversal geodesics. Thus in the successful case the
copy through the specified pair is unique. Every thin copy has opposite
lines, which stay opposite in the ambient geometry by girth. The last count
below is a double count of line–copy incidences in this model.

```text
$ python thin/check_opposite_seed.py
line orbit 2925 opposite partners checked 2048
three extensions and six centered transversals 2048
verified thin copies through line 0 256 opposite partners per copy 8
unordered opposite pairs 2995200 thin copies 24960
```

[theta_profiles.py](theta_profiles.py) tests whether point-distance data can
force a seventh point. Each candidate distance vector obeys the triangle
inequalities on the seed, the projection rule on all eleven seed lines, and
uniqueness of the inward line through a seed point at distance one to three.
It excludes all vectors extending the six. The
[integer certificate](theta_no_extension.json) assigns multiplicities to the
remaining named profiles. Its compact form uses the symmetries of the
labelled diagram. The verifier expands it and checks **every individual**
ordered-pair intersection table, including repeated points, and every cubic
identity $\sum_z\phi(d(z,x))\phi(d(z,y))\phi(d(z,w))=0$ on the 21 points.
The cubic identities follow by taking the squared norm of
$\sum_z e_z^{\otimes3}$, which is
$1755\sum_{j=0}^4 k_j\phi(j)^3=0$ for shells $(1,10,80,640,1024)$.

```text
$ python thin/theta_profiles.py
admissible named profiles 756 orbits 79
occupied named profiles 393 points 1755
ordered pair tables 441, cubic entries 1771, exact residual 0
points at distance three from all six tips 0
```

This is an integer solution for point profiles only. Two profiles of positive
multiplicity still have to be realized together with consistent distances,
lines and geodesics between their points. Thus these marginal constraints
alone cannot force the seventh point. Neither existence of a grid in an
arbitrary octagon nor existence of a centered transversal in some such grid
has been proved here.

## How many direction patterns can avoid all six transversals?

Suppose a grid exists and all its transversals fail. At each tip $x$, the
four tips outside its row and column occupy two more lines through $x$.
On each line the paths to its two assigned tips take different neighbours
of $x$. A partition of four into two pairs has three choices, giving
$3^9=19683$ labelled assignments. Row and column permutations and transpose
act on the assignments with 374 orbits. These symmetries merely rename
the grid.

For an off-row, off-column pair $x,y$, write $v$ for the partner of $y$
at $x$, and $w$ for the partner of $x$ at $y$. The path from $x$ to $y$
has two inner points $a,b$. The profile of $a$ on the nine tips takes
values $1,2,3$ on $x,y,v$ and four elsewhere. The profile of $b$ is
obtained by exchanging $x,y$ and replacing $v$ by $w$. Apply the projection
rule to the third point $p$ of the line $ab$. If $v=w$, the three tips
$x,y,v$ would all be at distance two from $p$, which is impossible. For
$v\ne w$, the profile of $p$ is two on $x,y$, four on $v,w$ and three
elsewhere. In particular the eighteen points obtained from the eighteen
unordered off-row, off-column pairs are different and are adjacent to no tip.

The prescribed configuration now has 81 points and 51 lines. There are nine
tips, 18 neighbours on their original lines, 36 neighbours on their other
used lines, and the 18 third points just described. Its lines are the six
row and column centers, the 27 lines through tips, and the 18 lines $ab$.
These are distinct by their tip-distance profiles and incidences. The
additional ambient points and lines, including unused lines through the tips,
are unknown.

[grid_directions.py](grid_directions.py) computes necessary distance domains
on these 81 points. Paths of at most four point steps along the prescribed
lines are geodesic when taken shortest in this configuration. Triangle
inequalities against the tips further bound each distance. Repeatedly
projecting the three-point line rule and imposing symmetry eliminates 93
more classes. Two classes fail when the lines through each point are
considered jointly with uniqueness of the inward line at distances one
through three. One further class fails after branching on both possible
values of a single distance and applying the same rule. These are finite
exclusions using only integer and bit operations.

The class `001001220` (18 labelings) fails a further geodesic test. If
$d(r,b)=3$ and $d(r,c)=d(b,c)=2$, then **every neighbour $a$ of $r$**
satisfies $d(a,b)=d(a,c)+1$. Indeed $c$ is the third point on the middle
line of the geodesic from $r$ to $b$. Both geodesics from $r$ start with
the same line and the same next point $t$. For $a=t$ the two distances are
2 and 1. For the other neighbour on that first line they are 3 and 2.
For a neighbour on any other line, the paths via $r$ have lengths four
and three in the point graph and are geodesic by incidence girth.

In the zero-based numbering of `skeleton()`, the remaining distance from
28 to 37 is three or four. Impose each value and propagate the joint
distance patterns on the named lines through each point. The two branches
give these point distances.

| $d(28,37)$ | $(r,a,b,c)$ | $d(r,a)$ | $d(r,b)$ | $d(r,c)$ | $d(b,c)$ | $d(a,b)$ | $d(a,c)$ |
|---:|---|---:|---:|---:|---:|---:|---:|
| 3 | (35,40,49,78) | 1 | 3 | 2 | 2 | 3 | 3 |
| 4 | (36,44,50,76) | 1 | 3 | 2 | 2 | 3 | 3 |

Both rows contradict the neighbour identity. The verifier checks the
six displayed distances in both branches; it uses no SAT solver.

Each of the 25 classes surviving this test has a [saved distance matrix](grid_metrics.json).
The verifier checks the specified profiles, all lines and inward directions,
and every upper bound on the number of these 81 points in two given shells
about any two of them. The bounds are the intersection numbers in the full
octagon. In particular they enforce the triangle inequalities. For all four
nonprincipal point eigenspaces it checks that the cosine Gram matrix is
positive semidefinite by exact fraction-free elimination. The least
eigenspace Gram matrices have rank 30. A further
[geodesic check](grid_geodesics.py) introduces the unique path for every
pair of named incidence vertices at distance below eight. The paths of any
three named vertices with perimeter below 16 form a tree, hence their
overlapping segments are identified. For each saved matrix these
identifications preserve all 132 named vertices and introduce no path
shorter than their specified distances. Four old witnesses were replaced
when checking this condition. These certificates describe partial
metrics. They supply neither the other 1,674 points nor their incidences.

```text
$ python thin/grid_directions.py
directions 19683 symmetry classes 374
coincident partners 14179 classes 252
inconsistent line distances 4608 classes 93
inconsistent inward lines 72 classes 2
inconsistent after distance split 18 classes 1
inconsistent middle line 18 classes 1
verified partial metrics 788 classes 25
```

At this stage ruling out extensions of these 25 **direction types**, allowing all
distance matrices of each type, would suffice for the second step.
Eliminating only the displayed matrix of a type would leave its other
possible matrices open. The first step, existence of a four-cycle somewhere,
would still be required.

## Additional grids among the prescribed lines

Five surviving classes force other grids within their own 51 lines.
[grid_neighbors.py](grid_neighbors.py) enumerates the six center lines and
the nine intermediate lines in each such diagram. It requires 27 distinct
attachment points and tips. This gives length-eight nonbacktracking paths
between centers in the same family, so those centers stay opposite in the
ambient octagon. Only additional grids entirely exhibited in this diagram
are counted.

In each additional grid, four of the nine direction choices are already
fixed if no transversal is centered. A tip has two off-row, off-column
geodesics along the same initial line and through different next points.
The two geodesics have incidence length six inside the given diagram, so
they remain unique in every completion. They specify one pair of the
matching, hence also its other pair. Checking the 788 surviving labelled
words against the four fixed positions leaves 22 possibilities each time.

| Original class | Additional grids exhibited | Number of the 22 labelings in this class |
|---|---:|---:|
| `000000111` | 3 | 2 |
| `000011011` | 2 | 4 |
| `000011022` | 1 | 4 |
| `001002111` | 2 | 8 |
| `002002111` | 2 | 4 |

Thus in a hypothetical counterexample every grid obtained by these
transitions belongs to the same five-class list. A centered transversal
in any one of them would already give the thin suboctagon. This restriction
alone permits all five classes.

A [full-size incidence construction](rooted_grid.md) shows that even the
complete rooted path counts at all nine tips can coexist with zero centered
transversals. Its incidence graph has girth ten.

[Five-line neighbourhood refutations](grid_star_extension.md) reduce the 25 classes above to 16. [Edge labels and the 21-point support test](opposite_edge_labels.md) give further necessary conditions for a grid-free octagon.

[Middle-line identities in completed stars](grid_middle_extension.md) exclude four
more classes, leaving 12 classes (254 labeled words). The strengthened
[point and line marginal equations](theta_line_margins.md) have an exact
nonnegative rational witness without a seventh point. Neither step of the
uniqueness argument is yet complete.

[Coupled completed stars](grid_coupled_stars.md) exclude three further
classes, leaving nine classes (164 labelled words). [Incidence directions,
star partitions and integrality congruences](theta_star_margins.md) sharpen
the no-extension marginal system and provide a rational witness satisfying
all new star equations, integral point orbit totals and the derived parity
conditions. Both steps of the uniqueness argument remain open.

[Tree patterns in coupled stars](grid_tree_stars.md) exclude `000000222`,
leaving eight classes (158 labelled words). The same note shows that the
numbers of centered, 16-cycle and 18-cycle point triples at mutual distance
three are fixed by the parameters, and why this alone does not force a grid.
