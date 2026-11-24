# Point and line margins on an opposite pair with no grid

The 21-point profile equations in [opposite_seed.md](opposite_seed.md)
permit an integer solution without a seventh point, but the five-line
support test rules out that particular solution. One can strengthen the
linear equations by counting **both points and lines** around all 32 named
incidence vertices (21 points and 11 lines). The resulting system still
has a nonnegative rational solution with no seventh point. This is an
obstruction to using just these linear equations and nonnegativity; it is
not an integral point-line configuration.

## Line profiles and distance tables

Start with the no-extension profiles on the 21 points given by
`theta_profiles.py`. Repeatedly discard profiles that cannot participate
in a five-line star, using the exact support test of
`theta_incident_lines.py`. In this case the process deletes just three
profiles, leaving 753 named point profiles in 78 orbits. A possible line
type is a multiset of three such profiles whose values on each of the 21
coordinates are a permutation of $(h,h+1,h+1)$. There are 8,982 named
types in 824 orbits.

For either kind of observer, the incidence distance from a named point to
a prospective line (or a prospective point to a named line) is twice the
minimum point distance plus one. For two different lines the incidence
distance is one plus the minimum of the distances from the three points
of one line to the other. A prospective line with the three named points
of a named line must equal that line and has distance zero instead.
The profiles therefore give incidence distances to all 32 named vertices.

To calculate the necessary joint counts, fix a vertex $u$ of either type.
The distance partition about $u$ is equitable: a vertex at distance $i$
has one neighbour at distance $i-1$ for $1\le i\le7$, all its neighbours
at distance seven for $i=8$, and all remaining neighbours at distance
$i+1$. The degree $d_i$ at level $i$ is five or three according to whether
its vertices are points or lines. Let $Q$ be this $9\times9$ radial
adjacency quotient and $D=\operatorname{diag}(d_i)$. The matrices of
nonbacktracking path counts on the distance partition are

$$W_0=I,\quad W_1=Q,\quad W_2=Q^2-D,\qquad
W_j=QW_{j-1}-(D-I)W_{j-2}\quad (3\le j\le8).$$

For the recurrence, prepend one edge to a path and subtract the paths
whose first new edge immediately reverses the next edge. For $j=2$ there
are $d_i$ such choices; for $j\ge3$ there are $d_i-1$ choices excluding the
first edge of the shorter path. Girth sixteen means a nonbacktracking
path of length at most eight is geodesic and is unique below eight.
Between vertices at distance eight all $d_i$ first edges yield a geodesic:
the neighbour on that edge is at odd distance at most seven from the other
vertex and cannot be closer than seven. Thus divide row $i$ of $W_8$ by
$d_i$ to get a distance-relation quotient $R_8$; put $R_j=W_j$ for $j<8$.
For a second named vertex $v$ at incidence distance $k$ from $u$, the
entry $(R_j)_{k,h}$ counts all vertices at distances $h$ from $u$ and $j$
from $v$. Their type is determined by the relevant distance parity.

The verifier cross-checks the point-point-point entries against the
independently calculated point-graph intersection numbers, and the
shell sizes for each base type against the breadth-first counts.

## Linear equations and certificate

Use a nonnegative multiplicity for each named point profile and each
named line type. A completion without a seventh point must satisfy:

- For each pair of the 32 named vertices, including repetitions, the
  joint distance tables for both types of observer just derived.
- The cubic least-eigenspace cosine identities on all triples of the 21
  point coordinates, including repetitions.
- For each point profile, the number of line places carrying that
  profile is five times its multiplicity. For each coordinate at distance
  one through three, exactly one line place per point goes inward in that
  coordinate. Line places are counted with multiplicity when a line has
  two points of the same profile.

The arm permutations and exchange of endpoints preserve the equations.
For rational feasibility there is no loss in averaging over these twelve
relabellings. The certificate assigns a common rational multiplicity to
each **named** profile or line type in an orbit. This averaging would
not justify requiring integral common multiplicities: an integral
completion could have asymmetric counts. The verifier generates each
named equation, merging only identical coefficient rows with the same
right-hand side. There are 1,531 distinct rows in 902 columns.

The [certificate](theta_line_margins.json) lists the nonzero numerators
in orbit order and a common denominator (26 digits). It was found by
solving the rational relaxation numerically and reconstructing the basic
solution with modular arithmetic. [Verification](theta_line_margins.py)
does not use the numerical solver or rational reconstruction: it uses
Python integers to multiply every row by the numerators and compares to
the integer right-hand side times the common denominator.

```text
$ python thin/theta_line_margins.py
profiles 753 78 line types 8982 824
matrix (1531, 902) 147766
exact rational witness 406 occupied orbits
all equations satisfied, total points 1755, total lines 2925
```

This strengthens the known marginal obstruction but does not provide an
integer solution to the strengthened equations, let alone a generalized
octagon. The inward counts are aggregated over each point profile, not
a decomposition of its occupied line places into actual five-line
stars. Those individual star partitions, integrality, compatibility
of distances between unnamed vertices and the remaining girth constraints
are possible further restrictions. The existence of a grid somewhere in
an arbitrary generalized octagon of order $(2,4)$ remains open.
