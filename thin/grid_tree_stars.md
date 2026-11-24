# Tree patterns in coupled stars

The coupled-star test leaves nine direction classes. Imposing, at every
neighbour of a completed star, the pattern forced by a small geodesic tree
excludes one more, `000000222`. This leaves eight classes (158 labelled
words). The exclusion is of the whole class. No saved distance matrix is used.

## The tree patterns

Let $r$ be a point and $b,c$ points with $i=d(r,b)$, $j=d(r,c)$, $k=d(b,c)$,
where $1\le i,j\le3$ and $i+j+k\le7$. The three geodesics have total
incidence length at most 14, so their union is a tree. If $i+j+k$ is even,
it branches at a point. The arm from $r$ then has point length
$\alpha=(i+j-k)/2$. If it is odd, it branches at a line and
$\alpha=(i+j-k-1)/2$ is the point distance from $r$ to its projection.

For a neighbour $a$ of $r$ put $\Delta=(d(a,b)-i,\ d(a,c)-j)$. If $a$ is on
a line through $r$ outside the tree, both paths through $r$ have incidence
length at most eight. They are nonbacktracking and hence geodesic, giving
$\Delta=(1,1)$. The other neighbours lie on the tree lines at $r$. There
are three cases.

| Case | Tree at $r$ | Possible $\Delta$ |
|---|---|---|
| $\alpha\ge1$ | both geodesics start with the same line and point | $(-1,-1)$, $(0,0)$, $(1,1)$ |
| $\alpha=0$, $i+j+k$ even | different first lines | $(-1,1)$, $(0,1)$, $(1,-1)$, $(1,0)$, $(1,1)$ |
| $\alpha=0$, $i+j+k$ odd | one first line, different next points | $(-1,0)$, $(0,-1)$, $(1,1)$ |

In the first row the second point of the shared line is at distance $i$ and
$j$ from $b$ and $c$. In the last row $r$ lies on the branch line. The
middle-line identity of [grid_middle_extension.md](grid_middle_extension.md)
is the first row with $(i,j,k)=(3,2,2)$.

Every completion of the configuration assigns values to all the variables
below, so each clause is a necessary condition.

## The formula

Start from the [coupled-star formula](grid_coupled_stars.md) for roots
$r,s$. For each root, take as targets the 80 other named points and the
unnamed neighbours of the other root. For every pair of targets $b,c$ and
every value of $(d(r,b),d(r,c),d(b,c))$ in the table above, add clauses
forbidding the other values of $\Delta$ at each of the ten neighbours of
$r$. A target may coincide with a neighbour; distance zero is then allowed
and is consistent with the patterns. The distances between unnamed
neighbours of different roots are the existing coupled variables.

The tree clauses against the other star's unnamed neighbours are essential.
With roots $29,73$ the formula with tree clauses on named targets only is
satisfiable. It is also satisfiable with no tree clauses at all.

## Search and result

For each of the nine classes, root pairs at named distance one to three
were tested in order of increasing distance: all of them for eight classes
and the first 407 of 1,662 for `000000222`, where the scan stopped after
three refutations, for 12,384 formulas in all. Only class `000000222`
failed, at the pairs $(29,73)$, $(30,69)$, $(33,74)$. Every other formula was
satisfiable; none reached the 120-second limit. The
[certificate](grid_tree_refutations.json) for $(29,73)$ is a Glucose proof
with deletions dropped, replayed by [grid_rup.py](grid_rup.py) without a
SAT solver.

```text
$ python thin/grid_tree_stars.py
000000222 excluded at roots 29 73 RUP additions 418
remaining classes 8 labelings 158
000000000 000000022 000000111 000011022 001001222 001002111 001012102 002002002
```

Class `000000222` is not one of the classes in the additional-grid table,
so the restriction there is unchanged: `000000111`, `000011022`,
`001002111`. A single root with the tree clauses refutes none of the nine
classes; all 81 roots were tried for each.

[Verifier](grid_tree_stars.py)

## Triples at mutual distance three

This part concerns the first step, forcing a grid. Let $x,y,z$ be points at
mutual distance three. As in [opposite_edge_labels.md](opposite_edge_labels.md),
the three geodesics form either a tree with a center line ("centered"), a
16-cycle or an 18-cycle. Call the triple's shape at $x$ (a), (b) or (c) if
the geodesics from $x$ to $y,z$ start with different lines, the same line
and different points, or the same line and point. Shape (c) is the centered
case, at all three points. A 16-cycle triple has shape (b) at exactly one
point and (a) at the others. An 18-cycle triple has shape (a) everywhere.

**The number of triples of each type is the same in every generalized
octagon of order $(2,4)$.**

- Centered: a line $L$ and three tips at incidence distance three from it,
  through different points of $L$, give $(4\cdot2)^3$ triples. The paths
  through $L$ have length six and are geodesic. The center is unique, so
  there are $2925\cdot512$.
- 16-cycle: fix $x$, a line $H=\{x,p,q\}$, and require $y$ to be reached
  through $p$ and $z$ through $q$. Then $y$ is at incidence distance five
  from $H$ through $p$ and $z$ similarly through $q$. If $d(y,z)=3$,
  the path $y\ldots p\,H\,q\ldots z$ of length ten and the geodesic $y z$
  close an apartment. Its element opposite $H$ is the line $w$ in the
  middle of the $y z$ geodesic. Conversely, a line $w$ opposite $H$ gives
  $y,z$ at distance five from $H$ on the geodesics through $p$ and $q$.
  The pairs are therefore in bijection with the 2048 lines opposite $H$.
  Counting the shape (b) point gives $1755\cdot5\cdot2048$ triples.
- 18-cycle: the remainder of $1755\cdot640\cdot p^3_{33}/6$, with
  $p^3_{33}=200$ from the intersection array.

The counts are $1{,}497{,}600$, $17{,}971{,}200$ and $17{,}971{,}200$.

In the notation of [opposite_seed.md](opposite_seed.md), an opposite line
pair $U$ lies in a grid exactly when its unused points $t_1,t_2,t_3$ form a
centered triple. A grid would follow if the number of pairs $U$ over a
given triple depended only on the triple's type, but **it does not**. In the
Ree–Tits octagon every opposite pair lies in a grid. Half of the centered
triples come from four pairs and the rest come from none.
[distance3_triples.py](distance3_triples.py) checks the three totals in the
model, together with this 0/4 split over all centered triples through one
point. This rules out that counting route, and neither step of the
uniqueness argument is closed.

```text
$ python thin/distance3_triples.py
p33 200 triples 37440000
centered 1497600 predicted 1497600
16-cycle 17971200 predicted 17971200
18-cycle 17971200 predicted 17971200
opposite pairs over centered triples through point 0 {0: 1280, 4: 1280}
```
