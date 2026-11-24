# Coupling two completed five-line stars

The [middle-line star test](grid_middle_extension.md) leaves 12 direction
classes. Completing two stars at once and imposing their mutual distance
relations excludes three more, leaving nine classes (164 labelled words).
The exclusions are of entire direction classes, not of particular saved
81-point distance matrices.

## Necessary distances between neighbours

Start from the usual 81-point metric variables and complete the stars at
two named points $r,s$ using the five-line formula and middle-line identity.
For each star, retain its named neighbours and introduce slots for the
unnamed ones. Distinct slots in the same star must have point distance one
if they are on the same line through the root and two otherwise. These
distances follow from the nonbacktracking incidence paths of lengths two
and four, respectively.

Between the two stars, introduce a single symmetric distance variable for
each pair of unnamed slots. A slot at one root may represent the same
ambient point as a slot at the other root, so distance zero is allowed.
Their distances to any named point must satisfy the triangle inequalities.

There is a further restriction when the roots are not opposite. Let $a$
be a neighbour of $r$, $b$ a neighbour of $s$, and put

$$d=d(r,s),\qquad p=d(a,s),\qquad q=d(b,r),\qquad h=d(a,b).$$

If $d=2$ or $3$, the unique geodesic from $r$ to $s$ determines whether a
neighbour is on its first point, on the first line at the other point,
or on another line. These three cases give $p=d-1,d,d+1$ (and similarly
for $q$). The resulting path between $a,b$ has point length
$\ell=p+q-d$, obtained by deleting, replacing or extending the ends of the
root geodesic. In the case $d=2$ and $p=q=1$ the two inward neighbours
are the same point. Otherwise the resulting path is nonbacktracking in
the incidence graph; its length is at most ten. Therefore

$$h=\ell\quad(\ell\le4),\qquad h\ge8-\ell\quad(\ell>4).$$

Indeed, if a path of incidence length $2\ell>8$ had a shortcut of length
$2h$ with $2h+2\ell<16$, their union would contain a cycle of length less
than sixteen.

For $d=1$, both roots lie on one line. If $p=0$ or $q=0$, one neighbour
is the other root and $h=1$. If $p=q=1$, both neighbours are the third
point of the root line and $h=0$. In the remaining cases $p,q\in\{1,2\}$
and $h=p+q-1$, by the short paths on or through the root line. No extra
condition of this form is imposed for $d=4$.

Finally, each unnamed neighbour of either root is also an observer of
the other root's five full lines. Its distances to the three points of
each line must have the usual projection pattern $(j,j+1,j+1)$. If its
distance to the other root is one through three, exactly one neighbour
of that root is a step closer. A completion of the ambient octagon gives
all the variables in this formula; no condition assumes that unnamed
slots at different roots are distinct.

## Refutations

| Class | Roots | Labelled words | Saved RUP additions |
|---|---:|---:|---:|
| `000002002` | 32, 37 | 18 | 99 |
| `001110221` | 33, 36 | 36 | 775 |
| `002012102` | 39, 52 | 36 | 38 |

The [certificates](grid_coupled_refutations.json) were obtained by SAT
search and trimmed to clause additions verified by reverse unit
propagation. The replay uses only the original Boolean formula and the
saved clauses, and requires no SAT solver.

```text
$ python thin/grid_coupled_stars.py
000002002 excluded at roots 32 37 RUP additions 99
001110221 excluded at roots 33 36 RUP additions 775
002012102 excluded at roots 39 52 RUP additions 38
remaining classes 9 labelings 164
```

The remaining representatives are

```text
000000000  000000022  000000111  000000222
000011022  001001222  001002111  001012102
002002002
```

None of the newly excluded classes were among the three remaining
additional-grid types, so the previous additional-grid restriction is
unchanged. The nine classes must still be addressed, and a separate
argument is needed to force a grid in an arbitrary octagon.

[Verifier](grid_coupled_stars.py)

[Tree patterns](grid_tree_stars.md) exclude `000000222`, leaving eight classes.
