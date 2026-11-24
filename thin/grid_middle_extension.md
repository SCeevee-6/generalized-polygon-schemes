# The middle-line identity on all ten neighbours

The five-line neighbourhood test leaves 16 direction classes. Four more of
these are impossible, leaving 12 classes (254 of the original 19,683 labelled
words). The exclusion applies to every metric of each of the four classes.

## Necessary condition

Use the notation and numbering of [opposite_seed.md](opposite_seed.md) and
the five-line neighbourhood formula of [grid_star_extension.md](grid_star_extension.md).
Suppose the named points $r,b,c$ satisfy

$$d(r,b)=3,\qquad d(r,c)=d(b,c)=2.$$

The incidence paths between them have lengths six, four and four. Their
union is a tree by girth sixteen. Thus $c$ is the third point on the middle
line of the $r,b$ geodesic. The paths from $r$ towards $b$ and $c$ use the
same first line and the same next point $t$.

For every neighbour $a$ of $r$ the distances to $(b,c)$ are $(2,1)$ if
$a=t$, $(3,2)$ for the other point of the first line, and $(4,3)$ if $a$
is on another line through $r$. In the last case the paths via $r$ are
nonbacktracking of lengths eight and six in the incidence graph, hence
geodesic. In particular $d(a,b)=d(a,c)+1$.

The earlier 81-point test used this identity on named neighbours. It also
applies to each **unnamed** neighbour in a completed five-line star. Build
the previous star formula with unknown distances to all 81 named points.
For each named $b,c$ for which the three antecedent distances are possible,
and each of the ten rows for neighbours $a$ of $r$, forbid every choice of
$d(a,b),d(a,c)$ violating the identity when the antecedent holds. This is a
clause with the three negated antecedent literals and two negated distance
literals. The old star formula and these new clauses are necessary in every
generalized octagon containing the prescribed grid. No entries of the saved
partial metric for the direction word are fixed.

## Exact exclusions

| Direction class | Root $r$ | Labelled words | RUP additions |
|---|---:|---:|---:|
| `000000002` | 36 | 18 | 58 |
| `000001220` | 59 | 36 | 32 |
| `001010222` | 29 | 36 | 65 |
| `002002111` | 52 | 18 | 77 |

Each formula has a short [clause certificate](grid_middle_refutations.json).
The independent [checker](grid_rup.py) computes unit propagation directly:
for every proposed clause it temporarily negates all its literals and
requires a contradiction using the initial formula and previous additions.
It then adds the clause and finally requires a contradiction without any
assumptions. No SAT solver is needed to replay these proofs. The clauses
were obtained from a SAT search, but the verification is by integer
operations alone. The original star refutations continue to use python-sat.

```text
$ python thin/grid_middle_extension.py
000000002 excluded at root 36 RUP additions 58
000001220 excluded at root 59 RUP additions 32
001010222 excluded at root 29 RUP additions 65
002002111 excluded at root 52 RUP additions 77
remaining classes 12 labelings 254
```

The remaining representatives are

```text
000000000  000000022  000000111  000000222
000002002  000011022  001001222  001002111
001012102  001110221  002002002  002012102
```

The additional-grid calculation in [opposite_seed.md](opposite_seed.md)
originally allowed five classes and 22 labelings for each additional grid.
Excluding `000011011` and now `002002111` leaves three of those classes:
`000000111`, `000011022`, `001002111`, with respectively 2, 4 and 8
labelings (14 in total). This does not exclude any of the three classes.

It remains to eliminate the 12 direction classes or otherwise force a
centered transversal in at least one grid. The existence of a grid in an
arbitrary generalized octagon remains a separate open step.

[Coupling two completed stars](grid_coupled_stars.md) excludes three of
these twelve classes, leaving nine.
