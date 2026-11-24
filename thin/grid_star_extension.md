# Extending the five lines at a grid point

Of the 25 direction classes left by the 81-point tests, nine cannot occur
in a generalized octagon of order $(2,4)$. This leaves 16 classes, or 362
of the original 19,683 labelled words. Two further exclusions apply to
particular distance assignments within surviving classes.

| Excluded direction class | Roots whose neighbourhoods give a contradiction |
|---|---|
| `000011011` | 39 |
| `000012102` | 77, 36, 37 |
| `001010100` | 29 |
| `001010102` | 31 |
| `001010122` | 31 |
| `001010200` | 29 |
| `001010202` | 72, 30 |
| `001020120` | 31 |
| `001021120` | 31 |

## Why a completion would satisfy the formula

Use the point numbering of `skeleton()`. For every pair of its 81
points introduce a single distance variable shared by both orders. Start
with the domains from the line and inward-line propagation in
[grid_directions.py](grid_directions.py). For each named observer and
each named line, require that the three distances are a permutation of
$(h,h+1,h+1)$.

At a root $r$ put its ten neighbours in five pairs, one for each line.
Fill the first pairs with the neighbours on its named lines. The other
pairs have fresh distance variables to every named point. A completion
can order the remaining lines and their two neighbours arbitrarily. For
each named target $y$ require the following.

- The distances from $y$ to $r$ and the two members of each pair obey the
  three-point line rule.
- If $0<d(r,y)<4$, exactly one of the ten neighbours is a step closer to
  $y$. The geodesic of incidence length at most six is unique.
- Each new neighbour sees a valid distance triple on every named line.
  Its distances to any two named targets obey the triangle inequalities
  with the variable distance between the targets.

Several roots share the same named distance variables. Distances between
two fresh vertices at different roots are unnecessary for these
contradictions. The construction allows fresh vertices from different
roots to denote the same ambient point. It uses no assumption about
automorphisms of a completion. Every completion assigns values to the
formula, so a refutation excludes it.

For the nine rows of the table, all named distances range over their
domains. For each of the next two words add just the displayed
conjunction. Each conjunction is impossible and occurs in the corresponding
saved matrix. Alternative assignments remain open.

| Class | Root | Forbidden simultaneous distances |
|---|---:|---|
| `000000111` | 27 | $d(27,33)=4$, $d(27,37)=3$, $d(27,41)=3$, $d(33,56)=3$ |
| `001002111` | 37 | $d(27,33)=4$, $d(27,37)=4$, $d(29,37)=4$, $d(37,60)=3$ |

Of the five classes in the earlier additional-grid table, `000011011`
is excluded. Hence a further grid exhibited by a surviving class in that
table has at most 18 of the original 22 labelled assignments: the four
assignments of class `000011011` are impossible.

## Reproduction

Run `python thin/grid_star_extension.py` with NumPy, SciPy and
`python-sat` installed. For every contradiction the verifier builds a
Boolean formula, obtains a Glucose 4 refutation, and replays its clause
additions by reverse unit propagation with MiniSat 2.2, through the empty
clause. It then counts the excluded orbits by row and column permutations
and transposition. The nine whole-class proofs use no entries from
`grid_metrics.json`. The original 81-point witnesses remain available to
reproduce the weaker tests.

[Verifier](grid_star_extension.py)

[The next star test](grid_middle_extension.md) uses the middle-line identity
at all ten neighbours and excludes four further whole direction classes.
