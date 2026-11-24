# A nine-tip grid forces a full thin suboctagon

**Conditional theorem (computer-assisted).** Let $\mathcal O$ be a
generalized octagon of order $(2,4)$. If an opposite-line pair in
$\mathcal O$ admits the nine-tip grid of [opposite_seed.md](opposite_seed.md),
the grid has a centered transversal. The thirty-line construction there
gives a full suboctagon of order $(2,1)$. De Bruyn's Main Theorem then
identifies $\mathcal O$ as the Ree–Tits octagon. The remaining existence
problem on this route is to force one such grid in an arbitrary octagon.
[De Bruyn (2015)](https://doi.org/10.1016/j.disc.2015.04.022)

Assume a grid has no centered transversal. At each of its nine tips the
four off-row, off-column geodesics are paired by their initial lines. There
are three matchings per tip. The seed argument in
[opposite_seed.md](opposite_seed.md) rules out coincident partners and,
otherwise, prescribes 51 full lines and the distances of 81 distinct
points to the nine tips. The verifier imports the enumeration and this
seed construction from [grid_directions.py](grid_directions.py). It
initializes every other off-diagonal distance among the 81 points with all
four positive values. It uses neither the stored partial metrics nor the
earlier star refutations.

## The 39 ordered line-pair patterns

For two ordered full lines, record the nine point distances in a square.
The incidence graph has girth 16 and diameter eight. Every nonbacktracking
path of length at most eight is geodesic, uniquely so below eight.

- Coincident lines give the six reorderings of the square with zeros on
  the diagonal and ones elsewhere.
- Lines at line-graph distance $h+1$, for $h=0,1,2$, have a unique
  connecting geodesic. Its attachment points $a,b$ give the square with
  entry $h+[i\ne a]+[j\ne b]$ in position $(i,j)$. The longest path used
  has incidence length $2h+4\le8$. This gives 27 squares.
- Opposite lines give six squares with threes on a matching and fours
  elsewhere. Each point on the first line projects uniquely to the second
  at incidence distance seven from the line. Reversing the geodesic
  between the two lines pairs their attachment points bijectively.

All 39 possibilities remain available for two newly named lines until
distance information restricts them.

## A test on three full lines

Let $a,b,c$ denote candidate squares on $AB,BC,AC$, respectively. Together
they give all incidence distances among twelve labels: three lines and
their nine point positions. Within each square, the distance from a point
to the other line is one plus twice its row or column minimum. The
line-line distance is zero for coincident lines, and otherwise two plus
twice the minimum of the square.

Discard a candidate triple if it violates a triangle inequality. For
every four labels form the three sums of distances obtained by pairing
them in the three possible ways. Write the sums in order $s_1\le s_2\le s_3$.
A second necessary condition is

$$s_1+s_2<16\quad\Longrightarrow\quad s_2=s_3.$$

Indeed the shortest of the three cyclic tours of these four labels has
length $s_1+s_2$. Realize each leg by a geodesic. A closed walk of length
less than 16 has tree support, since its support contains no cycle. It
traverses every edge of that tree at least twice. The tree therefore has
at most seven edges, and all its path distances agree with the ambient
distances. In a tree the two largest pairing sums agree: each edge
separating one label from three contributes equally to all sums, and
all two-against-two separating edges induce the same split. This argument
also permits repeated labels.

The checker enumerates the $39^3$ candidate triples once. A pattern on
$AC$ must have support through **each** middle line $B$: some patterns
still allowed on $AB,BC$ must form a retained triple with it. The checker
alternates these restrictions with projections between line-pair squares
and individual point-distance domains. A real configuration would retain
its actual square and distance in every domain at every iteration.

## Completing a forced path

Each certificate consists of batches of endpoint pairs. Before a batch,
the verifier requires each requested pair to have a singleton distance
$d=2$ or $3$. It inserts its unique geodesic and a full three-point line
on every edge. A previously named point forced to distances $k,d-k$ from
the endpoints lies at position $k$ on that geodesic and can be reused.
Otherwise the position gets a fresh label. The checker sets the distances
between positions along the path and repeats the three-line test.

A fresh label starts with distances zero through four to the other
points, subject to its defining incidences. Different fresh labels may
therefore denote the same ambient point. Two paths may share a point or
line even when the verifier gives them separate labels. An empty domain
rules out every completion of the seed.

## Replication

The independent replay of the last eight classes uses at most 105 point
labels and 66 full lines. Counts include the original 81 points and 51
lines; labels need not remain distinct outside that seed.

| Direction word | Forced paths | Point labels | Full lines |
|---|---:|---:|---:|
| `000000000` | 1 | 86 | 54 |
| `000000022` | 2 | 87 | 55 |
| `000000111` | 1 | 86 | 54 |
| `000011022` | 1 | 86 | 54 |
| `001001222` | 4 | 99 | 62 |
| `001002111` | 4 | 97 | 61 |
| `001012102` | 3 | 96 | 60 |
| `002002002` | 6 | 105 | 66 |

The same replay covers the earlier cases, including `000000222` with
two distance-two paths. Orbit counts and sizes are recomputed from the
row and column permutations and transpose. These merely rename the
diagram; no symmetry of an ambient octagon is assumed.

| Contradiction | Direction classes | Labelled words |
|---|---:|---:|
| Coincident partners | 252 | 14,179 |
| Three-line test on the seed | 96 | 4,698 |
| Three-line test after adding forced paths | 26 | 806 |
| **Total** | **374** | **19,683** |

Run from the repository root with Python, NumPy and SciPy.

```sh
python thin/grid_path_closure.py
python thin/grid_path_closure.py 002002002
python thin/grid_path_closure.py --control
```

The positive control builds the Ree–Tits model from the stored ATLAS
generators. With one line fixed it chooses a second line at each of the
five distances, orders both lines in all 36 ways, and takes all 2,925
choices of a third line. The 526,500 triples cover all 39 first-pair
patterns and retain the true third square. A separate path test uses a
reordered copy of a line and two forced geodesics sharing an unnamed
interior point. All true distances survive, including zero between the
two fresh labels of that point. The model supplies this control only.
The exclusions use the universal rules proved above.

[Checker](grid_path_closure.py) · [Path certificates](grid_path_certificates.json)
