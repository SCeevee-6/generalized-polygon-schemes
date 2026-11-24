# Involutions and failed opposite pairs

Start with the 21 points in theta_profiles.py. Adjoin all three middle-tip
geodesics with extend_paths, allowing repeated point slots, and propagate
the line constraints of grid_path_closure.py. Branch over every value of
one undecided distance and repeat. This yields 33 labelled matrices. The
twelve permutations of the arms and the endpoints act on all new slots
because a short geodesic and the third point of a line are unique.

The edge labels in opposite_edge_labels.md distinguish one centered matrix,
two 18-cycle orbits of sizes 6 and 2, and three 16-cycle orbits of sizes
12, 6 and 6. In each failed orbit every stabilizer element swaps the
endpoints exactly when its arm permutation is odd. Run
`python thin/failed_pair_types.py` to reproduce the enumeration.
The checker also compares every two-anchor shell count of distinct named
points, and dually of distinct named lines, against the ambient intersection
number. All five orbits pass these upper bounds. Simultaneous ambient
multiplicities remain to be determined.

## Fixed flags and exchanged opposite pairs

Let an involution fix P points, F lines and I incident point-line flags.
Write S=(5P-I)/2 for the total number of transpositions of incident lines
at fixed points, and T=(3F-I)/2 for the corresponding number of point
transpositions on fixed lines. Let B count unordered exchanged opposite
line pairs for which the involution preserves all three arms individually.

At incidence distances 2, 4 and 6 the geodesic between exchanged lines
has a unique midpoint. Conversely choose a transposition at a fixed
midpoint and extend on one specified side. The image supplies the other
half. The respective numbers of unordered pairs are S, 4T and 8S.
At incidence distance 8 this same construction yields 32T marked pairs:
on one side of the middle line there are 4 choices of next line, 2 of
next point and 4 of final line. Every resulting path is geodesic by girth.
The marking selects an invariant arm. An involution on the three arms fixes
one or three, so the number of unordered exchanged opposite pairs is
32T-2B. Counting all 2925 lines gives

$$4B=F+18S+72T-2925=109F+45(P-I)-2925.$$

Every pair counted by B violates the stabilizer parity of a failed pair.
It therefore yields a grid. In particular a grid-free octagon must satisfy
the exact fixed-flag equality with B=0 and hence 45 divides F. The converse
divisibility implication need not hold. No involution is assumed to exist
in an arbitrary octagon.

For a point-centered central collineation the radius-four fixed ball gives
P=91, F=45 and I=135. The equality yields B=0 for an involution of this
type, so this count alone gives no grid. A line-centered fixed ball is
already excluded by Proposition 4.5 applied to the dual parameters (4,2).
See [Parkinson, Temmermans and Van Maldeghem](https://cage.ugent.be/~hvm/artikels/225.pdf).
