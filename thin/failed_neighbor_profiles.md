# Simultaneous point profiles at four neighboring opposite pairs

At a failed pair, take the [36-slot middle-tip metric](failed_pair_types.md).
It has 34 or 36 distinct points. It contains the 21-point seeds of the
original pair and each of its three neighbors in the cubic graph. For
a neighbor, start with its two endpoint lines, match their points at
distance three, and read the internal path points and the third point
of each path line. The verifier checks the entire 21-by-21 distance
matrix after each such reconstruction.

The following two necessary tests survive in **all five** failed
metric orbits. They were designed to seek a contradiction at a
single opposite pair together with all three of its neighbors.

For each of the four seeds start from the 753 point profiles surviving
the [five-line test](opposite_edge_labels.md). Join the four tables by
exact agreement on all shared points. A profile specifies distances to
the named lines as well, by taking the minimum point distance, doubling,
and adding one. Impose all triangle inequalities involving the observer
and the short four-point test with three named incidence vertices.
For the latter, if the two smallest of the three cross sums add to
less than sixteen, the two largest must be equal.

1. Assign multiplicities modulo two to the surviving point profiles.
   Impose the complete point-observer distance margins at every two named
   points or lines, using the [incidence quotients](theta_line_margins.md).
   Also impose the cubic least-eigenspace identities on all triples of
   named points. The scaled cosine has odd value only at distance four,
   so the cubic row is the intersection of three opposite shells.
   Binary elimination with the right-hand side attached checks consistency.
2. Give each profile in the closed point neighborhood of the named
   points an integer multiplicity. A profile adjacent to two named points
   has multiplicity at most one by girth. For every named point choose
   five incident lines with integer multiplicities. The profiles on
   each candidate line must obey the incidence-distance rule at every
   named point and line. At each named point the selected lines use
   exactly the selected neighbor multiplicities, and provide exactly one
   inward line to each named target at incidence distance one through
   seven. When two named points have point distance three, match their
   selected internal path neighbors. A match requires consecutive profiles
   compatible with collinearity and a surviving profile for the third
   point on their line.

The first check is over the binary field and the second uses ordinary
nonnegative integers. The verifier solves the second system and checks
all rows and bounds again in integer arithmetic. The two checks are
separate. Neither assigns incidence throughout the octagon.

Run `python thin/failed_neighbor_profiles.py` from the repository root.

The forced-line computation remains feasible in all five orbits. This
checks a failed pair and all three neighbors jointly, including all lines
named in that 36-slot metric. It does **not** check all lines of an entire
twelve-cycle diagram at once: distances between slots from different
36-slot metrics can remain undetermined. In particular, it does not
classify the 13 saved twelve-cycle witnesses under the proposed full
two-seed diagram test, and it does not rule out other assignments.
