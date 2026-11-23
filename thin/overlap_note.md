# Reading O(2) off its thin suboctagons

O(2) has 24,960 suboctagons of order (2,1); they form one orbit of 2F4(2). Join two of them when they share 13
points. The 13 shared points are always the radius-two ball of a point p inside both copies, and p is the centre of
the edge. The resulting graph is 45-regular with two components of 12,480 vertices, swapped by the outer
automorphism. Shell sizes from a vertex are 1, 45, 1170, 9300, 1964, and the intersection sizes by graph distance
are: distance 1: 13; distance 2: 2 or 5; distance 3: 0, 1, 3 or 9; distance 4: 0 or 5.

The automorphism group of one component, as an abstract graph, has order 17,971,200 = |T| (BLISS via igraph). So
every automorphism of the bare graph comes from the octagon and the centre labelling is invariant.

`overlap_recognition.py` recovers the labelling from the uncoloured component alone:

1. In a K4, opposite edges have the same centre.
2. Take a chordless 4-cycle B, B1, C1, C in which both diagonals have exactly 4 common neighbours. Opposite edges
   have the same centre exactly when the vertices outside the cycle that are adjacent to exactly two cycle vertices
   span 8 edges (otherwise they span 4).
3. Rules 1 and 2 never join edges with different centres. On all 280,800 edges they produce 3,510 classes of 80
   edges, two per centre. Two classes belong to the same centre exactly when no edge of the graph joins their vertex
   sets (for other pairs there are at least 45 such edges).
4. Two centres are collinear exactly when they label two edges of a common triangle. This gives all 8,775 collinear
   pairs, so the whole octagon.

Rules 1–4 were found by looking at the labelled graph and then checked against every label, so this is a
recognition procedure for the known component, not yet a characterisation of abstract graphs with these local
properties. The 4-cycle counts per edge are 2340 in total: 720 chordless with diagonal counts (2,2), 360 + 360 with
one chord, 360 chordless rectangles and 360 chordless non-rectangles with counts (4,4), and 180 in K4s.
