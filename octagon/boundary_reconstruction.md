# Boundary reconstruction and joint outer-shell constraints

These results apply to any generalized octagon of order (2,4), without a
regularity assumption on its Clebsch covers. They do not prove uniqueness.
They strengthen the finite tests in [boundary_span.md](boundary_span.md).

## Reconstruction from the boundary metric

Fix a root x, write S = Gamma_3(x), O = Gamma_4(x), and use point-graph
distances. Form the graph T on S with edges at ambient distance two.
Its 5-cliques are canonically the sets of five boundary neighbours of
the 1024 opposite points.

Indeed, every opposite point has five such neighbours on distinct lines.
For three pairwise distance-two points, the union of their unique
length-four incidence geodesics is a tree: any cycle in the union would
have length at most 12, contrary to girth 16. The equal-length paths
therefore form a tripod centered at a point, two incidence steps from
each endpoint. For any further member of a clique the same argument
with the first two points forces the same center (the unique midpoint
of their geodesic). A center at distance two from x has only four
outward lines, a center at distance three has just one neighbour in S,
and a closer center has none. Thus every 5-clique has an opposite
center, and each center supplies exactly one clique.

Distinct cliques Q_v and Q_w intersect in at most one point by uniqueness
of the length-four geodesic between two boundary neighbours. Adjacent
opposite points share the boundary point on their common line. If two
cliques intersect at u, choose a in Q_v minus u and b in Q_w minus u.
When v and w are adjacent, the path from a to b through their common
line has incidence length six. Otherwise the path through u has
incidence length eight. In both cases the path is simple (the respective
boundary neighbours occupy different lines), and no shorter path exists:
it would give a cycle of length less than 16. Thus v and w are
adjacent exactly when d(a,b)=3; otherwise d(a,b)=4. Disjoint cliques
are not adjacent.

The boundary metric therefore determines the complete uncoloured outer
graph, including its number of components. Every boundary point occurs
in eight cliques; these eight outer points induce four disjoint edges.

## Exact Laplacian identities

Use H for the 640 by 1024 incidence matrix of these cliques, A for the
outer adjacency matrix, L=5I-A, and B for the boundary Gram matrix of
the previous note. Put F=(g_{d(v,u)}) with v in O and u in S, and
G_O=(g_{d(v,w)}) on O; g=(64,32,8,0,-1).

The ambient Gram matrix satisfies adjacency eigenvalue 5. Taking its
outer rows gives

    L F = H^T B,
    L G_O = H^T F^T,
    L G_O L = H^T B H.

Also 1^T F=0: the distance-four matrix is a polynomial in ambient
adjacency, so its product with the full Gram matrix is a scalar
multiple of that Gram matrix. The root-to-boundary entries vanish.

If the shell is connected (equivalently rank B=350), the unique
mean-zero solution is F=L^+ H^T B. In particular candidate boundary
metrics do not require independent profile choices. Entries of F
must lie in {32,8,0,-1}. For every prime p and z satisfying Lz=0 mod p,
integrality requires BHz=0 mod p. This follows by multiplying LF by
z^T and transposing; it is a rejection test, not a sufficient condition.

The earlier Schur residual identity gives G_O=F B^+ F^T+J/64 in the
connected case. Substitution gives

    G_O = L^+ H^T B H L^+ + J/64.

Hence M=64 L^+ H^T B H L^+ must have diagonal 4095 and off-diagonal
entries in {2047,511,-1,-65}. These are exact rational arithmetic tests.

## Joint spectral and girth constraints

Write C for ambient point adjacency. Its spectrum is
10,5,1,-3,-5. For z on O orthogonal to constants, extend z by zero.
Then ||Cz||^2=||Az||^2+||Hz||^2 <=25||z||^2. Thus

    H^T H <= 25I-A^2   on the orthogonal complement of 1.

Write H^T H=5I+A+R. Distinct opposite points share at most one boundary
neighbour. Adjacent points share exactly one, and every other pair
sharing a boundary point belongs to different outward lines there.
Any outer path of length at most five joining that pair would give
a cycle of incidence length at most 14. (The path through the boundary
has interior disjoint from any shortest outer path, including its
line vertices unless a shorter path to an endpoint already exists.)
Equivalently the unique ambient geodesic of point length two uses
the boundary, and an outer path of length at most five cannot form
a second incidence path without violating girth. Hence R is a symmetric
zero-diagonal 0-1 matrix supported at outer distance >=6, including
infinite distance, and R1=30*1.

For j=0,...,5, trace(R A^j)=0, so

    trace(H^T H A^j)=5 trace(A^j)+trace(A^{j+1}).

The outer graph has girth at least eight, so its first seven nontrivial
walk moments are those of the 5-regular tree. For j=0,...,7 these
traces are 1024,0,5120,0,46080,0,496640,0. The corresponding mixed
moments at j=0,...,5 are 5120,5120,25600,46080,230400,496640.

Thus the following simultaneous constraints hold:

* 5I+A+R is positive semidefinite and has rank at most 640;
* 20I-A-A^2-R is positive semidefinite on 1-perp;
* R is a 30-regular 0-1 graph supported on outer distances >=6;
* A is a 5-regular cover of Clebsch with girth at least eight;
* H has row sum 8, column sum 5, with four disjoint outer edges per row.

These constraints are necessary for all component patterns. The open
connected case additionally requires a single 64-sheeted cover.
No infeasibility certificate is currently known. Showing these
constraints impossible in the connected case would force a disconnected
root, but still would not by itself recognize a component or prove
unrestricted uniqueness.

## Reproduction

Run python octagon/boundary_reconstruction.py from the repository root.
The checker reconstructs both the outer point set and its edges from
boundary distances, compares to the true Ree--Tits metric, and checks
the exact Laplacian identities and mixed moments. It is a positive
control in the disconnected model, not an exhaustive connected-case
search.
