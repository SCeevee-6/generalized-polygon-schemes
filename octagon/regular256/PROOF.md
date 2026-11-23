# A regular 16-sheeted trace component is impossible in GO(2,4)

Computer-assisted research proof for independent verification. Here regular refers to a deck group acting transitively on each fibre. The underlying graphs of all the covers under consideration have valency five.

## An unconditional spectral bound

Let A be point adjacency and O the opposite layer at a point x. The intersection array is {10,8,8,8;1,1,1,5}, and the eigenvalues are 10,5,1,−3,−5. The least eigenspace has dimension 78. Its normalized point vectors e_y have Gram matrix G with entries (−1/2)^d(y,z), a scalar multiple of the orthogonal spectral projector. Their total frame operator is (1755/78)I. On each three-point line the vectors sum to zero. This representation is used by Cohen, O'Brien and Shpectorov in Proposition 7.4. [Paper](https://www.math.auckland.ac.nz/~obrien/research/octagon.pdf)

Write S_i for the sum of e_y e_y^T over points y in shell i about x. On the −5 eigenspace the distance matrices have eigenvalues (1,−5,20,−80,64), by the three-term recurrence. Their sum weighted by 4^(−i) is zero on this eigenspace. This weighted matrix is G entrywise squared. Thus T=sum_y G_xy e_y e_y^T vanishes: its squared Frobenius norm is the quadratic form of the entrywise square at column x of G. It follows that

    S_0 − S_1/2 + S_2/4 − S_3/8 + S_4/16 = 0.

Let P_0 project on e_x, P_1 on the five normalized sibling differences at depth one and P_2 on the forty at depth two. These 46 vectors are orthonormal. For points on distinct lines at the same depth the four cross distances of their child points agree, since the inward tree realizes paths of point length at most four. At different depths a point on the shallower line is equidistant from the two children on the deeper line. This also gives orthogonality to the root. Write P_3 for the complementary 32-dimensional projector. Summing line by line using the line sum gives S_0=P_0, S_1=(5/2)P_0+(3/2)P_1 and S_2=2S_1+(3/2)P_2. Using the total frame and the identity above,

    S_4 = 15I − 6S_0 + 2S_1 − 2S_2 = 4P_0 + 12P_1 + 12P_2 + 15P_3.

In particular a unit vector supported on O has squared norm at most 15/(1755/78)=2/3 after projection onto the −5 eigenspace. This follows from the equality of the norms of a restriction map and its transpose. Every other eigenvalue is at least −3. The adjacency Rayleigh quotient on O is at least −13/3. Thus no outer component can be bipartite, since a bipartite 5-regular component has eigenvalue −5. The proof uses no group assumption.

## Spectrum for a hypothetical 256-point component

For a component C put f=1_C−(|C|/1024)1_O. For a proper component its Rayleigh quotient is 5, the largest eigenvalue on the complement of constants. Hence Af=5f. On any of the 640 boundary points in shell three, the two outer endpoints of each of the four outward lines belong to the same component. Evaluation of Af=0 there shows exactly |C|/256 outward lines are assigned to C. The same formula holds directly when O is connected. Thus for |C|=256 each boundary point corresponds to one of the 640 outer edges.

Let H be the adjacency of C, N the unsigned edge-by-vertex incidence (rows indexed by boundary points), and M the involution swapping boundary siblings on their inward line. Then N^T N=5I+H. For a unit H-eigenvector g of eigenvalue l, extended by zero, Ag equals lg on C and Ng on the boundary, with zero elsewhere. Its spectral measure for A has mean l and variance 5+l. The polynomials (t+3)(t−1) and (t+5)(t+3) are nonnegative on the five-point ambient spectrum. Integrating them gives (l+1)(l+2)>=0 and (l+5)(l+4)>=0. The preceding lower bound removes l=−5.

The ambient distance-two and distance-three matrices restricted to C are respectively H²−5I and H³−9H+N^T M N. A two-step geodesic between points of C cannot leave C, since both endpoints on the boundary point's unique C-line are already adjacent. A three-step geodesic leaving C must go through two adjacent siblings in the boundary. Indeed a geodesic with one outer step and a boundary shortcut would contain the preceding impossible two-step route. A route using a different outer component in its interior cannot both leave and return in three steps. Conversely the sibling route gives a length-six nonbacktracking incidence path, hence a geodesic by girth 16. Outer nonbacktracking walks of three steps are likewise unique geodesics, giving H³−9H. Expanding the Gram matrix now gives

    16 G restricted to C = J + 3(−H³+H²+6H−N^T M N).

For l other than 5 use g orthogonal to constants and |g^T N^T M N g| <= 5+l. Positivity of the Gram form implies −(l+1)(l²−2l−5)>=0. Together these give the necessary intervals [−4,−2] and [−1,1+sqrt(6)] for every nonprincipal eigenvalue. Eigenvalues of any connected quotient cover lift to eigenvalues upstairs. All tests below use the relaxed upper bound 69/20, which exceeds 1+sqrt(6).

## Exhaustive cover calculation

Use the even five-bit words with edges at Hamming distance four for the Clebsch base. The general trace construction follows by projecting an opposite point to the five root lines. Adjacent outer points change four bits, and each outer vertex has one edge of each colour. Thus a component of size 256 gives a connected 16-sheeted cover. For a regular such cover the fibre group G has order 16. There is a normal subgroup chain giving a regular quotient of order four followed by arbitrary double-cover lifts at degrees eight and sixteen. For example take a central series with successive factors of order two, which exists by induction through the nontrivial centre of a finite 2-group. The order-four quotient is C4 or C2 x C2. Our last two lifts range over all double covers regardless of regularity over the original base.

A second test uses the 170 vertices in the two depth-three trees branching away from each end of a fixed edge. These vertices are distinct at the top by girth at least eight. If a quotient has n vertices, the number of visits to any one projected vertex is at most m=256/n. For a possible double lift of the quotient, these visits must split into two groups of at most m/2. The scripts apply the test to every edge.

**First double quotient.** Switching along a spanning tree leaves 25 free chord signs. We enumerate all 2^25 signings modulo base translations and coordinate permutations (19,132 orbits). A nonzero switching class cannot have signed eigenvalue 5: equality at the spectral maximum on a connected base forces all signs to equal products of vertex signs. Compute the signed characteristic polynomial by integer Newton sums, count its real roots across the interval endpoints by Descartes after a rational shift, and keep 245 orbit records. The zero eigenvalue of a shifted polynomial is handled separately. Signed 64-bit arithmetic is sufficient for Newton sums since 16 times 2^16 times 5^16 is below 2^63. The polynomial shifts use arbitrary precision.

If H_s denotes the signed adjacency, the signed visit counts for the edge tree are the row (e_u + sign(uv)e_v)^T (H_s³−8H_s). To verify the formula, the rows for the four successive levels are obtained from the starting row by I, H_s−I, H_s²−H_s−4I and H_s³−H_s²−8H_s+4I. These count the outward nonbacktracking paths, measured relative to the sheet at u. Thus the two loads are (unsigned count ± signed count)/2. The capacity 8 per vertex leaves 48 orbit types, comprising 27,979 signings.

**Elementary order-four quotient.** Enumerate sign spaces {0,a,b,a+b} whose three nonzero vectors pass the first test. Move a vector of least orbit index to the stored representative a. Enumerate every b with b<a+b and orbit indices of b,a+b at least that of a. The four loads are given by Fourier inversion using the three signed counts and the unsigned count. Requiring each load at most four leaves 16 records.

**Cyclic order-four quotient.** For each first signing a in the 48 representative classes, orient base edges lexicographically and take their Z/4 voltage equal to a+2b. Tree voltages are zero and b ranges over all 25-bit words. The 32-vertex first quotient has sheets j=0,1. Relative to binary coordinates t=j+2k, crossing an oriented edge changes the upper bit by b+(j AND a) over F2. The binary switch is the same on reversing that lifted edge. Thus every projected tree visit has an affine bit label in b. If there are eight visits at a quotient vertex then exactly four must receive each label. The XOR of those eight labels must vanish. Solve those equations, enumerate every affine solution, and check the full capacity inequalities also at vertices with five to seven visits. Exactly 152 cyclic records remain. The faithful complex character of C4 gives a Hermitian block with entries 1,i,−1,−i. Integer real/imaginary matrix arithmetic and the forbidden open interval (−2,−1) exclude 136, leaving 16 records.

**Arbitrary final double covers.** For each remaining four-sheeted quotient set edge signs to zero on a spanning tree (equivalently keep variables only on chords). Each visit in an edge tree gives the parity of a walk's chord mask, relative to the sheet at its first endpoint. Group by destination vertex. Normalize each group's labels by XOR with the first, which preserves the capacity condition up to swapping the two sheets. Whenever exactly m visits reach a vertex the XOR of their labels is m/2 modulo two. Solve these linear equations, enumerate the affine solutions, then impose the two separate capacity bounds. Repeat on the next quotient. The elementary order-four branch gives 40 degree-eight records and 40 degree-sixteen records. An exact two-colouring test finds that all 40 final graphs are bipartite. The cyclic branch gives 16 degree-eight records and zero degree-sixteen records. These counts allow repetitions across quotient descriptions. Consequently a regular 256-point component is impossible.

## Scope and remaining gap

This excludes all deck groups of order 16. It also excludes any 256-point component, regular or otherwise, factoring through one of the regular order-four quotients followed by two double lifts. A nonregular degree-sixteen action need not admit such a chain. In particular the argument does not yet remove 256 from the general list of possible component sizes. The known Ree–Tits model has two 512-point components. Recognizing one of them suffices by the previous half-reconstruction and De Bruyn's theorem on thin suboctagons. [De Bruyn](https://doi.org/10.1016/j.disc.2015.04.022)

The 512-point branch is treated separately in `../README.md`. An [independent implementation](README.md) now reproduces the exclusion and bypasses the cyclic spectral filter in the terminal check.

## Extension to every tower of four double covers

**Claim.** A 256-vertex opposite component cannot factor as four successive
connected double covers of its trace base. In particular, the monodromy
of a hypothetical 256 component must have order divisible by an odd prime.

Take one of the 48 possible first double quotients. For its next double
cover, fix zero voltage on a spanning tree and vary all 80−32+1=49 binary
chord voltages. There is no condition of regularity over the 16-vertex
base. Group the two height-three edge trees by destination in the 32-vertex
graph. Each walk determines the parity of its chord mask. At the 64-vertex
stage each bin splits into two of size at most four. Forbid all five-element
monochromatic subsets in each bin and enumerate every satisfying assignment.
This gives 632 four-sheet records, with possible repetitions up to isomorphism.

For such a signing, the signed adjacency M on the 32 vertices acts on the
antisymmetric subspace of the double lift. Its eigenvalues occur in the final
graph. Search numerically for a unit vector near an eigenvalue in one of the
forbidden intervals (−∞,−4), (−2,−1), (69/20,5). Round the centre to t/1000
and the vector to an integer vector w. Let d be the distance from integer t
to the complement of the selected interval scaled by 1000. Reject only if
the following **integer** inequality holds:

    ||1000 M w − t w||² < d² ||w||².

The spectral theorem puts an eigenvalue strictly inside the interval. Failed
rounding simply keeps a candidate. This is an exclusion certificate even with
a different numerical eigensolver. The remaining four-sheet records are
expanded through the same exhaustive capacity-based double lifts as above.
A two-colouring of every final graph contradicts the universal lower bound
−13/3.

    four-sheet capacity records                       632
    retained after certified spectral exclusions      66
    eight-sheet capacity records                     106
    sixteen-sheet capacity records                    64
    bipartite sixteen-sheet records                   64

To relate this calculation to monodromy, let G be the permutation group on
the 16 sheets induced by closed walks and H the stabilizer of one sheet.
If G is a 2-group, choose successive maximal overgroups from H up to G
with each smaller subgroup maximal in the next. Every such index is two.
The corresponding equivariant maps of coset spaces give four connected
double covers. Conversely the permutations preserving a tower of binary
fibres lie in the automorphism group of a depth-four rooted binary tree,
which is a 2-group. This covers nonregular examples as well as regular
ones. Transitive degree sixteen alone does not imply 2-group monodromy.

Run `OPENBLAS_NUM_THREADS=1 python octagon/regular256/double_towers.py`.
The [first-quotient breakdown](double_towers_output.txt) records every
count. The table counts switching classes above a fixed representative at
each earlier stage, rather than isomorphism classes of final graphs.
