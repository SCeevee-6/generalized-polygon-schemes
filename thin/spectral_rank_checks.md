# Spectral and rank routes that do not close the gap

These are checks of four other ways to force a grid. None of them does,
and the numbers below say how far each falls short. The last section
records an exact consequence of the homology of the incidence graph that
has not been used yet.

## The least eigenspace is only a 3-design

Normalized vectors of the eigenvalue $-5$ (multiplicity 78) give 1,755
unit vectors in $\mathbb R^{78}$ with inner products $(-1/2)^{d}$. The
Gegenbauer sums $\sum_i k_i\,C^{(38)}_j((-1/2)^i)$ over the distance
shells vanish for $j=1,2,3$ and equal $4101435/2$ at $j=4$. So in every
GO(2,4) the point set is a spherical 3-design and not a 4-design. It has
four inner products other than one, so the Delsarte–Goethals–Seidel
rigidity statements, which need strength at least $2\cdot4-2$, do not
apply.

## Positive semidefiniteness on the twelve-cycle witnesses

For every eigenvalue $\theta$ of the point graph, the cosine matrix
$[\omega_\theta(d(x,y))]$ on any set of points is positive semidefinite
of rank at most the multiplicity. For each of the 13 joint
[twelve-cycle witnesses](indexed_cycle_paths.md), the verifier takes the
propagated distance domains, merges slots at distance zero and searches
for large sets of slots with all pairwise distances determined (45 to 65
of the 207 to 295 points). All four nontrivial cosine matrices are
positive semidefinite there, and the rank at $-5$ is 23 to 33, far below
78. Determined parts of these diagrams are too small for the rank bound
to matter.

## Rank of the local matrix over GF(2)

For the [local rank criterion](opposite_cycles.md) note that lines at
line distance 2, 3 and 4 have 1, 0 and 0 common concurrent lines. Hence
$M_{A,C}=N_{S_1}^{\mathsf T}A_1N_{S_2}$ exactly, where $N$ is the
point–line incidence matrix, $A_1$ point collinearity, and $S_1,S_2$ the
two sides of $J_{A,C}$. Since the block only involves line distances
2, 3 and 4, the matrices in the line-distance algebra that agree with
$M$ on it modulo 2 are exactly $L_2+\beta I+\gamma L_1$. In the
Ree–Tits octagon their 2-ranks are 1594, 2925, 2924 and 869. For points,
$\operatorname{rank}_2(A_1+I)=1025$ and $\operatorname{rank}_2 N=1675$,
while the opposition matrix has 2-rank only 26. A bound
$\operatorname{rank}M<48$ cannot come from any of these global matrices.

The Ree–Tits octagon sits at the other extreme. At a line $A$, all 96
indexed systems $J_{A,C}$ consist of 32 four-cycles, so $M_{A,C}$ has
rank 32 and every opposite pair lies in a grid.

## An identity from the cycle space

Let $\Gamma$ be the incidence graph. Its binary cycle space $Z$ has
dimension $8775-4680+1=4096$ in every GO(2,4), and apartments span it.
An edge $UV$ of the cubic graph $J$ determines one apartment, the union
of the two geodesics between the lines of $U$ through the lines of $V$.
Each apartment contains four opposite pairs $P_0,\dots,P_3$ in cyclic
order and arises from exactly the two edges $P_0P_2$ and $P_1P_3$. Write
$\sigma$ for this pairing of edges. At a vertex $U$, the three geodesics
between its lines give the relation that its three apartments sum to
zero in $Z$.

The map from $\mathbb F_2^{E(J)}$ to $Z$ is therefore onto and vanishes
on the star space $B$ (of dimension $N-c$ for $c$ components) and on the
$\sigma$-invariant vectors $S$ (of dimension $3N/4$). Comparing
dimensions gives
\[
 4096\ \le\ \tfrac{3N}{2}-(N-c)-\tfrac{3N}{4}+\dim(B\cap S),
\]
that is,
\[
 \dim K\ \ge\ 4096+\tfrac N4,\qquad
 K=\{x\in\mathbb F_2^{V(J)}:\ x(P_0)+x(P_1)+x(P_2)+x(P_3)=0
 \text{ for every apartment}\}.
\]
Equivalently, the $N\times\frac{3N}{4}$ incidence matrix between opposite
line pairs and apartments has 2-rank at most $\frac{3N}{4}-4096$. This
holds in every GO(2,4), with or without grids. Turning it into a grid
criterion would need an independent computation of $\dim K$.

## Status

The uniqueness problem is still open in the literature. The known
characterizations need extra hypotheses: an automorphism group regular
on the points opposite a fixed point and fixing its lines (Cohen,
O'Brien and Shpectorov, J. Algebra 421 (2015) 369–393), or a full
suboctagon of order $(2,1)$ (De Bruyn).

Run `python thin/spectral_rank_checks.py` from the repository root
(about a minute; it rebuilds the Ree–Tits octagon from
`data/generators.npy`).
