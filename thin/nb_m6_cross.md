# Six-step endpoint returns in the opposite-pair graph

In a grid-free generalized octagon of order $(2,4)$, every
non-backtracking six-step walk of opposite pairs whose endpoint pairs
share a line retains that line at **each even step**. No line of the
initial pair can be at line distance one from a line of the final pair.
If all pairs in the component have the 18-cycle shape, distance two is
also impossible and seven-step walks have no common endpoint line.

## Computing the sixth walk distribution

Use $D,R,S,J,Q_k$ as in the [boundary identities](opposite_boundary.md)
and [non-backtracking moment calculation](nb_moment_sdp.md).
At $U=\{A,B\}$ let $C_i$ be its three middle lines and $V_i$ the neighbor
formed by the two $C_j$ with $j\ne i$. Write $Z_i$ for the third middle of
$V_i$ beyond $A,B$. The other neighbors of $V_i$ are $\{A,Z_i\}$ and
$\{B,Z_i\}$. Write $Y_{iA},Y_{iB}$ for their third middles beyond the two
ends of $V_i$. Then $S(V_i)$ consists of $C_i,Y_{iA},Y_{iB}$. If $T$ is
the opposite-pair-by-line matrix recording the six $Y$'s, we have
\[
  JS=R^{\mathsf T}+T,\qquad
  J^3D^{\mathsf T}=14R^{\mathsf T}+2T .
\]
Indeed $J^2D^{\mathsf T}=6D^{\mathsf T}+2S$ and $JD^{\mathsf T}=2R^{\mathsf T}$.

Extend the tip geodesics on each of the 32 labelled failed metrics by
the six geodesics locating the $Y$ lines, as in the fifth-moment verifier.
For each pair of lines retain all 39-pattern candidates agreeing with
the resulting point domains and read their line distances. In every case
$Y_{iA}$ has distance two to $A$, four to $B$, and four to every $C_j$
with $j\ne i$. It has distance two or three to $C_i$. Exchange $A,B$
for the other side. Two different $Y$'s on one side are opposite.
The distance $d(Y_{iA},Y_{iB})$ is two or three. For the six pairs
$(Y_{iA},Y_{jB})$ with $i\ne j$ the possible distances are contained in
$\{2,3,4\}$. The following are counts of *available* distance-two positions.

| Shape of $U$ | Labelled metrics | Available positions among the six |
|---|---:|---:|
| 18 | 8 | 0 |
| 16 | 12 | 2 |
| 16 | 12 | 3 |

Let $u,v$ be the means over $U$ of the numbers of those six pairs
*actually* at distance two and three. If $p$ is the fraction of
16-shape pairs, then $0\le u\le3p$, $v\ge0$ and $u+v\le6$.
[Distance enumeration and exact checks](nb_m6_cross.py)

Among the three pairs within $S(V)$, precisely one has distance two for
a 16-shape $V$, and none does for an 18-shape $V$. As $U$ runs over the
three neighbors of $V$, the distinguished $C_i$ runs over all three
entries of $S(V)$. The close pair hence contributes twice to the
$C_i$–$Y_{i*}$ count and once to the $Y_{iA}$–$Y_{iB}$ count. Averaged
per $U$, these counts are $2p$ and $p$.

Take the Gram matrix of $14R^{\mathsf T}+2T$. The following are the mean
row counts entering it, before multiplication by 1,024 opposite pairs
per ambient line. The middle row counts one orientation.

| Lines | Weight | Distance 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|---:|
| $C,C$ | 196 | 3 | 0 | 0 | 0 | 6 |
| $C,Y$ | 56 | 0 | 0 | $2p$ | $6-2p$ | 12 |
| $Y,Y$ | 4 | 6 | 0 | $2p+2u$ | $6-2p+2v$ | $24-2(u+v)$ |

Subtracting the shorter walk terms $87Q_0+47Q_2+11Q_4$ gives
\[
 Q_6=1024\bigl(24,\ 0,\ 32p+8u,\ 96-32p+8v,\ 264-8(u+v)\bigr).
\]
A walk that keeps one designated endpoint at every even step has three
choices of first edge and then two choices at steps three and five.
Each even step is uniquely determined by the endpoint and the preceding
step. This gives $2\cdot3\cdot2^2=24$ marked returns per start pair,
exactly accounting for $Q_6(0)$. The arguments apply within each cubic
component by restricting all pair sums to it and rescaling by its size.

## Seventh and eighth steps

The diagonal contributions in $(14R^{\mathsf T}+2T)^{\mathsf T}J
(14R^{\mathsf T}+2T)$ come only from $4T^{\mathsf T}JT$. Indeed for
adjacent $U,V$, a middle of $U$ either belongs to $V$ or is opposite
both endpoints of $V$. It cannot be a $Y$ of $V$. It also cannot be a
middle of $V$, by the cross-distance calculation for adjacent triples.
The lower odd moments have zero diagonal. If a $Y$ of $U$ occurs in $T(V)$
it must be at distance two from an endpoint of $V$. Among the three
$C$'s only its own $C_i$ can qualify, and this $C_i$ belongs to two
neighbors of $U$. There are exactly two qualifying $Y$'s per 16-shape
metric and zero per 18-shape metric. Thus $Q_7(0)\le16384p$.

An endpoint retained at the even steps of an eight-step walk gives
$2\cdot3\cdot2^3=48$ marked walks per start. It is unique along the
entire walk since the pair at step two shares only that endpoint with
the start. Its two cross comparisons with the other end at the start
and finish are oppositions. An indexed walk on an eight-cycle gives one
additional coincident endpoint at step eight. With $t$ the fraction of
indexed incidences on longer cycles this gives
\[
 Q_8(0)\ge1024(60-12t), \qquad Q_8(4)\ge96\cdot1024 .
\]
For a long indexed walk the partner eight steps later is a different
leaf of the same branch, at distance $3-k$ where $k=0,1,2$ is the shared
prefix length. At adjacent starts $k+k'\le2$ by the
[matching prefix bound](opposite_leaf_maps.md). Thus
$Q_8(1)+Q_8(2)+Q_8(3)\ge12288t$ and $Q_8(2)+Q_8(3)\ge6144t$.
At step seven three endpoint comparisons on every indexed walk have distance
two. The fourth compares leaves on opposite sides. The leaf at step minus
one is the unique match in the branch containing the step-seven leaf, so
the fourth distance is two on an eight-cycle and three or four on a long
cycle. This gives $Q_7(2)\ge36864+12288(1-t)$ and
$Q_7(3)+Q_7(4)\ge12288t$. [Branch blocks](indexed_cycle_blocks.md)

## What the projected moment tests still allow

At $p=0$, $t=1/2$, take $s=2$ in the fifth-moment formula and $(u,v)=(0,3)$
in the sixth. Complete the distributions with
\[
 Q_7=(0,0,250244,458591,77597),\qquad Q_8=(55296,0,0,12288,1505280).
\]
These numbers satisfy the bounds above. For each nonconstant line
eigenvalue $7,3,-1,-3$, rational symmetric elimination gives strictly
positive pivots for the order-eight Hankel matrix and all three
localizing matrices for $3+\lambda$, $3-\lambda$ and $9-\lambda^2$.
The constant eigenvalue gives the point mass at three. Thus the
strengthened moment tests still allow an all-18 distribution. This
certificate establishes feasibility of these projected inequalities only.
Run `python thin/nb_m6_cross.py` from the repository root for the geometry
and rational certificate.
