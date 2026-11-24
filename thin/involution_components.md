# Involutions on opposite components

Let $p$ be the center of an involution in a grid-free generalized octagon
of order $(2,4)$. For a point $y$ opposite $p$, let $w(y)$ count the
geodesics from $p$ to $y$ whose midpoint is fixed. The seven possible
branch distributions and their global weight counts were obtained in
[the parity argument](involution_branch_parity.md).

## Each component has the same mean

Let $C$ be one collinearity component on the points opposite $p$ and
write $n=|C|$. Each point at incidence distance six from $p$ has
exactly $n/128$ neighbors in $C$, by the
[boundary count](../octagon/README.md).
A line $H$ at distance three from $p$ has sixteen descendant points
at distance six. Each $y$ in $C$ whose path uses $H$ uses exactly one
such descendant, by uniqueness of the path in that branch. Thus every
one of the forty lines $H$ occurs for exactly $n/8$ points of $C$.

If $Q$ of the forty lines are selected, summing over them gives
\[
\sum_{y\in C}w(y)=nQ/8.
\]
In the 31-point case $w$ takes only the values one and five. In the
71-point case it takes only three and five. The mean therefore gives
the following exact counts for every component, whether or not it is
preserved by the involution.

| Fixed points | Weight one | Weight three | Weight five |
|---:|---:|---:|---:|
| 31 | $15n/16$ | 0 | $n/16$ |
| 71 | 0 | $5n/8$ | $3n/8$ |
| 91 | 0 | 0 | $n$ |

## Four-cycles in a quotient component

Suppose the involution preserves $C$ and form the quotient graph by
its pairs of vertices in $C$. Since each point of $C$ is opposite
its image, their distance in the outer graph is at least four.
Consequently the quotient is simple, and projection is a graph
cover with fibers of size two. The outer graph has girth at least
eight. Indeed, each ambient line contains at most two opposite
points, so an outer cycle gives an incidence cycle of twice its
length. A quotient cycle of length below eight therefore cannot
lift to a closed cycle of that length. It must lift to a path from
a point to its image. In particular the quotient has girth at least
four.

A quotient four-cycle through the orbit of $y$ lifts in its two
directions to two distinct four-edge outer paths from $y$ to its
image. Each path is an ambient geodesic. The involution interchanges
the two paths after reversal. The geodesics counted by $w(y)$ have
fixed midpoints in the inner shell, so none lies in the outer graph.
The other $5-w(y)$ geodesics therefore give at most $(5-w(y))/2$
quotient four-cycles through the orbit. In particular:

- In the 31-point case the $n/32$ orbits of weight five lie on no
  four-cycle. Each other orbit lies on at most two.
- In the 71-point case the four-cycles are vertex-disjoint and there
  are at most $\lfloor5n/64\rfloor$ of them.
- In the 91-point case the quotient has girth at least five.

These are necessary tests for a component and its involution.

## A check on the binary parity condition

The matrix with rows indexed by the 1024 points opposite $p$ and
columns indexed by the 40 lines at distance three records whether a
geodesic uses the line. A selector vector must map to the all-ones
vector over $\mathbb F_2$. The matrix and its restrictions to the two
opposite components all have binary rank 29 in the Ree–Tits model.
Enumerating the resulting 2048 parity solutions gives the following
counts for the total selected counts admitted by the fixed-tree theorem.

| Branch distribution | Parity solutions in the model |
|---|---:|
| $(2,2,2,2,2)$ | 16 |
| $(0,4,4,4,8)$ | 240 |
| $(4,4,4,4,4)$ | 288 |
| $(6,6,6,6,6)$ | 48 |
| $(8,8,8,8,8)$ | 1 |

The computation uses only the parity equations on the model, so the
counts certify feasible selector vectors for those equations. Whether
an involution realizes one of those vectors is a separate condition.
Run [the verifier](involution_parity_model.py) with
`python thin/involution_parity_model.py` from the repository root.

The component constraints leave the four fixed-point counts open.
