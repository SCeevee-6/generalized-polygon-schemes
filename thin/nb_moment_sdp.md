# Non-backtracking moments and a semidefinite test

This extends the [projected moments](opposite_boundary.md) beyond $m_4$.
The setup is exact. With the constraints known so far the resulting
semidefinite program is feasible for every long-cycle fraction, for a
structural reason explained below, and the note says which local data
would change that.

## Moments as walk distributions

For the cubic graph $J$ let $NB_k$ be the matrix of non-backtracking walks
of length $k$. Then $NB_2=J^2-3I$ and $NB_{k+1}=J\,NB_k-2NB_{k-1}$, so
$J^j=\sum_k c_{j,k}NB_k$ with integer coefficients, in any cubic graph.
For instance $J^4=NB_4+7NB_2+15I$ and
$J^8=NB_8+15NB_6+95NB_4+319NB_2+543I$.

Write $Q_k(d)$ for the number of pairs $(L,F)$ of lines at line distance
$d$ with $L$ in the first and $F$ in the last pair of a non-backtracking
walk of length $k$, summed over walks and divided by 2925. For a line
eigenvalue $e$ with cosines $a_d$, the projected moment is
\[
 m_j(e)=\sum_k c_{j,k}\sum_d a_d\,Q_k(d).
\]
The identities of [opposite_boundary.md](opposite_boundary.md) give, in a
grid-free octagon with a fraction $p$ of 16-cycle pairs,

| $k$ | $Q_k(0),\dots,Q_k(4)$ |
|---|---|
| 0 | 2048, 0, 0, 0, 2048 |
| 1 | 0, 0, 12288, 0, 0 |
| 2 | 6144, 0, 0, 0, 18432 |
| 3 | 0, 0, 36864, 0, 12288 |
| 4 | 12288, 0, $8192p$, $24576-8192p$, 61440 |

The last row is the inversion of the known $m_4$ over the five line
eigenspaces. The verifier checks that the walk expansion reproduces
$m_0,\dots,m_4$ exactly.

## Walks along indexed systems

Each opposite pair $\{X,Y\}$ lies in six indexed systems, three with
$X$ as root and three with $Y$. Walking $k$ steps along the cycle of one
of them in either direction gives 12 walks per pair. For $k\ge3$ such a
walk lies in only one system, because two systems with the same root
would need two intermediate vertices containing both of their second
lines. For even $k$ the last pair shares the root with the first. Hence,
per line,
\[
 Q_k(0)\ge12288,\quad Q_k(4)\ge24576\qquad(k=4,6),
\]
and on indexed eight-cycles the walk of length eight closes up, so
\[
 Q_8(0)\ge 12288\,(2-t),
\]
where $t$ is the fraction of indexed edge incidences on cycles of length
at least twelve. At $k=4$ the first bound is an equality, so every
non-backtracking walk of length four that ends on a pair sharing a line
runs along an indexed system. For odd $k$ the three line pairs other than
the last-partner pair are at distance two, giving $Q_k(2)\ge36864$.

## The semidefinite program

The spectrum of $J$ lies in $[-3,3]$, so for each nonconstant line
eigenspace the sequence $m_0,\dots,m_8$ must be a moment sequence of a
positive measure on that interval. The verifier imposes the Hankel
condition and the three localizing conditions for $3\pm\lambda$ and
$9-\lambda^2$. The unknowns are $Q_5,\dots,Q_8$ (nonnegative, with their
fixed totals and the bounds above), $p$ and $t$. The program is feasible
for every $t\in[0,1]$.

This cannot change while $t$ enters only through the bound on $Q_8(0)$.
Raising $t$ only weakens that bound, so feasibility at $t=0$ gives
feasibility everywhere. An upper bound on $t$ needs a constraint that
grows with $t$. On an indexed cycle of length at least twelve the walk of
length eight ends at $\{X,Y''\}$ with $Y''\ne Y$, so it moves weight
$12288\,t$ to $Q_8(d^*)$ with $d^*=d(Y,Y'')$. For a twelve-cycle this is
the four-step relation read backwards, so $d^*\in\{2,3\}$. For longer
cycles $d^*$ is not known.

## What would make it bite

Two local computations would add real constraints.

1. The line distance between the partners $Y_i$ and $Y_{i+4}$ of an
   indexed cycle of length at least sixteen, from the cycle diagrams.
2. The exact value of $m_5$. Writing $J^2D^{\mathsf T}=6D^{\mathsf T}+2S$
   gives $m_5$ from $DJD^{\mathsf T}$, $RS$ and $S^{\mathsf T}JS$. The first
   two are known: per opposite pair, $RS$ has six entries at distance two
   and three at distance four. For an edge from $U$ to $V=\{C_1,C_2\}$ the
   nine distances between the third lines $Z_{12},Z_{13},Z_{23}$ of $U$
   and the third lines $C_3,Y_A,Y_B$ of $V$ are $4,2,2$ in the first
   column and $d(Z_{12},Y_A)=d(Z_{12},Y_B)=2$. Only the four distances
   between $Z_{13},Z_{23}$ and $Y_A,Y_B$ are not yet determined. Here $Y_A$
   is the third middle line of $\{A,Z_{12}\}$, at distance two in the
   cubic graph from $U$.

Run `python thin/nb_moment_sdp.py` from the repository root. The exact
checks use the standard library; the semidefinite part needs `cvxpy`
with the Clarabel solver and is skipped otherwise.
