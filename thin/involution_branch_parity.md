# Parity on the points opposite an involution center

Suppose a grid-free generalized octagon of order $(2,4)$ has an
involution. By [the fixed-tree theorem](fixed_involution.md) it fixes all
elements within incidence distance three of a center point $p$, and each
of the 40 lines at distance three has either two fixed children or two
swapped children at distance four. Call a line of the first kind selected.

Divide the distance-three lines into five branches according to the first
line through $p$. Each branch has eight lines. Write $q_i$ for its
selected count and $Q=\sum_i q_i$. The total fixed point count is
$11+2Q$. The following table lists the branch counts that satisfy the
necessary conditions proved below, up to permutation. The last three
columns count points opposite $p$ by their number of geodesics from $p$
using selected lines.

| Fixed points | Selected lines in the five branches | One | Three | Five |
|---:|---|---:|---:|---:|
| 31 | $(2,2,2,2,2)$ | 960 | 0 | 64 |
| 51 | $(0,4,4,4,8)$ | 256 | 768 | 0 |
| 51 | $(2,2,4,6,6)$ | 288 | 704 | 32 |
| 51 | $(2,4,4,4,6)$ | 304 | 672 | 48 |
| 51 | $(4,4,4,4,4)$ | 320 | 640 | 64 |
| 71 | $(6,6,6,6,6)$ | 0 | 640 | 384 |
| 91 | $(8,8,8,8,8)$ | 0 | 0 | 1024 |

## The parity condition

For $y$ opposite $p$, let $b_i(y)=1$ if the geodesic from $p$ to $y$
in branch $i$ uses a selected line, and zero otherwise. Each of the five
branches supplies exactly one geodesic. The point $y$ is exchanged with
an opposite point $y'$ by the fixed-tree argument. Its five geodesics
to $y'$ form an involution set, so an odd number are preserved.

The midpoint of each preserved geodesic from $y$ to $y'$ is fixed.
It is at distance four from both $y$ and $p$, so it is the midpoint
of the geodesic from $p$ to $y$ in some branch and its parent line is
selected. Conversely, a fixed midpoint on a geodesic from $p$ to $y$
has exactly one fixed incident line, namely its inward line. The
geodesic segment toward $y$ and its image therefore go through
different lines and concatenate to a nonbacktracking path of length
eight from $y$ to $y'$. This is a preserved geodesic by girth. Consequently
$\sum_i b_i(y)$ is odd for every $y$ opposite $p$.

## Even branch counts

Fix a point $z$ at distance six from $p$ in branch $i$. Its eight
outward point neighbors are opposite $p$. The bit $b_i$ is constant
on these eight neighbors. For each other branch $j$, the eight neighbors
use its eight distance-three lines exactly once.

To see the latter assertion, take a distance-three line $H$ in
branch $j$. The nine-edge path from $z$ through $p$ to $H$ is
nonbacktracking. A path of length at most five would give a cycle of
length at most fourteen. Thus $d(z,H)=7$. The unique seven-edge
geodesic starts outward from $z$: the inward line through $z$ is
opposite $H$ by the eight-edge path through $p$. Its first point
after $z$ is an outward neighbor $y$ of $z$, with $d(y,H)=5$. The
path from $y$ to $H$ uses a child of $H$, since going through its
parent would give $d(p,y)\le6$. Thus $H$ occurs on the geodesic to $y$
in branch $j$. Uniqueness of the geodesic from $z$ to $H$ prevents a
second such neighbor from using $H$.

Sum the odd parity relation across the eight neighbors of $z$.
The contribution of branch $i$ is even, so $Q-q_i$ is even. Since
$Q\in\{10,20,30,40\}$, each $q_i$ is even.

## Counts of individual patterns

A given distance-three line $H$ occurs on geodesics to 128 opposite
points, by the choices along the five remaining edges. Two such lines
in different branches occur together for 16 opposite points. Indeed,
choose one of the two children of each line. The two children are
opposite via $p$. Of their five geodesics, one has midpoint $p$ and
the other four have midpoints opposite $p$: any two of the geodesics
form a 16-cycle, whose half-cycles are geodesic by girth. The four
choices of child pairs give disjoint sets of four midpoints by uniqueness
of the geodesic from $p$ in each branch.

Thus the single and pair counts of bits equal to one are $128q_i$ and
$16q_iq_j$. Put $s_i=(-1)^{b_i}$, so $\prod_i s_i=-1$. On the 16 odd
patterns, the constant, five signs, and ten sign products are an
orthogonal basis. To see orthogonality, a product of two distinct basis
functions is the product of a nonempty proper subset of the five signs.
Its sum over odd patterns vanishes by flipping one sign inside and
one outside that subset. The resulting inversion gives the count of
the pattern $s$ as
\[
\frac{1}{16}\left(
1024 + \sum_i s_i(1024-256q_i)
+ \sum_{i<j}s_is_j(1024-256(q_i+q_j)+64q_iq_j)
\right).
\]
Checking nonnegativity for each odd pattern gives the table. This can
be reproduced with integer arithmetic.

~~~python
from itertools import combinations, combinations_with_replacement, product

for q in combinations_with_replacement(range(0, 9, 2), 5):
    if sum(q) not in (10, 20, 30, 40):
        continue
    counts = [0] * 6
    for bits in product((0, 1), repeat=5):
        if sum(bits) % 2 == 0:
            continue
        sign = [1 - 2 * bit for bit in bits]
        value = 1024 + sum(sign[i] * (1024 - 256 * q[i]) for i in range(5))
        value += sum(
            sign[i] * sign[j] * (1024 - 256 * (q[i] + q[j]) + 64 * q[i] * q[j])
            for i, j in combinations(range(5), 2)
        )
        if value < 0 or value % 16:
            break
        counts[sum(bits)] += value // 16
    else:
        print(11 + 2 * sum(q), q, counts[1::2])
~~~

These conditions leave seven distributions and four total fixed-point
counts. They apply to counterexamples with an involution.
