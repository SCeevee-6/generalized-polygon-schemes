# An involution in a grid-free octagon must be point-centered

**Theorem.** Assume that a generalized octagon of order $(2,4)$ has no
nine-tip grid. Every nonidentity involution fixes precisely the incidence
ball of radius four about a point. In particular its fixed structure has
91 points, 45 lines and 135 flags. Equivalently, the existence of a
nonidentity involution with any other fixed structure forces a grid and
hence the Ree–Tits octagon.

The [fixed-tree theorem](fixed_involution.md) puts the fixed incidence
graph in the radius-four ball about a point $o$, and fixes its entire
radius-three ball. In particular two fixed lines cannot be opposite:
their three geodesics would be fixed pointwise by the
[failed-pair parity rule](failed_pair_types.md), producing a cycle in the
fixed tree. Write $g$ for the involution.

Suppose a line $M$ in the radius-three shell interchanges its two outer
points. Its only fixed point $p$ is at distance two from $o$. Choose one
of the two exchanged points and extend the edge from $M$ through it to a
nonbacktracking path of length four ending at a line $A$. The path and its
image concatenate at $M$ without backtracking, so $A,gA$ are opposite.
The path through $M$ is one of their three arms. Since $g$ swaps the
endpoint lines, the parity rule makes its arm permutation a transposition.
Thus $M$ is the only fixed middle line; call the other middle lines
$B,gB$. The unused points on the three middle lines are $p,b,gb$:
on $M$ the arm used the two interchanged points, leaving $p$.

The three unused points have mutual point distance three. Every line
through $p$ is fixed, since it lies at incidence distance one or three
from $o$. Hence the geodesics from $p$ to $b$ and $gb$ begin on the
same line $H$, which $g$ fixes. A centered triple already gives a grid.
In a failed pair the [cycle classification](opposite_edge_labels.md)
therefore forces the 16-cycle shape, with its shared initial edge at $p$.
In that shape $H$ is opposite the middle line $W$ of the geodesic between
$b$ and $gb$: the 16-cycle has an eight-edge arc from $H$ to $W$.
But $g$ reverses the unique geodesic between $b$ and $gb$, so $W$ too is
fixed. The opposite fixed lines $H,W$ contradict the fixed-tree theorem.

Consequently all 40 radius-three lines fix both their outer points. This
fixes the 80 points in the radius-four shell, in addition to $o$ and its
ten neighbors. The radius-four containment gives equality. In the notation
of the fixed-tree count this proves $T=0$ and eliminates $P=31,51,71$.

The argument admits point-centered central involutions. It places no
existence requirement on the automorphism group. In particular an
octagon with no involution is outside the scope of this fixed-element
argument.
