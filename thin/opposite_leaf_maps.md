# A repeat among eight outputs forces a grid

This test applies with no automorphism hypothesis. All distances in the
construction are incidence distances. For a point $a$ on a line $A$, let
$X(a,A)$ be the lines reached by a three-edge path from $a$ whose first
line differs from $A$. Girth makes these paths unique, giving
$4\cdot2\cdot4=32$ leaves.

Fix a line $A=\{u,a,a'\}$. Choose a line $K\ne A$ through $u$, then a
point $v\ne u$ on $K$. There are eight choices of a point $c$ collinear
with $v$ off $K$, two on each of four lines $C$ through $v$.

**The bijections.** For each $D\in X(a,A)$ there is a unique
$f_{a,c}(D)\in X(c,C)$ at line-graph distance two from $D$.
Indeed $D$ and $C$ are opposite via $A,K$; the eight-edge path uses $v$
on $C$. The unique arm starting at $c$ has a middle line $B$ at
incidence distance three from $c$ and four from $D$. Conversely a line
in $X(c,C)$ at line-graph distance two from $D$ supplies precisely this
arm: the two connections from $B$ towards $D,C$ must use distinct points
on $B$, since $D,C$ are opposite. Exchanging the two sides gives the
inverse map. The same construction at $a'$ gives $f_{a',c}$.

There is also a constraint on prefixes of these bijections. Suppose the
paths to distinct source leaves share $k$ initial edges from $a$, and
the image paths share $l$ from $c$. Then $k+l\le2$. Join the two
source leaves along their tree, join their images along the other tree,
and include the two four-edge connecting paths. The resulting cyclic
walk has length $(6-2k)+(6-2l)+8$. At every join, the tree path uses the
attachment toward $a$ or $c$, whereas the connecting path uses another
attachment (the opposite-line arms through $A,C$ guarantee this).
The walk is cyclically nonbacktracking and hence has length at least 16.
In particular every choice of a length-$k$ source prefix and a
length-$(3-k)$ target prefix occurs exactly once. There are 32 such
prefix pairs for each $k=0,1,2,3$.

**The repeat test.** Fix $D\in X(a,A)$. For each $c$ form the output
$$E_c=f_{a',c}^{-1}(f_{a,c}(D))\in X(a',A).$$
If any two of the eight outputs agree, the octagon has a grid.

To see this, write $B_c=f_{a,c}(D)$. Two distinct $c,d$ yield distinct
$B_c,B_d$: the length-five geodesic from $v$ to $B_c$ specifies $c$.
If $E_c=E_d=E$, then $D,E$ are opposite via $A$, and their three
middle lines are $A,B_c,B_d$. The unused point on $A$ is $u$. For
$A,B_c$, the third middle line besides $D,E$ is $C$. Its arm begins
$u,K,v,C,c$ and ends at $B_c$ through the unused point on $B_c$.
Thus the geodesics from $u$ to the two other unused points of the pair
$D,E$ share their first two edges $u,K,v$. Of the three shapes of a
mutually distance-three triple, only the centered one allows those two
shared edges. This gives the grid for $D,E$. Its center line is the
third vertex after $u$ along either geodesic, so the two points $c,d$
also belong to the same $C$. [Middle-tip criterion](opposite_seed.md)

**Exact count across blocks.** Varying $K$ and $v$ gives eight blocks of
eight outputs for the same source leaf $D$. Every one of the 32 possible
$E$ occurs exactly twice in these 64 places. For $D,E$ the two middle
lines other than $A$ each supply one place: for such a middle line $B$,
the third arm between $A,B$ (besides those through $D,E$) starts at
$u$, and uniquely specifies $K,v,C,c$. The preceding construction
recovers $B$ and $E$ in this place.

In particular a grid-free octagon must have eight distinct outputs in
every block. The across-block count alone permits this. In the Ree–Tits
octagon every opposite pair has a grid, so each of its two places lies
in the same block as the other, and each block has four doubled outputs.
Run `python thin/opposite_leaf_maps.py` for a model check at one fixed
flag, including all 32 source leaves, all eight blocks, and all 128
bijections and their prefix constraints. [Model verifier](opposite_leaf_maps.py)
