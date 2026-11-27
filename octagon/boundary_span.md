# Detecting outer components from the distance-three sphere

Fix a point $x$ in an arbitrary GO(2,4). Let $S$ be its 640 points at
distance three and let $k$ be the number of connected components on its
1024 opposite points. Put $g=(64,32,8,0,-1)$ and form the integer
matrix $B=(g_{d(u,v)})_{u,v\in S}$ using ambient distances. Then
\[
                     \operatorname{rank}_{\mathbb R} B=351-k .
\]
In particular rank 350 means connected, and ranks 349, 348 and 347 mean
two, three and four components. No distances between opposite points
of the root are needed. Here and below the proofs use the
[boundary-count theorem](README.md).

## Proof

The matrix $G=(g_{d(u,v)})$ on all points equals 320 times the orthogonal
projector onto the adjacency eigenspace $V_5$. The five normalized
cosines follow from the intersection-array recurrence at eigenvalue 5.
The projector diagonal is $351/1755=1/5$, giving the stated scaling.

The rank deficit $351-\operatorname{rank} B$ is the dimension of
the functions in $V_5$ vanishing on $S$. On the radius-two ball such
a function obeys the eigenvalue equation of the induced graph.
The positive radial vector with values $1,1/2,1/8$ on these three
layers has eigenvalue five there (the last layer has one parent and
one sibling). The ball is connected, so Perron–Frobenius gives a
one-dimensional eigenspace. Subtract a multiple of the full radial
eigenfunction $r(y)=g_{d(x,y)}/64$. It too vanishes on $S$ and agrees
with that vector on the ball.

What remains is supported on the outer shell. Its restriction to each
connected 5-regular outer component is a 5-eigenfunction and hence
constant. Orthogonality to constants gives one relation on the $k$
values. Conversely all $k-1$ component contrasts are ambient
5-eigenfunctions by the boundary-count theorem. Together with $r$
these exhaust the kernel.

## Reading a component size without its edges

For an outer point $v$, set $b=(g_{d(u,v)})_{u\in S}$ and let $B^+$
denote the Moore–Penrose inverse. If its component has size $s$, then
\[
                 64(64-b^T B^+b)=20480/s-19 .
\]
Indeed the Schur residual on the left before multiplication by 64 is
320 times the diagonal of the orthogonal projector onto the kernel
space just described. On outer points the component-contrast projector
has diagonal $1/s-1/1024$. The independent radial summand is
orthogonal to those contrasts, has squared norm 5 and takes value
$-1/64$ there. This adds $1/20480$ to the diagonal.

**Target for a subsequent proof.** To force a disconnected root, it
suffices to rule out simultaneous rank 350 at every root. In that
case every outer column at every root has the same Schur residual
$1/64$. For any nonsingular principal 350-row minor, adjoining any
one outer point therefore divides its determinant by 64 (the chosen
rows span the boundary). This gives an integer divisibility test on
boundary metrics. More strongly, with connected outer shell the residual
projector is purely radial. Its scaled off-diagonal entries between
outer points are also $1/64$. Thus any two distinct candidate outer
profiles $b,c$ must satisfy
\[
             64b^TB^+c\in\{2047,511,-1,-65\}.
\]
These four values recover their mutual distance exactly. This is a
finite compatibility test: a proposed rank-350 boundary must admit
1024 distinct profiles satisfying the diagonal equation and forming
a clique for this four-value test. Testing the diagonal equation
first can discard individual profiles before any pairs are compared.
A low modular rank alone cannot certify a low real
rank. A disconnected root would still need component recognition to
finish the uniqueness argument.

Run `python octagon/boundary_span.py` for a Ree–Tits check. At root zero
the rank modulo 101 is 349, matching the real upper bound from the two
components. The checker also tests the bordered identity modulo 101 at
two outer points. The universal formulas above follow from the proof.

## Further reduction

[Boundary reconstruction](boundary_reconstruction.md) recovers the entire
outer graph from boundary distance-two cliques, gives exact Laplacian and
modular profile tests, and couples outer adjacency to the boundary incidence
matrix via mixed moments and positive-semidefinite constraints.
