# Maximum flag families in O(2)

The flag opposition graph has 8775 vertices, valency 4096 and least eigenvalue −64, so a family of pairwise
non-opposite flags has at most 135 members. Write N for the point–flag incidence matrix, N_L for the line–flag one,
V_θ for the θ-eigenspace of the collinearity graph (θ = 10, 5, 1, −3, −5) and W_θ for the θ-eigenspace of the line
graph (θ = 12, 7, 3, −1, −3, multiplicities 1, 351, 650, 675, 1248).

**The −64 eigenspace.** `flag_structure.py` checks that N^T V_5, N^T V_−3, N_L^T W_−1 and N_L^T W_7 are all
contained in the −64 eigenspace E, with dimensions 351, 675, 675, 351 and total rank 2052. Similarly N^T V_1 and
N_L^T W_3 lie in the 64-eigenspace (1300), N^T V_−5 in the 256-eigenspace (78), N_L^T W_−3 in the 16-eigenspace
(1248), and the all-ones vector spans the 4096-eigenspace. That accounts for 1 + 2052 + 1300 + 78 + 1248 = 4679
dimensions. On the remaining 4096 let n_λ be the multiplicities. From tr O = 0 and tr O² = 8775·4096,

    Σ n_λ = 4096,   Σ λ n_λ = 4096,   Σ λ² n_λ = 4096,

so Σ (λ−1)² n_λ = 0 and every remaining eigenvalue is 1. In particular E is exactly the 2052-dimensional space above.

**Consequence.** If F is a family of 135 pairwise non-opposite flags, Hoffman equality puts 1_F in ⟨1⟩ ⊕ E, so

    1_F(p, L) = α_p + β_L

for some real α on points and β on lines. Along a line the values α_p differ by integers, and the incidence graph is
connected, so after shifting α by a constant and β by its negative both α and β are integer valued. Conversely every
0/1 vector of this form with α ∈ ⟨1⟩ ⊕ V_5 ⊕ V_−3 and β ∈ ⟨1⟩ ⊕ W_−1 ⊕ W_7 and 135 ones is a maximum family
(its quadratic form x^T O x vanishes).

**Occupancy identities.** Let u_p be the number of flags of F at p and z_L the number on L. Then for every point
and every line

    S_2(p) − S_1(p) − 5 u_p = 5,      Z_2(L) − 3 Z_1(L) + 5 z_L = 3,

where S_i(p) is the sum of u over points at distance i from p and Z_i(L) the sum of z over lines at distance i from
L in the dual. These are (A−5)(A+3)u = 5·1 and (B+1)(B−7)z = 3·1. Both canonical families satisfy them (checked).

The two canonical families are the cases β constant (all flags at the 27 points of a line neighbourhood) and α
constant (all flags on the 45 lines of a dual point neighbourhood). A mixed family would need both α and β
non-constant. `flag_mixed_ilp.py` searches for one: fix a point p0 with one flag in and one out (the stabiliser of a
point is 2-transitive on its five lines), integer α, β in [−4, 5], the two identities, and α non-constant. HiGHS ran 3000 s without finding a solution or
proving there is none, so the classification is still open; the potential form is the reduction it should start from.

## Exact argument for every octagon of order (2,4)

The pullback calculation also follows from unique geodesics. Write O for flag
opposition, A₄ for point opposition and B₄ for line opposition. For a flag
(p,L) and a point q opposite p, the unique point of L nearest q is at distance
three. The first line on the unique geodesic from q to that point is at dual
distance three from L. The other four lines through q are opposite L, by girth
16. Exchanging points and lines gives the dual count of two. Thus

    O Nᵀ = 4 Nᵀ A₄,         O N_Lᵀ = 2 N_Lᵀ B₄.

The intersection arrays give the following eigenvalues in the adjacency
orderings above.

| | Principal | Second | Third | Fourth | Last |
|---|---:|---:|---:|---:|---:|
| A₄ | 1024 | −16 | 16 | −16 | 64 |
| B₄ | 2048 | −32 | 32 | −32 | 8 |

The two full pullback spaces intersect in the constants. Indeed, if an
assignment to points and an assignment to lines agree on all flags,
connectedness forces one common value. Consequently the nonconstant blocks
have the dimensions stated above, exactly. The trace and trace-square
calculation gives the last 4096 eigenvalues. This proves the potential
description for every GO(2,4). For the occupancy identities, point–line
incidence intertwines B with A+2I on corresponding line and point modules:
if J has rows indexed by points and columns by lines, then JJᵀ=A+5I
and JᵀJ=B+3I, so JB=(A+2I)J. Applying the distance polynomials to the
occupancies gives the constants 5 and 3. Existence of a mixed maximum
family remains open in this investigation.
