# Terwilliger algebra of the GO(2,4) collinearity graph at a point

Script: `terwilliger_wedderburn.py`. Output: `out_orbits/terwilliger_wedderburn_output.txt` and `out_orbits/terwilliger_modules.csv`.

Status: computed here for the Ree-Tits octagon O(2). I haven't checked whether the Terwilliger algebra of this graph has been worked out in the literature.

## Method

Everything lives in the centraliser algebra C of H = 2F4(2)_x. C has dimension 258, with a basis indexed by the H-orbits on ordered pairs; products are computed exactly from the orbit labels.

T(x) = ⟨A, E*_0, …, E*_4⟩ is built in three steps:
1. Spanning words are found modulo a prime.
2. They are put in RREF over Q.
3. Closure under A and every E*_i is checked exactly.

This reproduces dimension 174.

The centres:
- The centre of T(x) is the exact nullspace of the commutators with A and the E*_i.
- The centre of C comes from a random integer combination of the commutator conditions. Each basis vector is then checked exactly against all 258 orbit matrices.

A central element with squarefree characteristic polynomial on the centre gives the central idempotents exactly. These are polynomials in it, one for each Q-irreducible factor.

For a simple component W with idempotent e_W:
- dim W comes from the rank of left multiplication by e_W, computed modulo p = 67108859. A rank mod p is at most the rank over Q, and the ranks add up to exactly 174 (resp. 258), so every one of them is the rank over Q.
- dim E*_i W, dim E_j W and the multiplicity m of W in C^1755 come from exact traces on C^1755: trace(e_W) = d·m, trace(E*_i e_W) = m·dim E*_i W, trace(E_j e_W) = m·dim E_j W. Here E_j are the primitive idempotents of A.

## T(x): Wedderburn decomposition

The centre has dimension 15, and all 15 components are defined over Q. So T(x) is a sum of 15 full matrix algebras M_d(Q), with Σ d² = 174, and Σ d·m = 1755.

| d = dim W | multiplicity m | endpoint | dual endpoint | diameter | dual diameter | dim E*_i W, i = 0..4 | dim E_j W, j = 0..4 | thin | dual thin |
|---:|---:|---:|---:|---:|---:|---|---|---|---|
| 5 | 1 | 0 | 0 | 4 | 4 | 1 1 1 1 1 | 1 1 1 1 1 | yes | yes |
| 3 | 4 | 1 | 1 | 2 | 2 | 0 1 1 1 0 | 0 1 1 1 0 | yes | yes |
| 4 | 5 | 1 | 1 | 3 | 3 | 0 1 1 1 1 | 0 1 1 1 1 | yes | yes |
| 3 | 30 | 2 | 1 | 2 | 2 | 0 0 1 1 1 | 0 1 1 1 0 | yes | yes |
| **7** | 40 | 2 | 1 | 2 | 3 | 0 0 1 **3 3** | 0 **2 2 2** 1 | **no** | **no** |
| 2 | 10 | 3 | 1 | 1 | 2 | 0 0 0 1 1 | 0 1 0 1 0 | yes | yes |
| 2 | 40 | 3 | 1 | 1 | 2 | 0 0 0 1 1 | 0 1 0 1 0 | yes | yes |
| **3** | 20 | 3 | 1 | 1 | 2 | 0 0 0 1 **2** | 0 1 1 1 0 | **no** | yes |
| **4** | 32 | 3 | 1 | 1 | 3 | 0 0 0 1 **3** | 0 1 1 1 1 | **no** | yes |
| **5** | 128 | 3 | 1 | 1 | 2 | 0 0 0 **2 3** | 0 1 **2 2** 0 | **no** | **no** |
| 1 | 42 | 3 | 2 | 0 | 0 | 0 0 0 1 0 | 0 0 1 0 0 | yes | yes |
| 2 | 80 | 3 | 2 | 1 | 1 | 0 0 0 1 1 | 0 0 1 1 0 | yes | yes |
| 1 | 1 | 4 | 1 | 0 | 0 | 0 0 0 0 1 | 0 1 0 0 0 | yes | yes |
| 1 | 100 | 4 | 2 | 0 | 0 | 0 0 0 0 1 | 0 0 1 0 0 | yes | yes |
| 1 | 117 | 4 | 3 | 0 | 0 | 0 0 0 0 1 | 0 0 0 1 0 | yes | yes |

(The dual eigenspaces E_j are ordered by eigenvalue 10, 5, 1, −3, −5.)

So the collinearity graph of GO(2,4) is **not thin** with respect to a point: 4 of the 15 irreducible T(x)-modules are not thin, with total multiplicity 40 + 20 + 32 + 128 = 220. These are:
- one with endpoint 2 and dimension 7;
- three with endpoint 3, of dimensions 3, 4 and 5.

Two of them, of dimensions 7 and 5, are also not dual thin. Every module with endpoint ≤ 1 is thin, and so is every module of dimension ≤ 2.

## The centraliser algebra C and the 84-dimensional gap

The centre of C has dimension 26. There are 20 Q-irreducible factors, 6 of them of degree 2, so 26 simple components over C. Each component is listed as (number of Galois conjugates, m_ρ = dimension of the simple C-module, n_ρ = dimension of the corresponding irreducible of H):

(1,1,1), (1,1,5), (1,1,10), (2,1,10), (2,1,16), (2,1,16), (1,2,10), (1,2,40) ×3, (1,3,4), (1,3,10), (2,3,10), (1,3,20), (1,4,5), (2,4,16), (1,4,40), (1,5,1), (2,5,64), (1,7,40).

Checks: Σ m_ρ² = 258 and Σ m_ρ n_ρ = 1755, counted with conjugates.

Restricting each simple C-module to T(x) (computed as trace(e_ρ e_W)/(n_ρ d_W)) gives:
- 24 of the 26 simple C-modules stay irreducible under T(x).
- Two split:
  - one C-module of dimension 2 (n = 40) becomes the two 1-dimensional T(x)-modules with endpoint 4 and dual endpoints 2 and 3;
  - one C-module of dimension 4 (n = 40) becomes the 2-dimensional module (endpoint 3, m = 40) plus those same two 1-dimensional modules.
- Many C-blocks carry isomorphic T(x)-modules. For example, the 3-dimensional module with endpoint 2 (m = 30) sits in three C-blocks, (1,3,10) and the conjugate pair (2,3,10). The 1-dimensional module with multiplicity 117 sits in five blocks.

So the gap 258 − 174 = 84 comes from two sources:
- T(x) cannot tell apart C-blocks that carry the same T(x)-module (several copies of M_d in C become one M_d in T(x));
- the two blocks above split.

T(x) does not see the H-irreducible labels n_ρ at all. Nothing here identifies the gap with the missing triple relations of `triple_relations_note.md`.
