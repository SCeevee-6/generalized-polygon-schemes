# Coherent configurations of T = 2F4(2)' on cosets of its maximal subgroups

For a transitive group the orbitals already form a coherent configuration, so the orbital configuration is its own
coherent closure. On cosets of a maximal subgroup M it is an association scheme in the wide sense. It is commutative
iff the permutation character 1_M^T is multiplicity-free, its rank is <pi,pi> = sum m_chi^2, and the number of
self-paired orbitals is sum nu(chi) m_chi, where nu is the Frobenius-Schur indicator.

How it was computed (`gap/tits_maxes.g`, GAP 4.16.1, ctbllib 1.3.11, atlasrep 2.1.12):

1. M comes from the ATLAS straight-line program TF42G1-max{i}W1 applied to the standard generators of
   TF42G1-p1755B0, and |M| is checked.
2. pi(g_k) = |C_T(g_k)| |g_k^T ∩ M| / |M|. The class of each element of M is fixed by order, centralizer order
   and an explicit conjugacy test against the ATLAS class word.
3. pi is decomposed against the library table, and also matched against the induced characters of the library
   tables of maxes.
4. Independent group count: rank = number of (M,M) double cosets; MgM is self-paired iff g^-1 is in MgM.

In every case the double-coset rank equals <pi,pi>, and the self-paired count equals the indicator sum.

| max | M | degree | rank | commutative | max mult. | self-paired | symmetric | constituents |
|---|---|---|---|---|---|---|---|---|
| 1 | L3(3):2 | 1600 | 4 | yes | 1 | 4 | yes | 1 + 351 + 624 + 624' |
| 2 | L3(3):2 | 1600 | 4 | yes | 1 | 4 | yes | 1 + 351 + 624 + 624' |
| 3 | 2.[2^8].5.4 | 1755 | 5 | yes | 1 | 5 | yes | 1 + 78 + 351 + 650 + 675 |
| 4 | L2(25) | 2304 | 7 | yes | 1 | 5 | no | 1 + 27 + 27' + 351 + 624 + 624' + 650 |
| 5 | 2^2.[2^8].S3 | 2925 | 6 | yes | 1 | 6 | yes | 1 + 351 + 624 + 624' + 650 + 675 |
| 6 | A6.2^2 | 12480 | 24 | no | 3 | 16 | no | 13 distinct, 351 three times, 675 twice |
| 7 | A6.2^2 | 12480 | 24 | no | 3 | 16 | no | same as max 6 |
| 8 | 5^2:4A4 | 14976 | 25 | no | 2 | 17 | no | 13 distinct, 351, 675, 1300, 1300' twice each |

The 351 in every row is the rational 351 (library X.9). The pairs 27/27' and 624/624' are algebraically conjugate.
The full data is in `gap/maxes_coherent_configurations.csv` (library character names X.i) and `gap/tits_maxes_output.txt`.

Summary:

- Five actions are commutative: degrees 1600 (two classes), 1755, 2304 and 2925.
- Four of those are symmetric. The L2(25) action on 2304 points is commutative but not symmetric: it has one pair of
  non-self-paired orbitals, which matches the complex pair 27, 27'.
- Three actions give non-commutative configurations: the two A6.2^2 classes (12480 points, rank 24, the algebra has
  centre of dimension 13) and 5^2:4A4 (14976 points, rank 25, centre of dimension 13).
- Only the 1755-point action is metric (the GO(2,4) distance scheme). The other actions were not tested for fusions.

Novelty: Breuer and Lux (Comm. Algebra 24 (1996) 2293-2316, doi:10.1080/00927879608825701) classified the
multiplicity-free permutation characters of the sporadic groups. I have not checked whether their tables include T,
so the commutative cases above are not claimed as new. The non-commutative cases are a direct computation from
standard data.
