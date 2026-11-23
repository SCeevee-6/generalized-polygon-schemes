# First stage of the regular 256-component exclusion, rerun

`check256_first_stage.py` (in `octagon/`) redoes the double-cover stage with separate code: signings are chord bits on
the BFS spanning tree used elsewhere in `octagon/`, the 1920 automorphisms act through their action on fundamental
cycles, and eigenvalues are computed in floating point (tolerance 1e-9) rather than by Descartes counting.

| quantity | PROOF.md | rerun |
|---|---:|---:|
| orbits of all 2^25 switching classes (Burnside) | 19,132 | 19,132 |
| orbits with every signed eigenvalue in [−4, −2] ∪ [−1, 69/20] | 245 | 245 (330,719 classes) |
| orbits passing the edge-tree capacity test (≤ 8 visits per double-cover vertex) | 48 | 48 |
| signings in those orbits | 27,979 | 27,979 |

The two unconditional lemmas behind the window were also checked in O(2): the Gram matrix of the 1024 opposite
points in the −5 representation has spectrum 4, 12^45, 15^32, both 512-point components have least eigenvalue
−3.9624 > −13/3, and every boundary point lies on exactly two lines of each component. The boundary count itself has a
short proof (octagon/README.md, §3).

This earlier floating-point check covered the first stage. The [subsequent independent implementation](README.md) covers the first stage with exact arithmetic and checks both later stages.
