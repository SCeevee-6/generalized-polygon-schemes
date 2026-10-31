# Fusions of the non-commutative orbital configurations of T (12480 and 14976 points), angle 4

Labels: **new** (not found in the sources I checked), **known**.

## Method

`search/fusion_closure.py`. For every transposition-closed union S of non-trivial orbitals (all 2^(U-1) - 1 of
them up to complement; U = 19 units for 12480, 20 for 14976) it computes the closure of the partition {0} | S | rest by
refinement inside the orbital algebra: blocks are split by the block-product coefficients and by the block of the
transposed class until nothing changes. This is the coherent closure of the graph S restricted to fusions, i.e. the
coarsest fusion in which S is a union of classes. Runs: 262143 sets (12480, each class) and 524287 sets (14976),
split in shards; outputs `max6_fusions.json`, `max7_fusions.json`, `max8_fusions.json`. On the commutative
actions (1600, 2304, 2925) the same code returns the known fusions that are closures of one symmetric relation (on
2304 points the non-symmetric rank-6 T.2 scheme is not of that form and is not returned, as expected).

## Results

- 12480 points (A6.2^2, both classes): every S gives back the full rank-24 configuration.
- 14976 points (5^2:4A4): every S gives either the full rank-25 configuration or one rank-16 fusion, which is not
  commutative (and not symmetric).

## Consequences

Let F be a fusion of rank at least 3. Some union S of F-classes is transposition-closed, non-empty and not
everything (if F has rank 3 and is not symmetric it would be a doubly regular tournament scheme, which needs an odd
number of points). F refines the closure of S, because F is a coherent configuration in which S is a union of classes.

- 12480: the closure of S is the full configuration, so F is the full configuration. The orbital configuration of T
  on 12480 points has no fusion strictly between itself and the trivial rank-2 scheme. **new**
- 14976: either F is the full configuration, or every such S has the rank-16 closure, so F refines the rank-16 fusion
  R16. A fusion of a commutative scheme is commutative, and R16 is not, so every fusion of rank >= 3 is
  non-commutative. 2F4(2) acts on these 14976 points with rank 16 (library permutation character), and its orbital
  scheme is a fusion of the T-configuration that therefore refines R16; equal ranks give R16 = the 2F4(2) orbital
  configuration. **new** as a statement; the 2F4(2) rank is **known** (ctbllib).
- In particular neither action has a commutative fusion other than the trivial one, which is consistent with the empty
  SRG search (an SRG would give a commutative rank-3 fusion).

## Assessment

A clean but small fact; it closes angle 4.
