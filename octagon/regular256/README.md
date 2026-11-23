# Reproduce the regular 256-component exclusion

Start with PROOF.md. The code reconstructs every search input from the abstract 16-vertex base. Python 3, NumPy and a little-endian C++17 toolchain with 32-bit int are sufficient. The temporary orbit table needs 128 MiB on disk. Run from the extracted directory.

    export OPENBLAS_NUM_THREADS=1
    python base.py
    g++ -std=c++17 -O3 orbits.cpp -o orbits
    ./orbits
    python spectral.py
    python towers.py
    g++ -std=c++17 -O3 fours.cpp -o fours
    ./fours cyc_*.dat
    python finish.py

Expected spectral counts are 19132, 245, 48 and 27979. The final counts are recorded in result.txt. A check using only the 16+152 four-sheeted records can start directly with finish.py using the supplied four.txt. The earlier intermediate and final counts are reproducible by starting at base.py. All final bipartiteness checks use graph traversal.

## Independent verification

From the repository root, with Python 3, NumPy and python-sat (CaDiCaL 1.9.5 backend), run

    OPENBLAS_NUM_THREADS=1 python octagon/regular256/check_quotients.py
    python octagon/regular256/check_terminal.py

The first program constructs the 1,920 base automorphisms and directly marks all
images of each unseen signing. It computes the signed characteristic polynomials
by the Berkowitz recurrence on trailing principal submatrices. Signed 64-bit
coefficients and products suffice since an order-j principal submatrix has
coefficients bounded by binom(j,k) 5^k in degree j−k, and the walk powers have
norm at most 5^k. The convolution sums are bounded by 16·2^16·5^16 < 2^63.
Endpoint shifts use Python integers. The double-quotient capacity is checked
by walking the two height-three trees themselves.

For each first signing, the four-sheet check associates an upper-bit parity
literal to each oriented walk. For cyclic quotients its constant is the upper
bit of the signed sum of the low edge voltages modulo four. In each bin of
fixed destination and low bit, the CNF forbids every five visits from having
the same upper bit. SAT enumerates assignments to the 25 chord bits, blocking
each complete assignment after it is read. The elementary case also imposes
the first-quotient test on its other two characters and the published ordering
convention. Only after enumeration does the program compare the stored records.

    double orbits, spectral, capacity 19132 245 48
    double signings passing capacity 27979
    four-sheet records [(0, 16), (1, 152)]
    both record lists agree

The second program reads the four-sheet records, uses a depth-first spanning
tree, solves the saturated-bin equations in GF(2), enumerates the free bits,
and recounts visits on each actual lift. It skips the cyclic spectral test,
so 152 cyclic eight-sheet records reach the last test. None lifts to a
16-sheet record. Counts allow repeated graphs from different quotient descriptions.

    quotient four eight sixteen bipartite
    elementary 16 40 40 40
    cyclic 152 152 0 0

The earlier floating-point double-stage check is documented in
[check_first_stage.md](check_first_stage.md). The original exact program and
these separate programs share the necessary spectral and capacity lemmas and
the input numbering. The independent enumeration does not import the original
search code. The original C++ orbit flag in fours.cpp now uses an explicit
one-time mark. Regenerating four.txt gave the same bytes.

## Nonregular double towers

`double_towers.py` allows every one of the 49 independent binary voltages
over each admissible first double quotient. It excludes every degree-16
trace cover with 2-group monodromy. The [proof](PROOF.md) gives the
subgroup argument and the exact eigenvalue certificate used in the search.

    OPENBLAS_NUM_THREADS=1 python octagon/regular256/double_towers.py
    # Last line: total 632 66 106 64 64

A second run with `--sat-terminal` uses breadth-first gauges and SAT for
the final splits in place of affine enumeration. All 48 rows agree.
