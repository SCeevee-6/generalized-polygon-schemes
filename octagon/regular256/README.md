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
