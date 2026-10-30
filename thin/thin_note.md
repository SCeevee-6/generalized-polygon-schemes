# Thin generalized octagons GO(1,t): fusions

Point graph of GO(1,t) = incidence graph of GQ(t,t); array {t+1,t,t,t;1,1,1,t+1}; integral iff 2t is a square (Feit-Higman).
Exact test (thin_fusions.py, exact Fractions) for t = 2, 8, 18, 32, 50: of the 15 partitions of {1,2,3,4}, the only closing ones are
the trivial ones plus {1,3}|{2,4} (bipartite parity) and {2}|{1,3}|{4} (two-sided halved SRG). Identical for all five t.
Not Q-polynomial for any of the five (0 of 24 orderings).
Explicit check (verify_thin_explicit.py): W(2) and W(8) incidence graphs (30 and 1170 vertices) have the predicted arrays exactly.
Halved graph = collinearity graph of GQ(t,t), SRG((t+1)(t^2+1), t(t+1), t-1, t+1): nothing beyond known objects; no primitive fusion.
