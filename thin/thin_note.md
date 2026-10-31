# Thin generalized octagons: fusions of the point-distance scheme

Run from the repo root: `python thin/thin_fusions.py` (exact Fractions, under 1 s). Output: `thin/out_thin/thin_fusions.json`.
`python thin/verify_thin_explicit.py` builds W(2) and W(8) and checks the arrays on actual graphs (about 40 s).

Arrays. Order (1,t): {t+1,t,t,t; 1,1,1,t+1}, the incidence graph of GQ(t,t). Order (s,1): {2s,s,s,s; 1,1,1,2}.
Integrality needs 2st to be a square, so t (or s) is in {2, 8, 18, 32, 50, ...}.

All 15 partitions of the classes {1,2,3,4} were tested for closure.

| order | n | proper fusions that close | what the new classes are |
|---|---|---|---|
| (1,t), t = 2, 8, 18, 32, 50 | 2(t+1)(t^2+1) | {1,3}\|{2,4} | the bipartition: K_{N,N} and 2K_N, N = n/2 |
| | | {2}\|{1,3}\|{4} | class 2 is two copies of SRG((t+1)(t^2+1), t(t+1), t-1, t+1), the collinearity graph of GQ(t,t) |
| (2,1) | 45 | {3}\|{1,2,4} | class 3 is SRG(45,16,8,4), primitive |
| | | {2}\|{3}\|{1,4} | 3-class scheme |
| (8,1) | 5265 | none | |

For t = 2 the halved graph is SRG(15,6,1,3); for t = 8 it is SRG(585,72,7,9). These come straight from the numbers
sum over i,j in B of p^h_ij in the JSON output.

None of the seven schemes is Q-polynomial (0 of 24 orderings).

Verdict on the claim that the thin octagons of orders (1,2) and (1,8) admit proper fusions: true, and it holds for every
t tested. Both fusions are imprimitive. One is the bipartition. The other is the halved graph, which is the collinearity
graph of GQ(t,t). So nothing new comes out of the (1,t) case. The order (2,1) case does have a primitive fusion. Its
graph has the parameters of the triangular graph T(10), and T(n) is determined by its parameters for n other than 8
(Chang 1959, Hoffman 1960), so its parameters match T(10). Whether the graph is isomorphic to T(10) was not checked. It has 45 vertices, one for each flag of GQ(2,2). Order (8,1) has no proper fusion.
This agrees with Remark "Sharpness of thickness" in `fusion/fusion_note.tex`.
