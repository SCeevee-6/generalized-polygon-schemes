"""Structure of maximum flag families in O(2) (135 pairwise non-opposite flags).
Checks numerically:
  * N^T maps the point eigenspaces V5, V-3 into the -64 eigenspace E of the flag opposition graph, and N_L^T maps the
    line eigenspaces W-1, W7 into E; together they span 2052 dimensions;
  * trace identities tr O = 0 and tr O^2 = 8775*4096 then force dim E = 2052 (see max_family_note.md);
  * both identities  S2(p) - S1(p) - 5u_p = 5  and  Z2(L) - 3Z1(L) + 5z_L = 3  hold for both canonical families.
Run: python flag_structure.py generators.npy"""
import sys, numpy as np, scipy.sparse as sps
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
sys.path.insert(0, '.')
from o2_build import octagon
def main(path):
    D, lines = octagon([g.astype(np.int64) for g in np.load(path)])
    nL = len(lines); I = csr_matrix((np.ones(3*nL), (np.repeat(np.arange(nL), 3), lines.ravel())), shape=(nL, 1755))
    LL = (I @ I.T).toarray(); np.fill_diagonal(LL, 0); DL = shortest_path(csr_matrix(LL > 0), unweighted=True).astype(np.int8)
    flags = np.array([(p, L) for L in range(nL) for p in lines[L]]); nF = len(flags); fp, fl = flags[:, 0], flags[:, 1]
    rows = []; cols = []
    for s in range(0, nF, 500):
        r, c = np.nonzero((D[fp[s:s+500]][:, fp] == 4) & (DL[fl[s:s+500]][:, fl] == 4)); rows.append(r+s); cols.append(c)
    O = csr_matrix((np.ones(sum(map(len, rows))), (np.concatenate(rows), np.concatenate(cols))), shape=(nF, nF))
    assert set(np.asarray(O.sum(1)).ravel().astype(int)) == {4096}
    N = csr_matrix((np.ones(nF), (fp, np.arange(nF))), shape=(1755, nF)); NL = csr_matrix((np.ones(nF), (fl, np.arange(nF))), shape=(nL, nF))
    A = (D == 1).astype(float); B = (LL > 0).astype(float)
    wA, VA = np.linalg.eigh(A); wB, VB = np.linalg.eigh(B)
    sub = lambda w, V, t: V[:, np.abs(w - t) < 1e-6]
    blocks = [N.T @ sub(wA, VA, 5), N.T @ sub(wA, VA, -3), NL.T @ sub(wB, VB, -1), NL.T @ sub(wB, VB, 7)]
    for X in blocks: assert np.allclose(O @ X, -64 * X)
    M = np.hstack(blocks); sv = np.linalg.svd(M, compute_uv=False)
    print('pulled-back -64 space: dims', [X.shape[1] for X in blocks], 'rank', int((sv > 1e-8*sv[0]).sum()))
    fid = {(int(p), int(L)): i for i, (p, L) in enumerate(flags)}
    def fam(pairs):
        x = np.zeros(nF, int); x[[fid[(int(p), int(L))] for p, L in pairs]] = 1; return x
    inc = lambda p: np.flatnonzero(I[:, p].toarray().ravel())
    fam_line = fam([(p, M) for p in np.flatnonzero((D[lines[0]] <= 1).any(0)) for M in inc(p)])
    fam_pt = fam([(p, M) for M in np.flatnonzero((D[0][lines] <= 1).any(1)) for p in lines[M]])
    for name, x in (('line-centred', fam_line), ('point-centred', fam_pt)):
        u = N @ x; z = NL @ x
        print(name, int(x.sum()), 'independent', int(x @ (O @ x)) == 0,
              'point identity', set(((D == 2).astype(int) @ u - (D == 1).astype(int) @ u - 5*u).astype(int)),
              'line identity', set(((DL == 2).astype(int) @ z - 3*((DL == 1).astype(int) @ z) + 5*z).astype(int)))
if __name__ == '__main__': main(sys.argv[1])
