"""Orbital schemes of 2F4(2)' on cosets of its maximal subgroups.

Input: a text file with lines "SIZE |G| |M| n", "GA [..]", "GB [..]" (images of the two standard
generators in the action on cosets of M, written by gap/coset_actions.g).
Everything below is recomputed from the two permutations alone:
  * stabiliser orbits of point 0 (Schreier generators, so the orbitals are exact),
  * the full orbital matrix L[x,y] and the integer intersection numbers p^k_ij (A_i A_j = sum_k p^k_ij A_k),
  * commutativity, pairing, centre dimension of the orbital algebra (exact, sympy),
  * for commutative schemes: eigenmatrix P, multiplicities, Krein parameters (floating point),
    P- and Q-polynomial tests, and closure of every fusion (exact integer test).
Usage: python coset_schemes.py cosets/max5.txt out_cosets/max5
"""
import sys, re, json, itertools, ast
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def read_action(path):
    txt = open(path).read()
    size = [int(x) for x in re.search(r"SIZE ([0-9 ]+)", txt).group(1).split()]
    perms = []
    for tag in ("GA", "GB"):
        m = re.search(tag + r" \[([^\]]*)\]", txt, re.S)
        perms.append(np.array([int(v) - 1 for v in m.group(1).replace("\n", " ").split(",")], dtype=np.int64))
    return size, perms


def orbitals(gens, n):
    # BFS tree from point 0 with transversal t_x (t_x[0] = x), stored as int16/int32 rows.
    dt = np.int16 if n < 32767 else np.int32
    T = np.empty((n, n), dtype=dt)
    Tinv = np.empty((n, n), dtype=dt)
    ident = np.arange(n)
    T[0] = ident; Tinv[0] = ident
    seen = np.zeros(n, bool); seen[0] = True
    order = [0]; head = 0
    while head < len(order):
        x = order[head]; head += 1
        for g in gens:
            y = int(g[x])
            if not seen[y]:
                seen[y] = True
                t = g[T[x].astype(np.int64)]           # t_x then g
                T[y] = t
                inv = np.empty(n, np.int64); inv[t] = ident; Tinv[y] = inv
                order.append(y)
    assert seen.all(), "action not transitive"
    # stabiliser orbits from all Schreier generators s = t_x * g * t_{xg}^{-1}
    comp = np.arange(n)
    nschreier = 0
    for x in range(n):
        tx = T[x].astype(np.int64)
        for g in gens:
            y = int(g[x])
            s = Tinv[y].astype(np.int64)[g[tx]]
            assert s[0] == 0
            nschreier += 1
            if np.array_equal(comp[s], comp):
                continue
            k = comp.max() + 1
            A = coo_matrix((np.ones(n), (comp, comp[s])), shape=(k, k))
            _, lab = connected_components(A, directed=True, connection="weak")
            comp = lab[comp]
    # canonical labels: class 0 = {0}, others by (size, min point)
    roots = {}
    for y in range(n):
        roots.setdefault(int(comp[y]), []).append(y)
    blocks = sorted(roots.values(), key=lambda b: (b != [0], len(b), b[0]))
    label = np.empty(n, np.int64)
    for i, b in enumerate(blocks):
        label[b] = i
    L = np.empty((n, n), dtype=np.int8)
    for x in range(n):
        L[x] = label[Tinv[x].astype(np.int64)]
    del T, Tinv
    return L, [len(b) for b in blocks], [b[0] for b in blocks], nschreier


def intersection_numbers(L, reps):
    r = len(reps)
    p = np.zeros((r, r, r), dtype=np.int64)          # p[k,i,j] = p^k_ij
    row0 = L[0].astype(np.int64)
    for k, y in enumerate(reps):
        col = L[:, y].astype(np.int64)
        p[k] = np.bincount(row0 * r + col, minlength=r * r).reshape(r, r)
    return p


def center_dim(p):
    import sympy
    r = p.shape[0]
    rows = []
    for j in range(r):
        for k in range(r):
            rows.append([int(p[k, i, j] - p[k, j, i]) for i in range(r)])
    M = sympy.Matrix(rows)
    return r - M.rank()


def wedderburn_types(rank, zdim):
    # multisets of block sizes m_1..m_z with sum m^2 = rank
    out = []
    def rec(left, parts, maxm):
        if len(parts) == zdim:
            if left == 0:
                out.append(tuple(parts))
            return
        for m in range(min(maxm, int(left ** 0.5)), 0, -1):
            rec(left - m * m, parts + [m], m)
    rec(rank, [], rank)
    return out


def eigenmatrix(p, n):
    r = p.shape[0]
    B = [p[:, l, :].astype(float) for l in range(r)]        # B_l[h,j] = p^h_lj
    rng = np.random.default_rng(1)
    C = sum(c * b for c, b in zip(rng.normal(size=r), B))
    w, V = np.linalg.eig(C)
    P = np.zeros((r, r), complex)
    for i in range(r):
        v = V[:, i]; piv = np.argmax(abs(v))
        for l in range(r):
            P[i, l] = (B[l] @ v)[piv] / v[piv]
            assert np.allclose(B[l] @ v, P[i, l] * v, atol=1e-7)
    i0 = int(np.argmin(np.abs(P - np.array([float(np.sum(p[0, l, :])) for l in range(r)])).sum(1)))
    P = P[[i0] + [i for i in range(r) if i != i0]]
    Q = n * np.linalg.inv(P)
    m = Q[0].real
    return P, Q, m


def krein(P, Q, n):
    r = P.shape[0]
    q = np.einsum("li,lj,kl->kij", Q, Q, P) / n             # q[k,i,j] = q^k_ij
    return q


def polynomial_orderings(T, tol):
    """T[k,i,j] structure constants (p or q). Return the j with which the scheme is P/Q-polynomial."""
    r = T.shape[0]
    good = []
    for c in range(1, r):
        layers = [[0], [c]]
        seen = {0, c}
        ok = True
        while len(seen) < r:
            nxt = sorted({k for k in range(r) if k not in seen for j in layers[-1] if abs(T[k, c, j]) > tol})
            if len(nxt) != 1:
                ok = False; break
            layers.append(nxt); seen.add(nxt[0])
        if not ok:
            continue
        order = [l[0] for l in layers]
        far_zero = all(abs(T[order[h], c, order[j]]) <= tol for h in range(r) for j in range(r) if abs(h - j) > 1)
        step = all(abs(T[order[j + 1], c, order[j]]) > tol for j in range(r - 1))
        if far_zero and step:
            good.append(order)
    return good


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def fusion_closed(p, blocks):
    bl = [[0]] + blocks
    for A in bl:
        for B in bl:
            for C in bl:
                vals = {int(sum(p[h, i, j] for i in A for j in B)) for h in C}
                if len(vals) != 1:
                    return False
    return True


def fused_info(p, part, n):
    bl = [[0]] + sorted(sorted(B) for B in part)
    r = len(bl)
    f = np.zeros((r, r, r), dtype=np.int64)
    for c, C in enumerate(bl):
        for a, A in enumerate(bl):
            for b, B in enumerate(bl):
                f[c, a, b] = sum(p[C[0], x, y] for x in A for y in B)
    P, Q, m = eigenmatrix(f, n)
    q = krein(P, Q, n)
    info = dict(fused_rank=r, fused_valencies=[int(f[0, a, :].sum()) for a in range(r)],
                fused_multiplicities=[int(round(x)) for x in m],
                fused_P_polynomial=polynomial_orderings(f.astype(float), 0.5),
                fused_Q_polynomial=polynomial_orderings(q, 1e-8))
    if r == 3:
        info["srg"] = [n, int(f[0, 1, :].sum()), int(f[1, 1, 1]), int(f[2, 1, 1])]
    return info


def main(inp, outp):
    (gord, mord, n), gens = read_action(inp)
    assert all(len(g) == n for g in gens) and gord // mord == n
    L, val, reps, nsch = orbitals(gens, n)
    r = len(val)
    p = intersection_numbers(L, reps)
    pair = [int(L[y, 0]) for y in reps]
    # sanity: p^0_ij = delta(j, i*) k_i ; sum_j p^k_ij = k_i
    for i in range(r):
        for j in range(r):
            assert p[0, i, j] == (val[i] if j == pair[i] else 0)
        for k in range(r):
            assert p[k, i, :].sum() == val[i]
    comm = bool(all(np.array_equal(p[:, i, j], p[:, j, i]) for i in range(r) for j in range(r)))
    zdim = center_dim(p)
    res = dict(input=inp, n=n, group_order=gord, subgroup_order=mord, rank=r, valencies=val,
               pairing=pair, self_paired=sum(1 for i in range(r) if pair[i] == i),
               symmetric=all(pair[i] == i for i in range(r)), commutative=comm, centre_dim=zdim,
               wedderburn_types=wedderburn_types(r, zdim), schreier_generators=nsch)
    with open(outp + "_intersection_numbers.csv", "w") as f:
        f.write("k,i,j,p_k_ij\n")
        for k in range(r):
            for i in range(r):
                for j in range(r):
                    if p[k, i, j]:
                        f.write(f"{k},{i},{j},{p[k, i, j]}\n")
    if comm:
        P, Q, m = eigenmatrix(p, n)
        q = krein(P, Q, n)
        res["multiplicities"] = [int(round(x)) for x in m]
        res["eigenmatrix_P"] = [[[round(z.real, 6), round(z.imag, 6)] for z in row] for row in P]
        res["min_abs_nonzero_krein"] = float(np.min(np.abs(q[np.abs(q) > 1e-6])))
        res["max_abs_imag_krein"] = float(np.abs(q.imag).max())
        res["krein_negative"] = int((q.real < -1e-8).sum())
        res["zero_krein"] = [[k, i, j] for k in range(r) for i in range(r) for j in range(r) if abs(q[k, i, j]) < 1e-8]
        res["P_polynomial_orderings"] = polynomial_orderings(p.astype(float), 0.5)
        res["Q_polynomial_orderings"] = polynomial_orderings(q, 1e-8)
        fus = []
        for part in set_partitions(list(range(1, r))):
            if len(part) in (1, r - 1):
                continue
            star = all(sorted(pair[i] for i in B) in [sorted(C) for C in part] for B in part)
            if fusion_closed(p, part):
                fus.append(dict(blocks=sorted(sorted(B) for B in part), transpose_closed=star, **fused_info(p, part, n)))
        res["proper_partitions_tested"] = sum(1 for part in set_partitions(list(range(1, r))) if len(part) not in (1, r - 1))
        res["closed_fusions"] = fus
    json.dump(res, open(outp + ".json", "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("n", "rank", "valencies", "self_paired", "commutative", "centre_dim")}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
