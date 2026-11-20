"""Distance-regular Cayley graphs on abelian groups for given intersection arrays.

An integral Cayley graph on an abelian group has its connection set a union of atoms
[g] = {x : <x> = <g>}. The character sums over atoms are integers, so the spectrum of the
array becomes a 0/1 model over atoms, solved with CP-SAT. Every solution is then checked
directly for distance-regularity and its array. All solutions are enumerated with no-good cuts.

usage: python abelian_cayley.py "52,42,16;1,6,28" 25,25 [--limit 500] [--tl 900]
"""
import sys, json, math, argparse, time
try:
    from ortools.sat.python import cp_model  # must load before numpy on Windows
except OSError:
    cp_model = None
import numpy as np
from fractions import Fraction
from itertools import product


def parse(s):
    b, c = s.split(';')
    return [int(x) for x in b.split(',')], [int(x) for x in c.split(',')]


def spectrum(b, c):
    """integer eigenvalues and multiplicities from the array (exact)."""
    d = len(b); k = b[0]; B = b + [0]; C = [0] + c
    a = [k - B[i] - C[i] for i in range(d + 1)]
    ks = [1]
    for i in range(d):
        ks.append(ks[-1] * B[i] // C[i + 1])
    v = sum(ks)
    T = np.zeros((d + 1, d + 1))
    for i in range(d + 1):
        T[i, i] = a[i]
        if i < d: T[i, i + 1] = B[i]
        if i > 0: T[i, i - 1] = C[i]
    th = sorted({int(round(x)) for x in np.linalg.eigvals(T).real}, reverse=True)
    assert len(th) == d + 1, 'non-integral spectrum'
    ms = []
    for t in th:
        u = [Fraction(1), Fraction(t, k)]
        for i in range(1, d):
            u.append(((t - a[i]) * u[i] - C[i] * u[i - 1]) / B[i])
        m = Fraction(v) / sum(ks[i] * u[i] ** 2 for i in range(d + 1))
        assert m.denominator == 1
        ms.append(int(m))
    return v, th, ms


def group(ns):
    E = np.array(list(product(*[range(n) for n in ns])), dtype=np.int64)
    return E


def atoms(E, ns):
    e = math.lcm(*ns)
    units = [u for u in range(1, e) if math.gcd(u, e) == 1]
    nsa = np.array(ns)
    idx = {tuple(x): i for i, x in enumerate(E)}
    lab = -np.ones(len(E), dtype=np.int64); reps = []
    for i, x in enumerate(E):
        if lab[i] >= 0: continue
        orb = {idx[tuple((u * x) % nsa)] for u in units}
        for j in orb: lab[j] = len(reps)
        reps.append(i)
    return lab, reps, e, idx


def ulm_type(x, ns):
    """Aut(G)-orbit invariant of x in the abelian group prod Z_n: for each prime p, the height
    sequence of the p-component (Baer: two elements of a finite abelian p-group are in the same
    Aut-orbit iff their height sequences agree)."""
    from sympy import factorint
    primes = sorted({q for n in ns for q in factorint(n)})
    key = []
    for q in primes:
        comp = []  # p-component of x, coordinate-wise (Z_n = Z_{q^e} x Z_{n/q^e} by CRT)
        mods = []
        for xi, n in zip(x, ns):
            e = 0; nn = n
            while nn % q == 0: nn //= q; e += 1
            if e == 0: continue
            mods.append(q ** e); comp.append(int(xi) % (q ** e))
        seq = []; y = comp[:]
        while any(y):
            h = min(_val(yi, q) for yi in y if yi)
            seq.append(h); y = [(q * yi) % mq for yi, mq in zip(y, mods)]
        key.append(tuple(seq))
    return tuple(key)


def _val(a, q):
    k = 0
    while a % q == 0: a //= q; k += 1
    return k


def is_drg(E, ns, Sidx, idx, b, c):
    v = len(E); nsa = np.array(ns)
    Svec = E[Sidx]
    nbr = np.array([[idx[tuple((x + s) % nsa)] for s in Svec] for x in E])
    dist = -np.ones(v, dtype=np.int64); dist[0] = 0; fr = [0]; dd = 0
    while fr:
        dd += 1; nf = []
        for x in fr:
            for y in nbr[x]:
                if dist[y] < 0: dist[y] = dd; nf.append(y)
        fr = nf
    D = dist.max()
    if (dist < 0).any() or D != len(b): return False, int(D)
    B = b + [0]; C = [0] + c
    for x in range(v):
        i = dist[x]; dn = dist[nbr[x]]
        if (dn == i - 1).sum() != C[i] or (dn == i + 1).sum() != B[i]: return False, int(D)
    return True, int(D)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('array'); ap.add_argument('ns')
    ap.add_argument('--limit', type=int, default=500); ap.add_argument('--tl', type=float, default=900)
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--out', default=None); ap.add_argument('--symbreak', action='store_true')
    A = ap.parse_args()
    b, c = parse(A.array); ns = [int(x) for x in A.ns.split(',')]
    v, th, ms = spectrum(b, c); k = b[0]
    assert v == math.prod(ns), (v, ns)
    E = group(ns); lab, reps, e, idx = atoms(E, ns)
    na = len(reps); size = np.bincount(lab, minlength=na)
    M = np.zeros((v, na)); M[np.arange(v), lab] = 1
    w = (e // np.array(ns))
    P = ((E[reps] * w) @ E.T) % e
    Kf = np.cos(2 * np.pi * P / e) @ M
    K = np.rint(Kf).astype(np.int64)
    assert np.abs(Kf - K).max() < 1e-6
    t0 = time.time()
    d = len(b); B_ = b + [0]; C_ = [0] + c
    a_ = [k - B_[i] - C_[i] for i in range(d + 1)]
    ks = [1]
    for i in range(d): ks.append(ks[-1] * B_[i] // C_[i + 1])
    U = []  # U[j][i] = u_i(theta_j)
    for t in th:
        u = [Fraction(1), Fraction(t, k)]
        for i in range(1, d): u.append(((t - a_[i]) * u[i] - C_[i] * u[i - 1]) / B_[i])
        U.append(u)
    Pm = [[ks[i] * U[j][i] for i in range(d + 1)] for j in range(d + 1)]
    Qm = [[ms[j] * U[j][i] for j in range(d + 1)] for i in range(d + 1)]
    for row in Pm + Qm:
        for z in row: assert z.denominator == 1, 'P or Q not integral'
    Pm = [[int(z) for z in r] for r in Pm]; Qm = [[int(z) for z in r] for r in Qm]
    Kd = (K.T * size[None, :]) / size[:, None]   # Kd[a,b] = sum over chars in atom b at an element of atom a
    Kd2 = np.rint(Kd).astype(np.int64); assert np.abs(Kd - Kd2).max() < 1e-6
    m = cp_model.CpModel()
    R = range(1, d + 1)
    x = {(aa, i): m.NewBoolVar(f'x{aa}_{i}') for aa in range(1, na) for i in R}
    y = {(bb, j): m.NewBoolVar(f'y{bb}_{j}') for bb in range(1, na) for j in R}
    for aa in range(1, na): m.AddExactlyOne([x[aa, i] for i in R])
    for bb in range(1, na): m.AddExactlyOne([y[bb, j] for j in R])
    for i in R: m.Add(sum(int(size[aa]) * x[aa, i] for aa in range(1, na)) == ks[i])
    for j in R: m.Add(sum(int(size[bb]) * y[bb, j] for bb in range(1, na)) == ms[j])
    for bb in range(1, na):
        for i in R:
            m.Add(int(K[bb, 0]) * (1 if i == 0 else 0) + sum(int(K[bb, aa]) * x[aa, i] for aa in range(1, na) if K[bb, aa])
                  == sum(Pm[j][i] * y[bb, j] for j in R))
    for aa in range(1, na):
        for j in R:
            m.Add(sum(int(Kd2[aa, bb]) * y[bb, j] for bb in range(1, na) if Kd2[aa, bb]) + int(Kd2[aa, 0]) * 0
                  == sum(Qm[i][j] * x[aa, i] for i in R))
    xs = [x[aa, 1] for aa in range(1, na)]
    if A.symbreak:
        # one-level symmetry breaking: Aut(G) is transitive on the atoms of each type, so for the
        # largest type T we may assume: if S_i meets T for the smallest such i, then the first atom of T lies in S_i.
        types = {}
        for aa in range(1, na): types.setdefault(ulm_type(E[reps[aa]], ns), []).append(aa)
        T = max(types.values(), key=len); r0 = T[0]
        meet = {i: m.NewBoolVar(f'meet{i}') for i in R}
        for i in R:
            m.AddMaxEquality(meet[i], [x[aa, i] for aa in T])
        # i* = smallest class meeting T; force x[r0, i*]
        for i in R:
            prev = [meet[j] for j in R if j < i]
            m.AddImplication(meet[i], x[r0, i]).OnlyEnforceIf([pv.Not() for pv in prev])
    sols = []; status_last = None; drgs = []
    while len(sols) < A.limit:
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = A.tl
        s.parameters.num_search_workers = A.workers
        st = s.Solve(m); status_last = s.StatusName(st)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
        sel = [aa for aa in range(1, na) if s.Value(x[aa, 1])]
        Sidx = np.where(np.isin(lab, sel))[0]
        ok, D = is_drg(E, ns, Sidx, idx, b, c)
        sols.append({'atoms': sel, 'drg': ok, 'diam': D})
        if ok: drgs.append([E[i].tolist() for i in Sidx])
        m.AddBoolOr([x[aa, 1].Not() for aa in sel] + [x[aa, 1] for aa in range(1, na) if aa not in sel])
    out = dict(array=A.array, ns=ns, v=v, theta=th, mult=ms, natoms=na, final_status=status_last,
               n_spectral_solutions=len(sols), n_drg=len(drgs), complete=(status_last == 'INFEASIBLE'),
               seconds=round(time.time() - t0, 1), drg_connection_sets=drgs[:3],
               solutions=sols[:50])
    txt = json.dumps(out)
    if A.out:
        open(A.out, 'w').write(txt + '\n')
    else:
        print(txt)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        import traceback
        open('abelian_cayley_error.log', 'a').write(' '.join(sys.argv) + '\n' + traceback.format_exc())
        raise
