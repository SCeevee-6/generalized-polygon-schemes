"""Decomposition of the fusion conditions for regular near octagons with parameters (s, t2, t3, t).

c_i = t_i + 1, b_i = s(t - t_i), t_1 = 0, t_4 = t. Region: s >= 1, 0 <= t2 <= t3, t >= t3 + 1 (integers; the
sign arguments below hold for reals in this region as well).

For each of the 13 proper partitions of {1,2,3,4} the closure conditions are polynomials in Q[s,t2,t3,t].
Their common zeros in the region are found by a factor-and-split elimination:
  * every equation is factored over Q; a factor that cannot vanish in the region is dropped. "Cannot vanish" is
    decided by writing s = 1+r, t3 = t2+w, t = t2+w+1+u and checking that all coefficients have one sign and the
    constant term is nonzero. If all coefficients have one sign and the constant term is zero, the factor vanishes
    only where all its monomials do, which gives branches r=0 / t2=0 / w=0 / u=0;
  * a reducible equation splits the problem into one branch per factor;
  * an equation linear in a variable v, A*v + B, is solved for v when A cannot vanish in the region, otherwise the
    problem splits into {A = 0, B = 0} and {A != 0, v = -B/A};
  * if nothing is linear, a Groebner basis (grevlex, QQ) of what is left is computed and factored again; if that
    does not help the branch is reported as a residual system.
Each branch ends as a list of solved variables plus residual irreducible equations. All branches are substituted
back into the original conditions as a check.

Usage: python near_octagon_decompose.py        (writes near_octagon_decompose_output.txt and .json)
"""
import os, sys, json, time
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from near_polygon_fusions import scheme, set_partitions, conditions

s, t, t2, t3 = sp.symbols('s t t2 t3')
r, w, u = sp.symbols('r w u')
VARS = [s, t2, t3, t]
SHIFT = {s: 1 + r, t3: t2 + w, t: t2 + w + 1 + u}
BACK = {r: s - 1, t2: t2, w: t3 - t2, u: t - t3 - 1}


def shifted(f):
    return sp.Poly(sp.expand(f.subs(SHIFT, simultaneous=True)), r, t2, w, u)


def sign_class(f):
    """'never' if f has no zero in the region, 'orthant' if one-signed without constant term, else 'free'.
    Polynomials involving a variable outside (s, t2, t3, t) are always 'free'."""
    if not sp.sympify(f).free_symbols <= {s, t, t2, t3}:
        return 'free'
    P = shifted(f)
    cs = P.coeffs()
    if all(c > 0 for c in cs) or all(c < 0 for c in cs):
        return 'never' if P.eval({r: 0, t2: 0, w: 0, u: 0}) != 0 else 'orthant'
    return 'free'


def orthant_branches(f):
    """f one-signed in shifted variables, no constant term: zero iff every monomial vanishes.
    Returns a list of sets of shifted variables that must be zero (covering all possibilities)."""
    P = shifted(f)
    gens = P.gens
    out = []

    def rec(monos, zero):
        monos = [m for m in monos if not any(m[i] > 0 and gens[i] in zero for i in range(4))]
        if not monos:
            out.append(frozenset(zero))
            return
        m = monos[0]
        for i in range(4):
            if m[i] > 0:
                rec(monos, zero | {gens[i]})
    rec(P.monoms(), frozenset())
    mins = [z for z in set(out) if not any(o < z for o in out)]
    return mins


def clean(eqs, nonzero):
    """factor, drop non-vanishing factors; returns None if the branch is dead, else list of lists of factors."""
    res = []
    seen = set()
    for e in eqs:
        e = sp.expand(e)
        if e == 0:
            continue
        if e.is_number:
            return None
        e = sp.numer(sp.together(e))
        facs = []
        for f, _ in sp.factor_list(e, *VARS)[1]:
            f = sp.expand(f)
            if f.could_extract_minus_sign():
                f = sp.expand(-f)
            if f in nonzero or sp.expand(-f) in nonzero:
                continue
            if sign_class(f) == 'never':
                continue
            if f not in facs:
                facs.append(f)
        if not facs:
            return None
        key = tuple(sorted(map(str, facs)))
        if key not in seen:
            seen.add(key)
            res.append(facs)
    return res


def solve(eqs, subs=(), nonzero=frozenset(), depth=0, log=None):
    """returns list of branches (subs, residual_eqs, nonzero)."""
    if subs:
        m = resolve(subs)
        eqs = [sp.numer(sp.together(sp.sympify(e).subs(m))) for e in eqs]
    E = clean(eqs, nonzero)
    if E is None:
        return []
    if not E:
        return [(subs, [], nonzero)]
    E.sort(key=lambda fl: (len(fl), sum(sp.count_ops(f) for f in fl)))
    # orthant-type single factors
    for i, fl in enumerate(E):
        if len(fl) == 1 and sign_class(fl[0]) == 'orthant' and not any(
                sp.degree(fl[0], v) == 1 and sp.diff(fl[0], v).is_number for v in fl[0].free_symbols):
            rest = [sp.Mul(*g) for j, g in enumerate(E) if j != i]
            out = []
            for zs in orthant_branches(fl[0]):
                out += solve(rest + [BACK[z] for z in zs], subs, nonzero, depth + 1)
            return out
    # reducible: split
    for i, fl in enumerate(E):
        if len(fl) > 1:
            rest = [sp.Mul(*g) for j, g in enumerate(E) if j != i]
            out = []
            for f in fl:
                out += solve(rest + [f], subs, nonzero, depth + 1)
            return out
    irr = [fl[0] for fl in E]
    # linear with non-vanishing coefficient
    best = None
    for f in irr:
        for v in VARS:
            if sp.degree(f, v) == 1:
                A = sp.expand(sp.diff(f, v)); B = sp.expand(f - A * v)
                sc = 'never' if A.is_number else sign_class(A)
                score = (0 if sc == 'never' else 1, sp.count_ops(A) + sp.count_ops(B))
                if best is None or score < best[0]:
                    best = (score, f, v, A, B, sc)
    if best is not None:
        _, f, v, A, B, sc = best
        val = sp.cancel(-B / A)
        rest = [sp.numer(sp.together(g.subs(v, val))) for g in irr if g != f]
        out = []
        if sc != 'never':
            out += solve([g for g in irr if g != f] + [A, B], subs, nonzero, depth + 1)
            out += solve(rest, subs + ((v, val),), nonzero | {A}, depth + 1)
        else:
            out += solve(rest, subs + ((v, val),), nonzero, depth + 1)
        return out
    # nothing linear: Groebner and retry once
    free = sorted(set().union(*[f.free_symbols for f in irr]), key=lambda x: VARS.index(x))
    G = sp.groebner(irr, *free, order='grevlex', domain=sp.QQ)
    Gex = [sp.expand(g) for g in G.exprs]
    if Gex == [1]:
        return []
    E2 = clean(Gex, nonzero)
    if E2 is None:
        return []
    if any(len(fl) > 1 for fl in E2) or any(any(sp.degree(fl[0], v) == 1 for v in fl[0].free_symbols) for fl in E2) \
            or len(E2) < len(irr):
        if sorted(map(str, [sp.Mul(*fl) for fl in E2])) != sorted(map(str, irr)):
            return solve([sp.Mul(*fl) for fl in E2], subs, nonzero, depth + 1)
    return [(subs, [fl[0] for fl in E2], nonzero)]


def resolve(subs):
    """compose the substitution chain into one map variable -> expression in the free variables."""
    m = {}
    for v, val in subs:
        m[v] = val
    for _ in range(len(subs) + 1):
        m = {v: sp.cancel(val.subs(m)) for v, val in m.items()}
    return m


def main():
    t0 = time.time()
    p = scheme(4)
    out, js = [], {}
    for part in set_partitions([1, 2, 3, 4]):
        if len(part) in (1, 4):
            continue
        part = sorted(sorted(q) for q in part)
        eqs = conditions(p, part)
        br = solve(eqs)
        comps = []
        for subs, resid, nz in br:
            m = resolve(subs)
            # check: every original condition vanishes on the branch (modulo the residual equations)
            ok = True
            if not resid:
                ok = all(sp.cancel(e.subs(m)) == 0 for e in eqs)
            comps.append({'solved': {str(k): str(sp.factor(v)) for k, v in sorted(m.items(), key=lambda kv: str(kv[0]))},
                          'residual': [str(sp.factor(g)) for g in resid],
                          'nonzero': sorted(str(sp.factor(a)) for a in nz), 'checked': bool(ok)})
        # dedupe
        uniq = []
        for c in comps:
            if c not in uniq:
                uniq.append(c)
        js[str(part)] = uniq
        out.append(f'{part}: {len(eqs)} conditions, {len(uniq)} branches  [{time.time() - t0:.0f}s]')
        for c in uniq:
            out.append('   ' + json.dumps(c))
        print(out[-1 - len(uniq)], flush=True)
    open(os.path.join(HERE, 'near_octagon_decompose_output.txt'), 'w').write('\n'.join(out) + '\n')
    json.dump(js, open(os.path.join(HERE, 'near_octagon_decompose.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
