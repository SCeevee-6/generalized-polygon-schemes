"""Exact fusion test for point schemes of generalized 2d-gons of order (s,t), d = 3, 4.

The intersection numbers p^k_ij are computed as rational functions of s, t from the intersection array
{s(t+1), st, ..., st; 1, ..., 1, t+1} (no eigenvalues, no square roots). For each partition of the
non-trivial classes the closure conditions are polynomial equations in s, t; their ideal is reduced with a
Groebner basis over QQ (sympy.groebner, exact arithmetic) and the zero set is described by factoring the
basis. Output: fusion/certify_fusions_output.txt
"""
import itertools, sys
import sympy as sp

s, t, x = sp.symbols("s t x")


def scheme(d):
    b = [s * (t + 1)] + [s * t] * (d - 1)
    c = [None] + [sp.Integer(1)] * (d - 1) + [t + 1]
    k = b[0]
    a = [0] + [k - b[i] - c[i] for i in range(1, d)] + [k - c[d]]
    v = [sp.Integer(1), x]
    for i in range(1, d):
        v.append(sp.expand((x * v[i] - a[i] * v[i] - b[i - 1] * v[i - 1]) / c[i + 1]))
    minpoly = sp.Poly(sp.expand(x * v[d] - a[d] * v[d] - b[d - 1] * v[d - 1]), x)

    def coords(poly):
        rem = sp.Poly(poly, x).rem(minpoly)
        out = [0] * (d + 1)
        for kk in range(d, -1, -1):
            lc = sp.Poly(v[kk], x).coeff_monomial(x ** kk)
            co = sp.cancel(rem.coeff_monomial(x ** kk) / lc)
            out[kk] = co
            rem = sp.Poly(sp.expand(rem.as_expr() - co * v[kk]), x)
        assert rem.is_zero
        return out

    p = {}
    for i in range(d + 1):
        for j in range(d + 1):
            cs = coords(sp.expand(v[i] * v[j]))
            for kk in range(d + 1):
                p[kk, i, j] = sp.factor(cs[kk])
    return p


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def conditions(p, d, part):
    bl = [[0]] + part
    eqs = set()
    for A in bl:
        for B in bl:
            for C in bl:
                vals = [sp.together(sum(p[h, i, j] for i in A for j in B)) for h in C]
                for v_ in vals[1:]:
                    num = sp.numer(sp.together(v_ - vals[0]))
                    num = sp.expand(num)
                    if num != 0:
                        eqs.add(num)
    return sorted(eqs, key=sp.default_sort_key)


def main(out="certify_fusions_output.txt"):
    lines = []
    for d in (3, 4):
        p = scheme(d)
        # sanity: p^0_ii = k_i, and the d=4 values at (2,4) match the 1755-point scheme
        if d == 4:
            kv = [p[0, i, i].subs({s: 2, t: 4}) for i in range(5)]
            assert kv == [1, 10, 80, 640, 1024], kv
        lines.append(f"d = {d}  (generalized {2*d}-gon), valencies " + ", ".join(str(sp.factor(p[0, i, i])) for i in range(d + 1)))
        for part in set_partitions(list(range(1, d + 1))):
            if len(part) in (1, d):
                continue
            eqs = conditions(p, d, part)
            G = sp.groebner(eqs, s, t, order="lex", domain=sp.QQ)
            facs = sorted({str(f) for g in G.exprs for f, _ in sp.factor_list(g)[1]})
            # zero-dimensional ideals: exact real-root isolation (lex basis, last element univariate in t)
            gens = [sp.Poly(g, s, t) for g in G.exprs]
            zero_dim = G.exprs != [1] and any(g.degree(s) == 0 and g.degree(t) > 0 for g in gens) and any(
                sp.Poly(g, s, t).monoms()[0][0] > 0 and all(m[0] == sp.Poly(g, s, t).monoms()[0][0] or m[0] < sp.Poly(g, s, t).monoms()[0][0] for m in sp.Poly(g, s, t).monoms()) for g in gens)
            real_sols = []
            if G.exprs != [1] and G.is_zero_dimensional:
                zero_dim = True
                ut = [g for g in G.exprs if not g.has(s)]
                assert len(ut) == 1
                for f, _ in sp.factor_list(ut[0])[1]:
                    G2 = sp.groebner(list(G.exprs) + [f], s, t, order="lex", domain=sp.QQ)
                    if list(G2.exprs) == [1]:
                        continue
                    for t0 in sp.Poly(f, t).real_roots():
                        if t0.is_rational:
                            rest = [sp.Poly(g.subs(t, t0), s) for g in G2.exprs if g.has(s)]
                            h = rest[0]
                            for r_ in rest[1:]:
                                h = sp.gcd(h, r_)
                            for s0 in sp.Poly(h, s).real_roots():
                                real_sols.append((s0, t0, "exact"))
                        else:
                            # t0 irrational: s from the basis element of lowest degree in s, at 60 digits
                            g = min((g for g in G2.exprs if g.has(s)), key=lambda g: sp.Poly(g, s).degree())
                            t0n = t0.evalf(60)
                            for s0 in sp.Poly(g.subs(t, t0n), s).nroots(n=50):
                                if abs(sp.im(s0)) < 1e-40 and all(abs(gg.subs({s: sp.re(s0), t: t0n})) < 1e-30 for gg in G2.exprs):
                                    real_sols.append((sp.re(s0), t0n, "60-digit"))
            else:
                zero_dim = G.exprs == [1]
            real_ge1 = [(a_, b_) for a_, b_, _m in real_sols if bool(a_ >= 1) and bool(b_ >= 1)]
            sols = [((round(float(a_), 6), round(float(b_), 6)), m_) for a_, b_, m_ in real_sols]
            blocks = "|".join("{" + ",".join(map(str, B)) + "}" for B in sorted(part))
            lines.append(f"  {blocks:22s} #eqs={len(eqs):3d}  GB={[str(sp.factor(g)) for g in G.exprs]}")
            lines.append(f"  {'':22s} zero-dim={zero_dim}  real solutions={sols}  real s,t>=1: {real_ge1}")
        # direct integer cross-check, no Groebner: closing partitions for 1 <= s, t <= 30
        hits = {}
        for part in set_partitions(list(range(1, d + 1))):
            if len(part) in (1, d):
                continue
            eqs = conditions(p, d, part)
            f = sp.lambdify((s, t), eqs, "math")
            ok = [(a_, b_) for a_ in range(1, 31) for b_ in range(1, 31) if all(abs(e) < 0.5 for e in f(a_, b_))]
            hits["|".join("{" + ",".join(map(str, B)) + "}" for B in sorted(part))] = ok
        thick = {k_: [z for z in v_ if min(z) >= 2] for k_, v_ in hits.items()}
        lines.append(f"  integer scan 1<=s,t<=30 (exact integer polynomial values): closing (s,t) per partition = {hits}")
        lines.append(f"  thick (s,t>=2) closing pairs = {thick}")
    txt = "\n".join(lines) + "\n"
    open(out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main(*sys.argv[1:])
