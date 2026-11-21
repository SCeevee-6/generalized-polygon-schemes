"""The three residual surfaces of near_octagon_decompose.py have no admissible integer point.

Regular near octagon (s, t2, t3, t): c_i = t_i + 1, b_i = s(t - t_i), t_1 = 0, t_4 = t.
Admissible means integers with s >= 2 and 0 <= t2 <= t3 < t.

Part 1 is a direct proof for {1,2,4}|{3} (Gamma_3 strongly regular), written so it can be
followed by hand. Part 2 is a general region argument, run on all three surfaces:

  On a surface t2 = N/D, R(s, t3, t) = 0, the number of admissible real roots t of R(s, t3, .)
  is locally constant in (s, t3) away from the zeros of
      lc_t(R), disc_t(R), R(s, t3, t3), Res_t(R, N), Res_t(R, N - t3 D), Res_t(R, D).
  If every irreducible factor of these has constant sign on a connected region (checked by
  substituting a polynomial or rational parametrisation of the region and looking at the signs
  of all coefficients), one sample point gives the count on the whole region. Rows and columns
  with s or t3 fixed are handled the same way in one variable, plus an exact check of the
  integer points below the largest real root.

Usage: python near_octagon_residual.py [output file]
"""
import json
import os
import sys
import time

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 60

s, t2, t3, t = sp.symbols("s t2 t3 t")
n, u = sp.symbols("n u")
zz, ww, lam, v = sp.symbols("zz ww lam v", nonnegative=True)

LOG = []


def say(*a):
    line = " ".join(str(x) for x in a)
    LOG.append(line)
    print(line, flush=True)


# ---------------------------------------------------------------- intersection numbers
def intersection_matrices():
    c = [0, 1, t2 + 1, t3 + 1, t + 1]
    b = [s * (t + 1), s * t, s * (t - t2), s * (t - t3), 0]
    k = b[0]
    a = [k - b[i] - c[i] for i in range(5)]
    L1 = sp.zeros(5, 5)
    for i in range(5):
        L1[i, i] = a[i]
        if i > 0:
            L1[i - 1, i] = c[i]
        if i < 4:
            L1[i + 1, i] = b[i]
    I = sp.eye(5)
    L2 = sp.simplify((L1 * L1 - a[1] * L1 - k * I) / c[2])
    L3 = sp.simplify((L1 * L2 - b[1] * L1 - a[2] * L2) / c[3])
    return L1, L2, L3  # (L_i)[j, h] = p^h_{ij}


# ---------------------------------------------------------------- part 1
def part1():
    say("== Part 1: {1,2,4}|{3}, direct proof")
    _, _, L3 = intersection_matrices()
    p = [sp.together(L3[3, h]) for h in range(5)]
    E1 = sp.factor(sp.numer(sp.together(p[1] - p[2])))
    E2 = sp.factor(sp.numer(sp.together(p[4] - p[2])))
    say("p^1_33 - p^2_33 numerator:", E1)
    # E1 = s (t - t2) * (linear in t2); solve the bracket for t2
    # the factors s, t - t2 and t2 + 1 are nonzero on the region; keep the one with t3 in it
    br = [f for f, m in sp.factor_list(E1)[1] if f.has(t3)][0]
    T2 = sp.solve(br, t2)[0]
    E2b = sp.factor(sp.numer(sp.together(E2.subs(t2, T2))))
    say("p^4_33 - p^2_33 on that branch:", E2b)
    fl = sp.factor_list(E2b)[1]
    R = [f for f, m in fl if sp.Poly(f, t).degree() == 3][0]
    if sp.Poly(R, t).LC().subs(s, 2) < 0:
        R = -R
    ufac = sp.expand(s * (t3 + 1) - t - t3 - 1)
    assert any(sp.expand(f - ufac) == 0 or sp.expand(f + ufac) == 0 for f, m in fl)
    t2p1 = sp.expand(sp.numer(sp.together(T2 + 1)))
    rest = [f for f, m in fl if f.has(t) and sp.Poly(f, t).degree() < 3
            and sp.expand(f - ufac) != 0 and sp.expand(f + ufac) != 0]
    for f in rest:
        assert sp.expand(f - t2p1) == 0 or sp.expand(f + t2p1) == 0
    say("other factors: t3 + 1, and", len(rest), "factor(s) equal to the numerator of t2 + 1 (nonzero)")
    say("so either t = (s-1)(t3+1) (family C) or R = 0 with R cubic in t")
    # family C branch: t2 = ((s-1)^2 (t3+1) - 1)/s, t2 <= t3 forces s^2 - 3s + 1 < 0, s = 2
    T2C = sp.factor(T2.subs(t, (s - 1) * (t3 + 1)))
    say("t2 on t = (s-1)(t3+1):", T2C, "; t3 - t2 =", sp.factor(t3 - T2C))

    # coordinates n = t3 + 1, M = (s-1) n + 1, u = (s-1) n - t, so t + 1 = M - u
    M = (s - 1) * n + 1
    F = sp.expand(T2 * s * (t3 + 1))
    Fu = sp.expand(F.subs(t3, n - 1).subs(t, (s - 1) * n - u))
    check = sp.expand(Fu - (n * ((s - 1) ** 2 * n - 1) - s * u * (M - u)))
    say("s n t2 = n((s-1)^2 n - 1) - s u (t+1):", check == 0)
    pL = n * ((s ** 2 - 3 * s + 1) * n + s - 1) / s
    say("t2 <= t3  <=>  u (t+1) >= pL with pL = n((s^2-3s+1)n + s - 1)/s")
    say("for s >= 3, pL > 0, so u > 0 and u >= pL/M, hence t <= M - 1 - pL/M < T := M - pL/M")

    Rn = sp.Poly(sp.expand(R.subs(t3, n - 1)), t)
    lc = sp.factor(Rn.LC())
    say("leading coefficient in t:", lc)
    Rt3 = sp.factor(Rn.as_expr().subs(t, n - 1))
    say("(b) R(t3) =", Rt3)
    Na = sp.expand(8 * Rn.as_expr().subs(t, (n - 1) / 2))
    Tb = M - pL / M
    Nc = sp.expand(-sp.numer(sp.together(Rn.as_expr().subs(t, Tb))))
    den = sp.factor(sp.denom(sp.together(Rn.as_expr().subs(t, Tb))))
    say("R(T) = -Nc /", den)

    def shift_pos(P, sub, gens):
        Q = sp.Poly(sp.expand(P.subs(sub)), *gens)
        return min(Q.coeffs()) > 0

    say("(a) 8 R(t3/2) has positive coefficients for s = 3+zz, n = 2+ww:",
        shift_pos(Na, {s: 3 + zz, n: 2 + ww}, (zz, ww)))
    say("(c) Nc has positive coefficients for s = 5+zz, n = 5+ww:",
        shift_pos(Nc, {s: 5 + zz, n: 5 + ww}, (zz, ww)))
    say("    Nc(4, n) > 0 for n >= 7:", shift_pos(Nc.subs(s, 4), {n: 7 + ww}, (ww,)))
    say("    Nc(3, n) > 0 for n >= 122:", shift_pos(Nc.subs(s, 3), {n: 122 + ww}, (ww,)))
    say("Root count: R(t3/2) > 0, R(t3) < 0, R(T) < 0 and lc > 0 leave no root in [t3, T]:")
    say("  two roots in (t3, T) would force a third beyond T, and then R < 0 on (-inf, t3].")
    say("(b) needs n >= 2 for s >= 4 and n >= 3 for s = 3.")

    # s = 2: R(2, t3+1+x, t3) has positive coefficients
    x = sp.symbols("x", nonnegative=True)
    R2 = sp.Poly(sp.expand(R.subs({s: 2, t: t3 + 1 + x})), t3, x)
    say("s = 2: R(2, t3+1+x, t3) has positive coefficients:", min(R2.coeffs()) > 0)

    # t3 <= 3: t2 in {0..t3}, eliminate t
    found = []
    for T3 in range(0, 4):
        for TT2 in range(0, T3 + 1):
            e = sp.expand(F.subs(t3, T3) - s * (T3 + 1) * TT2)
            res = sp.resultant(e, sp.expand(R.subs(t3, T3)), t)
            for S in sp.roots(sp.Poly(res, s), filter="Z"):
                if S >= 2:
                    for r in sp.roots(sp.Poly(R.subs({s: S, t3: T3}), t), filter="Z"):
                        if r > T3 and e.subs({s: S, t: r}) == 0:
                            found.append((S, TT2, T3, r))
    say("t3 <= 3, all s: admissible integer points:", found)

    # finite part: s = 3, t3 <= 120 and s = 4, t3 <= 5
    found = []
    for S, top in [(3, 120), (4, 5)]:
        for T3 in range(0, top + 1):
            for r in sp.roots(sp.Poly(R.subs({s: S, t3: T3}), t), filter="Z"):
                if r > T3:
                    val = sp.Rational(T2.subs({s: S, t: r, t3: T3}))
                    if 0 <= val <= T3:
                        found.append((S, val, T3, r))
    say("s = 3, t3 <= 120 and s = 4, t3 <= 5: admissible integer points:", found)
    say("Part 1 conclusion: Gamma_3 is strongly regular only for (2, i, 2i, 2i+1).")
    return R, T2


# ---------------------------------------------------------------- part 2
def boundary(R, T2):
    N = sp.numer(sp.together(T2))
    D = sp.denom(sp.together(T2))
    P = sp.Poly(R, t)
    polys = {"lc": P.LC(), "R(t3)": R.subs(t, t3), "disc": sp.discriminant(P),
             "Res(R,N)": sp.resultant(R, N, t),
             "Res(R,N-t3D)": sp.resultant(R, sp.expand(N - t3 * D), t)}
    if sp.Poly(D, t).degree() > 0:
        polys["Res(R,D)"] = sp.resultant(R, D, t)
    else:
        polys["D"] = D
    facs = []
    for key, val in polys.items():
        for f, m in sp.factor_list(sp.expand(val))[1]:
            if f.free_symbols and not any(sp.expand(f - g) == 0 for g, _ in facs):
                facs.append((f, key))
    return N, D, facs


def definite(f, sub, gens):
    num = sp.numer(sp.together(f.subs(sub)))
    cs = sp.Poly(sp.expand(num), *gens).coeffs()
    return all(c > 0 for c in cs) or all(c < 0 for c in cs)


def count_admissible(R, N, D, S, T3):
    cnt = 0
    for r in sp.real_roots(sp.Poly(R.subs({s: S, t3: T3}), t)):
        rv = sp.N(r, 50)
        if rv <= T3:
            continue
        Dv = sp.N(D.subs({s: S, t3: T3, t: r}), 50)
        if Dv == 0:
            return None
        val = sp.N(N.subs({s: S, t3: T3, t: r}), 50) / Dv
        if 0 <= val <= T3:
            cnt += 1
    return cnt


def integer_points(R, N, D, S, T3):
    out = []
    for r in sp.roots(sp.Poly(R.subs({s: S, t3: T3}), t), filter="Z"):
        if r > T3:
            Dv = D.subs({s: S, t3: T3, t: r})
            Nv = N.subs({s: S, t3: T3, t: r})
            if Dv == 0:
                out.append((S, f"{Nv}/0", T3, r))
            else:
                val = sp.Rational(Nv, Dv)
                if 0 <= val <= T3:
                    out.append((S, val, T3, r))
    return out


def line(R, N, D, facs, fixed, val, lo):
    other = t3 if fixed == s else s
    rho = lo
    for f, _ in facs:
        g = sp.Poly(f.subs(fixed, val), other)
        if g.is_zero:
            return f"factor vanishes identically on {fixed} = {val}"
        if g.degree() > 0:
            rr = [float(r) for r in sp.real_roots(g)]
            if rr:
                rho = max(rho, max(rr))
    top = int(rho) + 1
    pt = (val, top + 1) if fixed == s else (top + 1, val)
    cnt = count_admissible(R, N, D, *pt)
    pts = []
    for x in range(lo, top + 1):
        pts += integer_points(R, N, D, *((val, x) if fixed == s else (x, val)))
    return f"largest boundary root {rho:.3f}; real admissible roots at {pt}: {cnt}; " \
           f"integer points with {other} <= {top}: {pts}"


def part2(surfaces):
    say("== Part 2: region argument")
    for name, (R, T2, regions, srows, trows) in surfaces.items():
        t0 = time.time()
        N, D, facs = boundary(R, T2)
        say(f"-- {name}: {len(facs)} boundary factors, degrees "
            f"{sorted(sp.Poly(f, s, t3).total_degree() for f, _ in facs)}")
        for rname, sub, gens, sample, skip in regions:
            bad = [k for f, k in facs if not (skip and skip(f)) and not definite(f, sub, gens)]
            say(f"   region {rname}: factors not sign-definite: {bad}; "
                f"real admissible roots at {sample}: {count_admissible(R, N, D, *sample)}")
        for S in srows:
            say(f"   s = {S}: " + line(R, N, D, facs, s, S, 0))
        for T3, lo in trows:
            say(f"   t3 = {T3}, s >= {lo}: " + line(R, N, D, facs, t3, T3, lo))
        say(f"   ({time.time() - t0:.0f} s)")


def main(out):
    R124, T2_124 = part1()
    dec = json.load(open(os.path.join(HERE, "near_octagon_decompose.json")))
    loc = dict(s=s, t=t, t2=t2, t3=t3)

    def residual(key, j):
        b = dec[key][j]
        return (sp.sympify(b["residual"][0], locals=loc), sp.sympify(b["solved"]["t2"], locals=loc))

    R14, T2_14 = residual("[[1, 4], [2, 3]]", 3)
    R13, T2_13 = residual("[[1, 3], [2, 4]]", 10)

    # {1,3}|{2,4}: the factor C = s(s-1)t3 - s - (t3+1)^2 (common zero of R, N and D) has its
    # large root in (s^2-s-3, s^2-s-2) for s >= 4, so split s >= 6 there. The degree-20
    # discriminant factor is tangent to zero at infinity along t3 = s, so B is split at t3 = s.
    S6 = 6 + zz
    C = sp.expand(s ** 2 * t3 - s * t3 - s - t3 ** 2 - 2 * t3 - 1)
    regA = ("A: s >= 6, t3 >= s^2-s-2", {s: S6, t3: S6 ** 2 - S6 - 2 + ww}, (zz, ww), (6, 28), None)
    regB1 = ("B1: s >= 6, s <= t3 <= s^2-s-3",
             {s: S6, t3: S6 + (S6 ** 2 - 2 * S6 - 3) * lam / (1 + lam)}, (zz, lam), (6, 10), None)
    regB2 = ("B2: s >= 6, 1 <= t3 <= s", {s: S6, t3: S6 - (S6 - 1) * lam / (1 + lam)}, (zz, lam),
             (6, 3), None)
    say("C on A at t3 = s^2-s-2+ww:", sp.factor(sp.expand(C.subs(t3, s ** 2 - s - 2 + ww))))
    say("C at t3 = s^2-s-3:", sp.factor(C.subs(t3, s ** 2 - s - 3)), "(> 0 for s >= 4; C < 0 on A, so C has no zero in A, B1, B2)")

    surfaces = {
        "{1,2,4}|{3}": (R124, T2_124,
                        [("s >= 4, t3 >= 2", {s: 4 + zz, t3: 2 + ww}, (zz, ww), (5, 3), None)],
                        [2, 3], [(0, 4), (1, 4)]),
        "{1,4}|{2,3}": (R14, T2_14,
                        [("s >= 4, t3 >= 1", {s: 4 + zz, t3: 1 + ww}, (zz, ww), (5, 2), None)],
                        [2, 3], [(0, 4)]),
        "{1,3}|{2,4}": (R13, T2_13, [regA, regB1, regB2], [2, 3, 4, 5], [(0, 6)]),
    }
    part2(surfaces)
    say("Points listed with N/0 lie on the base locus N = D = 0 of the branch; they are the")
    say("branch s = (t3+1)^2/t of near_octagon_decompose.py, whose only admissible point is H(4,4).")
    open(out, "w", encoding="utf8", newline="\n").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "near_octagon_residual_output.txt"))
