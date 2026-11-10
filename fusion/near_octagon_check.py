"""Checks for the near-octagon fusion classification (see near_octagon_note.md).

1. The four parameter families satisfy the closure conditions identically (sympy).
2. Box scan: for all integers 2 <= s <= S, 0 <= t2 <= t3 < t <= T every partition is tested exactly (Fractions, no
   integrality filter); every closing (partition, point) must be one predicted by the families.
3. The three residual surfaces (one for each of {1,4}|{2,3}, {1,3}|{2,4}, {1,2,4}|{3}): for all integers
   2 <= s <= S2, 0 <= t3 <= T3 the residual polynomial is factored exactly as a polynomial in t (python-flint), and
   every integer root t > t3 is tested for an integer t2 with 0 <= t2 <= t3.

Usage: python near_octagon_check.py [S T S2 T3]     (defaults 10 30 200 2000; writes near_octagon_check_output.txt)
"""
import os, sys, json, time
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from near_polygon_fusions import scheme, set_partitions, conditions, int_array, p_numeric, closed_num

s, t, t2, t3, m, i = sp.symbols('s t t2 t3 m i')

FAMILIES = {
    'A': ({t2: m, t3: s * (m + 1), t: (m + 1) * (s ** 2 + 1) - 1}, ['[[1, 2], [3, 4]]']),
    'B': ({t2: 0, t3: 0, t: s}, ['[[1, 2, 3], [4]]']),
    'C': ({s: 2, t2: i, t3: 2 * i, t: 2 * i + 1}, ['[[1, 2, 4], [3]]', '[[1, 4], [2], [3]]']),
    'C1': ({s: 2, t2: 1, t3: 2, t: 3}, ['[[1, 2, 4], [3]]', '[[1, 4], [2], [3]]', '[[1, 4], [2, 3]]', '[[1, 3, 4], [2]]']),
    'D': ({s: 3, t2: 1, t3: 2, t: 3}, ['[[1, 3], [2, 4]]']),
}


def predicted(sv, a, b, tv):
    out = set()
    if b == sv * (a + 1) and tv == (a + 1) * (sv * sv + 1) - 1:
        out.add('[[1, 2], [3, 4]]')
    if a == 0 and b == 0 and tv == sv:
        out.add('[[1, 2, 3], [4]]')
    if sv == 2 and b == 2 * a and tv == 2 * a + 1:
        out |= {'[[1, 2, 4], [3]]', '[[1, 4], [2], [3]]'}
    if (sv, a, b, tv) == (2, 1, 2, 3):
        out |= {'[[1, 4], [2, 3]]', '[[1, 3, 4], [2]]'}
    if (sv, a, b, tv) == (3, 1, 2, 3):
        out.add('[[1, 3], [2, 4]]')
    return out


def parts():
    return [sorted(sorted(q) for q in p) for p in set_partitions([1, 2, 3, 4]) if len(p) not in (1, 4)]


def main(S=10, T=30, S2=200, T3=2000):
    out = []
    t0 = time.time()
    P = scheme(4)
    for name, (sub, pl) in FAMILIES.items():
        for part in parts():
            if str(part) in pl:
                ok = all(sp.cancel(e.subs(sub, simultaneous=True)) == 0 for e in conditions(P, part))
                out.append(f'family {name} {sub}: {part} closes identically: {ok}')
    out.append(f'[{time.time() - t0:.0f}s]')
    bad = []; n = 0; hits = {}
    for sv in range(2, S + 1):
        for tv in range(1, T + 1):
            for b in range(0, tv):
                for a in range(0, b + 1):
                    bb, cc = int_array(4, sv, (a, b), tv)
                    pf = p_numeric(bb, cc)
                    got = {str(part) for part in parts() if closed_num(pf, part)}
                    n += 1
                    if got != predicted(sv, a, b, tv):
                        bad.append(((sv, a, b, tv), sorted(got), sorted(predicted(sv, a, b, tv))))
                    for g in got:
                        hits[g] = hits.get(g, 0) + 1
    out.append(f'box scan 2 <= s <= {S}, 0 <= t2 <= t3 < t <= {T}: {n} points, closing (partition, point) pairs '
               f'{sum(hits.values())}, disagreements with the families: {len(bad)} {bad[:5]}')
    out.append('  by partition: ' + json.dumps(hits, sort_keys=True))
    out.append(f'[{time.time() - t0:.0f}s]')
    # residual surfaces
    import flint
    js = json.load(open(os.path.join(HERE, 'near_octagon_decompose.json')))
    for key in ('[[1, 4], [2, 3]]', '[[1, 3], [2, 4]]', '[[1, 2, 4], [3]]'):
        comp = [c for c in js[key] if c['residual'] and len(c['residual'][0]) > 300 and set(c['solved']) == {'t2'}]
        assert len(comp) == 1, (key, len(comp))
        R = sp.Poly(sp.sympify(comp[0]['residual'][0]), t)
        N, D = sp.fraction(sp.together(sp.sympify(comp[0]['solved']['t2'])))
        cf = [sp.lambdify((s, t3), sp.expand(c), 'math') for c in R.all_coeffs()[::-1]]
        Nf = sp.lambdify((s, t, t3), sp.expand(N), 'math'); Df = sp.lambdify((s, t, t3), sp.expand(D), 'math')
        part = json.loads(key)
        found = []; nroots = 0; outside = 0
        for sv in range(2, S2 + 1):
            for b in range(0, T3 + 1):
                co = [int(f(sv, b)) for f in cf]
                if not any(co):
                    found.append(('identically zero', sv, b)); continue
                pol = flint.fmpz_poly(co)
                for fac, _ in pol.factor()[1]:
                    if fac.degree() == 1:
                        c0, c1 = int(fac[0]), int(fac[1])
                        if c0 % c1:
                            continue
                        tv = -c0 // c1
                        if tv <= b:
                            continue
                        nroots += 1
                        d = Df(sv, tv, b); nn = Nf(sv, tv, b)
                        cand = range(0, b + 1) if d == 0 else ([nn // d] if nn % d == 0 else [])
                        adm = [a for a in cand if 0 <= a <= b]
                        if not adm:
                            outside += 1
                        for a in adm:
                            bb, cc = int_array(4, sv, (a, b), tv)
                            if closed_num(p_numeric(bb, cc), part) and str(part) not in predicted(sv, a, b, tv):
                                found.append((sv, a, b, tv))
        out.append(f'residual surface of {key}: 2 <= s <= {S2}, 0 <= t3 <= {T3}: integer roots t > t3: {nroots}, '
                   f'of which t2 not an integer in [0, t3]: {outside}; new admissible points: {found}')
        out.append(f'[{time.time() - t0:.0f}s]')
    txt = '\n'.join(out)
    open(os.path.join(HERE, 'near_octagon_check_output.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main(*(int(a) for a in sys.argv[1:5]))
