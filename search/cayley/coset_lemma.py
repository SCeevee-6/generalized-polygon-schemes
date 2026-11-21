"""Checks for the coset lemma (cayley_note.md, 'A short argument').

Lemma. Let Gamma = Cay(G, S) be distance-regular on an abelian group G with integral eigenvalues
k = theta_0 > theta_1 > ... . Let p be a prime whose Sylow p-subgroup is Z_{p^e} + Z_p^r with e >= 2. Then
    (p-1)(k - theta_J)/p = |S \\ A| <= p^I_{11}
for some J >= 1 and I >= 1, where A is the index-p subgroup containing every element whose p-part has order <= p,
h generates the subgroup K of order p in p^{e-1}G, and h lies in S_I. In particular (p-1)(k - theta_1)/p <= max(a_1, c_2).

This script
1. checks the two group facts the lemma uses, by brute force, on Z_{p^e} + Z_p^r for small cases and on the groups
   that occur here: for s outside A, s + K lies in the atom of s, and lambda (kernel A) sums to -|atom|/(p-1) on it;
2. tests the lemma on the distance-regular Cayley graphs that abelian_cayley.py finds for some known arrays on
   such groups (identity and inequality, read off from the actual connection set);
3. evaluates the inequality for the five arrays that pass the Krein test, and lists for every abelian group of
   order v which argument excludes it.
Output: coset_lemma_output.txt
"""
import itertools, json, math, os, subprocess, sys
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import abelian_cayley as ac  # noqa: E402  (ortools import inside is optional)

OUT = []


def say(*a):
    OUT.append(' '.join(str(x) for x in a))


def parse(s):
    return ac.parse(s)


def params(arr):
    b, c = parse(arr)
    v, th, ms = ac.spectrum(b, c)
    k = b[0]; a1 = k - b[1] - c[0]; c2 = c[1]
    return dict(v=v, th=th, ms=ms, k=k, a1=a1, c2=c2, b=b, c=c)


def p11(P, i):
    """p^i_{11} for i = 1, 2 (0 for i >= 3)."""
    return {1: P['a1'], 2: P['c2']}.get(i, 0)


# ---------------------------------------------------------------- 1. group facts
def group_facts(p, e, r):
    ns = [p ** e] + [p] * r
    E = np.array(list(itertools.product(*[range(n) for n in ns])), dtype=np.int64)
    nsa = np.array(ns); n = p ** e
    h = np.zeros(len(ns), dtype=np.int64); h[0] = p ** (e - 1)
    K = [tuple((j * h) % nsa) for j in range(p)]
    # A: first coordinate divisible by p (index p, contains all elements with p-part of order <= p)
    inA = (E[:, 0] % p == 0)
    units = [u for u in range(1, n) if math.gcd(u, n) == 1]
    ok_coset = ok_sum = True
    for x in E[~inA]:
        atom = {tuple((u * x) % nsa) for u in units}
        if not all(tuple((x + np.array(k_)) % nsa) in atom for k_ in K): ok_coset = False
        # lambda(y) = zeta_p^{y_0 / p^{e-1}} ... kernel A: lambda(y) = exp(2 pi i y_0 / p)
        s = sum(np.exp(2j * np.pi * (y[0] % p) / p) for y in atom)
        if abs(s - (-len(atom) / (p - 1))) > 1e-9: ok_sum = False
    return ok_coset, ok_sum


# ---------------------------------------------------------------- 2. examples
def class_of(E, ns, Sidx, idx, target):
    """distance from 0 to target in Cay(G,S)."""
    nsa = np.array(ns); v = len(E); S = E[Sidx]
    dist = {tuple([0] * len(ns)): 0}; fr = [tuple([0] * len(ns))]; d = 0
    while fr:
        d += 1; nf = []
        for x in fr:
            for s in S:
                y = tuple((np.array(x) + s) % nsa)
                if y not in dist: dist[y] = d; nf.append(y)
        fr = nf
    return dist[tuple(target)]


def check_example(arr, ns, p):
    out = os.path.join(HERE, 'out_lemma', arr.split(';')[0].replace(',', '-') + '_' + 'x'.join(map(str, ns)) + '.json')
    if not os.path.exists(out):
        with open(out + '.err', 'w') as fe:
            subprocess.run([sys.executable, 'abelian_cayley.py', arr, ','.join(map(str, ns)), '--limit', '3',
                            '--tl', '120', '--workers', '4', '--out', out], cwd=HERE, stdout=fe, stderr=fe)
    d = json.load(open(out)); P = params(arr)
    if not d['drg_connection_sets']:
        say(f'  {arr} on {ns}: no graph found ({d["final_status"]})'); return
    E = ac.group(ns); nsa = np.array(ns); e = round(math.log(ns[0], p))
    idx = {tuple(x): i for i, x in enumerate(E)}
    for S in d['drg_connection_sets']:
        S = np.array(S); Sidx = [idx[tuple(x)] for x in S]
        beta = int((S[:, 0] % p != 0).sum())
        lam = sum(np.exp(2j * np.pi * (x[0] % p) / p) for x in S).real
        J = min(range(len(P['th'])), key=lambda j: abs(P['th'][j] - lam))
        h = [0] * len(ns); h[0] = p ** (e - 1)
        I = class_of(E, ns, Sidx, idx, h)
        lhs = Fraction((p - 1) * (P['k'] - P['th'][J]), p)
        say(f'  {arr} on {ns}: beta={beta}, (p-1)(k-theta_J)/p={lhs} (J={J}), I={I}, p^I_11={p11(P, I)},'
            f' identity {"ok" if lhs == beta else "FAIL"}, inequality {"ok" if beta <= p11(P, I) else "FAIL"}')


# ---------------------------------------------------------------- 3. the five arrays
def abelian_types(v):
    from sympy import factorint
    from sympy.utilities.iterables import partitions
    per = []
    for q, e in sorted(factorint(v).items()):
        per.append([(q, sorted([q ** a for a, cnt in part.items() for _ in range(cnt)], reverse=True))
                    for part in (dict(z) for z in partitions(e))])
    return list(itertools.product(*per))


def bercov_star(typ):
    """(*) of Bercov: H = A x B x C, A cyclic p-group of order p^a != 4, B of exponent < p^a, C of order prime to p."""
    for q, ns in typ:
        top = ns[0]
        if ns.count(top) == 1 and top != 4: return q
    return None


def valencies(P):
    b, c = P['b'], P['c']; ks = [1]
    for i in range(len(b)): ks.append(ks[-1] * b[i] // c[i])
    return ks[1:]


def main():
    os.makedirs(os.path.join(HERE, 'out_lemma'), exist_ok=True)
    say('1. group facts (s + K inside the atom of s; lambda sum over the atom = -|atom|/(p-1)), s outside A')
    for p, e, r in [(2, 2, 1), (2, 2, 2), (2, 3, 1), (3, 2, 1), (3, 2, 2), (3, 3, 1), (5, 2, 1), (3, 2, 4),
                    (5, 2, 2), (5, 3, 1), (2, 3, 0), (5, 3, 0)]:
        say(f'  p={p} e={e} r={r}:', *group_facts(p, e, r))
    say('2. lemma on graphs found by abelian_cayley.py')
    for arr, ns, p in [("3,2,1;1,2,3", [4, 2], 2), ("4,3,2,1;1,2,3,4", [4, 2, 2], 2), ("5,4;1,2", [4, 2, 2], 2),
                       ("6,3;1,2", [4, 2, 2], 2), ("6,4,2;1,2,3", [9, 3], 3), ("7,6,5;1,2,3", [4, 2, 2, 2, 2], 2)]:
        check_example(arr, ns, p)
    say('3. the five arrays that pass the Krein test')
    arrays = ["52,42,16;1,6,28", "56,42,20;1,6,28", "49,40,22;1,5,28", "54,42,11;1,3,44", "63,48,24;1,6,28"]
    for arr in arrays:
        P = params(arr); v = P['v']
        say(f'  {{{arr}}}: v={v}, theta={P["th"]}, a1={P["a1"]}, c2={P["c2"]}')
        from sympy import factorint
        bad = {}
        for q in factorint(v):
            lo = Fraction((q - 1) * (P['k'] - P['th'][1]), q)
            bad[q] = lo > max(P['a1'], P['c2'])
            feas = [(I, J) for J in (1, 2, 3) for I in (1, 2, 3)
                    if Fraction((q - 1) * (P['k'] - P['th'][J]), q) <= p11(P, I)]
            say(f'    p={q}: (p-1)(k-theta_1)/p = {lo} vs max(a1,c2) = {max(P["a1"], P["c2"])} ->',
                'inequality violated' if bad[q] else 'no contradiction', f'; feasible (I,J): {feas}')
        for typ in abelian_types(v):
            name = ' x '.join(f'Z{n}' for q, ns in typ for n in ns)
            why = None
            for q, ns in typ:
                big = [n for n in ns if n > q]
                if len(big) == 1 and bad[q]: why = f'lemma (p={q})'; break
            q0 = bercov_star(typ)
            hd = [d_ * (round(v ** (1 / d_)) - 1) for d_ in range(2, 13) if round(v ** (1 / d_)) ** d_ == v and round(v ** (1 / d_)) > 2]
            ber = (f'Bercov (*) with p={q0}; Hamming valencies {hd} vs class sizes {valencies(P)}: '
                   + ('excluded' if not set(hd) & set(valencies(P)) else 'NOT excluded')) if q0 else 'Bercov (*) fails'
            say(f'    {name:28s} {why or "not covered by the lemma":26s} | {ber}')
    open(os.path.join(HERE, 'coset_lemma_output.txt'), 'w').write('\n'.join(OUT) + '\n')
    print('\n'.join(l for l in OUT if 'FAIL' in l or 'False' in l) or 'all checks ok', len(OUT))


if __name__ == '__main__':
    main()
