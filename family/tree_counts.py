"""Exact triple counts that hold in EVERY generalized octagon of order (s, t): the free-tree count.

Setting as in triple/triple_zero_proofs.py: u, v, w points, (U, V, W) = (d(u,v), d(u,w), d(v,w)); incidence graph
of girth 16, so balls of radius 7 are trees and non-backtracking paths of length <= 8 are geodesics.

Base configuration B: the union of the geodesics (in the incidence graph) between the non-opposite pairs. If at
most one pair is opposite and B is a tree, B is a tripod with legs (l_u, l_v, l_w) (all of one parity, the leaf
distances l_a + l_b equal 2 d(a,b) for non-opposite pairs and are >= 8 for the opposite pair). The legs are fixed
by the distances except for one opposite pair, where the shared part g = l_c can be 0 .. V + W - 4 (c the third
point); these are the base shapes of a type.

Lemma. Let i, j, k <= 3 with  max_{a != b} (2 i_a + 2 i_b) + diam(B) < 16. Then every point x with
d(x,u) = i, d(x,v) = j, d(x,w) = k hangs off B: the union of B and the three geodesics from x is a tree (a circuit
in it would have length <= diam(T_x) + diam(B) < 16). Conversely, every way of hanging a path off B in the free
(s+1, t+1)-biregular tree that ends at a point with tree distances (2i, 2j, 2k) is realised by exactly one such x
(the path exists in the octagon, its distances are true distances because a shortcut would close a circuit of
length < 12, and x determines its attachment because geodesics are unique). So
  [ijk] = sum over nodes p of B and tail lengths l >= 0 with l + d_B(p, a) = 2 i_a for a = u, v, w and the end a
          point, of N(p, l), where N(p, 0) = 1 and N(p, l) = (deg p - deg_B p) * prod_{m=1}^{l-1} (deg(m-th node) - 1),
with deg = t + 1 at points and s + 1 at lines. This is a polynomial in s and t and depends only on the base shape.
Usage: imported by family/triple_family.py; run alone it checks the counts against the Ree-Tits octagon O(2).
"""
import itertools

def shapes(U, V, W):
    """tripod legs (l_u, l_v, l_w) for the base configurations of type (U, V, W) with B a tree; None if two pairs
    are opposite (w isolated). For all-non-opposite types with U + V + W >= 8 a circuit shape also exists; only the
    tripod is returned here."""
    D = {('u', 'v'): U, ('u', 'w'): V, ('v', 'w'): W}
    opp = [p for p, d in D.items() if d == 4]
    if len(opp) >= 2:
        return None
    if not opp:
        lu = U + V - W; lv = U + W - V; lw = V + W - U          # half-sums times 2 / 2: legs in incidence steps
        legs = (lu, lv, lw)
        if min(legs) < 0:
            return []
        return [legs]
    (a, b), = opp
    c = ({'u', 'v', 'w'} - {a, b}).pop()
    da, db = 2 * D[tuple(sorted((a, c)))], 2 * D[tuple(sorted((b, c)))]
    out = []
    for g in range(0, min(da, db) + 1):
        la, lb, lc = da - g, db - g, g
        if la + lb >= 8 and (la - lb) % 2 == 0:
            legs = {a: la, b: lb, c: lc}
            out.append((legs['u'], legs['v'], legs['w']))
    return out

def tree_count(legs, i, j, k, s, t):
    """number of points at tree distances (2i, 2j, 2k) from the leaves u, v, w of a tripod with the given legs,
    in the free (s+1, t+1)-biregular tree; leaves are points."""
    lu, lv, lw = legs
    # nodes of the tripod: ('c',) the centre (depth from u: lu) and ('a', r) at distance r from the centre on leg a
    nodes = [('c', 0)] + [(a, r) for a, la in zip('uvw', legs) for r in range(1, la + 1)]
    L = dict(zip('uvw', legs))
    def dist(n, a):
        leg, r = n
        if leg == 'c':
            return L[a]
        return L[a] - r if leg == a else r + L[a]
    def is_point(n):
        return dist(n, 'u') % 2 == 0
    def degB(n):
        leg, r = n
        if leg == 'c':
            return sum(1 for l_ in legs if l_ > 0)
        return 1 if r == L[leg] else 2
    target = {'u': 2 * i, 'v': 2 * j, 'w': 2 * k}
    total = 0
    for n in nodes:
        dn = {a: dist(n, a) for a in 'uvw'}
        ell = target['u'] - dn['u']
        if ell < 0 or any(target[a] - dn[a] != ell for a in 'uvw'):
            continue
        pt = is_point(n)
        if ell == 0:
            total += 1 if pt else 0
            continue
        if (ell % 2 == 0) != pt:
            continue
        deg = (t + 1) if pt else (s + 1)
        cnt = deg - degB(n)
        cur_pt = not pt
        for m in range(1, ell):
            cnt *= ((t + 1) if cur_pt else (s + 1)) - 1
            cur_pt = not cur_pt
        total += cnt
    return total

def admissible(legs, i, j, k):
    lu, lv, lw = legs
    diamB = max(lu + lv, lu + lw, lv + lw)
    ia = (i, j, k)
    return max(2 * ia[0] + 2 * ia[1], 2 * ia[0] + 2 * ia[2], 2 * ia[1] + 2 * ia[2]) + diamB < 16

def constraints(U, V, W, s, t):
    """{shape legs: {(i,j,k): exact count}} for all admissible (i,j,k) with i, j, k <= 3."""
    sh = shapes(U, V, W)
    if not sh:
        return {}
    out = {}
    for legs in sh:
        out[legs] = {(i, j, k): tree_count(legs, i, j, k, s, t)
                     for i, j, k in itertools.product(range(4), repeat=3) if admissible(legs, i, j, k)}
    return out
