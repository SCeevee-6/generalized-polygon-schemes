"""Non-backtracking moment framework of nb_moment_sdp.md."""

from fractions import Fraction as Q

EIGENVALUES = (12, 7, 3, -1, -3)
SHELLS = (1, 12, 96, 768, 2048)
PER_LINE = 1024  # opposite pairs per line


def cosines(e):
    c = [Q(1), Q(e, 12)]
    for _ in range(3):
        c.append(((e - 3)*c[-1] - c[-2]) / 8)
    return c


def walk_expansion(jmax):
    """Rows c[j][k] with J^j = sum_k c[j][k] NB_k in a cubic graph."""
    nb = [[1], [0, 1], [-3, 0, 1]]  # NB_k as coefficient lists in J
    for k in range(2, jmax):
        shifted = [0] + nb[k]
        prev = nb[k-1] + [0]*(len(shifted) - len(nb[k-1]))
        nb.append([a - 2*b for a, b in zip(shifted, prev)])
    rows = []
    for j in range(jmax + 1):
        rem, row = [0]*j + [1], [0]*(jmax + 1)
        for k in range(j, -1, -1):
            coef = rem[k] if k < len(rem) else 0
            row[k] = coef
            for i, v in enumerate(nb[k]):
                rem[i] -= coef*v
        assert not any(rem)
        rows.append(row)
    return rows


def q(dist, e):
    return sum(a*v for a, v in zip(cosines(e), dist))


def known_distributions(p):
    """Q_k(d) per line for k <= 4 in a grid-free octagon; p = 16-shape fraction."""
    return [[2048, 0, 0, 0, 2048],
            [0, 0, 12288, 0, 0],
            [6144, 0, 0, 0, 18432],
            [0, 0, 36864, 0, 12288],
            [12288, 0, 8192*p, 24576 - 8192*p, 61440]]


def repo_moments(e, p):
    a2, a3, a4 = cosines(e)[2:]
    g = 1024*(3 + 6*a3) + 2048*p*(a2 - a3)
    return [2048*(1 + a4), 12288*a2, 12288*(1 + 2*a4), 12288*(8*a2 + a4),
            36*2048*(1 + a4) + 24*6144*a4 + 4*g]


def check_exact():
    c = walk_expansion(8)
    assert c[4][:5] == [15, 0, 7, 0, 1] and c[8] == [543, 0, 319, 0, 95, 0, 15, 0, 1]
    for p in (Q(0), Q(1, 3), Q(1)):
        dist = known_distributions(p)
        for k, row in enumerate(dist):
            assert min(row) >= 0 and sum(row) == 4*3*2**(k-1)*PER_LINE if k else sum(row) == 4*PER_LINE
        for e in EIGENVALUES:
            mine = [sum(c[j][k]*q(dist[k], e) for k in range(j + 1)) for j in range(5)]
            assert mine == repo_moments(e, p)
    # Walks of length 4 along indexed systems: 12 per pair, each ending on a pair
    # that shares exactly one line. They account for all of Q_4(0).
    assert known_distributions(0)[4][0] == 12*PER_LINE
    print("walk expansion and Q_0..Q_4 reproduce the projected moments m_0..m_4")


def sdp_range():
    import cvxpy as cp
    c = walk_expansion(8)
    unknown = {k: cp.Variable(5, nonneg=True) for k in (5, 6, 7, 8)}
    p, t = cp.Variable(), cp.Variable()
    cons = [p >= 0, p <= 1, t >= 0, t <= 1]
    for k, var in unknown.items():
        cons.append(cp.sum(var) == 12*2**(k-1)*PER_LINE)
    cons += [unknown[6][0] >= 12*PER_LINE, unknown[6][4] >= 24*PER_LINE,
             unknown[8][0] >= 12*PER_LINE*(2 - t), unknown[8][4] >= 24*PER_LINE,
             unknown[5][2] >= 36*PER_LINE, unknown[7][2] >= 36*PER_LINE]
    for e in EIGENVALUES[1:]:
        a = [float(x) for x in cosines(e)]
        dist = [[float(v) for v in row] for row in known_distributions(0)[:4]]

        def qk(k):
            if k < 4:
                return sum(a[d]*dist[k][d] for d in range(5))
            if k == 4:
                row = [12288, 0, 8192*p, 24576 - 8192*p, 61440]
                return sum(a[d]*row[d] for d in range(5))
            return a @ unknown[k]
        qs = [qk(k) for k in range(9)]
        m = [sum(c[j][k]*qs[k] for k in range(j + 1)) / (3.0**j*2048) for j in range(9)]
        mats = ([[m[i+j] for j in range(5)] for i in range(5)],
                [[m[i+j] + m[i+j+1] for j in range(4)] for i in range(4)],
                [[m[i+j] - m[i+j+1] for j in range(4)] for i in range(4)],
                [[m[i+j] - m[i+j+2] for j in range(4)] for i in range(4)])
        for rows in mats:
            x = cp.bmat(rows)
            cons.append((x + x.T)/2 >> 0)
    out = []
    for objective in (cp.Minimize(t), cp.Maximize(t)):
        problem = cp.Problem(objective, cons)
        problem.solve(solver="CLARABEL")
        out.append((problem.status, round(float(t.value), 6)))
    return out


if __name__ == "__main__":
    check_exact()
    try:
        print("SDP range of the long-cycle fraction:", sdp_range())
    except ImportError:
        print("cvxpy not installed; skipped the SDP")
