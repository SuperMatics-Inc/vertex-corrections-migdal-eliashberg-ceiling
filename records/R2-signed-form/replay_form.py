"""Exact finite algebra checks; no continuum or Matsubara approximation claim."""

from fractions import Fraction as F
from itertools import product


GRID = 5
W = {(0, -1): F(1, 5), (0, 0): F(11, 7), (0, 1): F(1, 5),
     (1, -1): F(2, 9), (1, 0): F(1, 2), (1, 1): F(2, 9),
     (4, -1): F(2, 9), (4, 0): F(1, 2), (4, 1): F(2, 9)}
POINTS = tuple(product(range(GRID), range(-4, 4)))
ZERO = (F(0), F(0))


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def scale(z, c):
    return z[0] * c, z[1] * c


def mm(a, b):
    return tuple(tuple(add(mul(a[i][0], b[0][j]), mul(a[i][1], b[1][j]))
                       for j in range(2)) for i in range(2))


def sub(p, s):
    return (p[0] - s[0]) % GRID, p[1] - s[1]


def lines(p):
    k, n = p
    d = min(k, GRID - k)
    nf = max(n, -n - 1)
    sign = 1 if n >= 0 else -1
    u = sign * (2 * nf + 1 + F(d + 1, nf + 3))
    e = F((d + 1) * (nf + 1), nf + 2) - F(3, 5)
    y = sign * F((nf + 1) * (d + 2), 3) if nf <= 1 else F(0)
    return -e / (e * e + u * u), u / (e * e + u * u), y


def signed(mutation=None, static=False):
    total = F(0)
    for p, (s, v_s), (t, v_t) in product(POINTS, W.items(), W.items()):
        if static and (s[1] or t[1]):
            continue
        a, c = sub(p, s), sub(p, t)
        corners = (p, a, sub(a, t), c)
        r, q, y = zip(*(lines(w) for w in corners))
        alt = y[0] - y[1] + y[2] - y[3]
        term = -q[0] * q[1] * q[2] * q[3] * alt**2 / (2 if mutation == 'square' else 4)
        if mutation != 'drop-real':
            for j in range(1, 4):
                i, k = [v for v in range(1, 4) if v != j]
                term += q[0] * q[j] * r[i] * r[k] * (y[0] - y[j])**2 / 2
        total += v_s * v_t * term
    return total


def nambu(mutation=None, static=False):
    # Independent test: matrix multiply the literal vertex and extract the
    # anomalous derivative by exact odd interpolation at ±1 and ±2.
    sigma3 = (((F(1), F(0)), ZERO), (ZERO, (F(-1), F(0))))

    def crossed(points, epsilon):
        out = sigma3
        for i, p in enumerate(points):
            real, q, y = lines(p)
            sign = -1 if mutation == 'middle-sign' and i == 1 else 1
            h = -sign * epsilon * q * y, F(0)
            out = mm(mm(out, (((real, -q), h), (h, (-real, -q)))), sigma3)
        return out

    total = F(0)
    for p in POINTS:
        _, qp, yp = lines(p)
        if yp == 0:
            continue
        for (s, vs), (t, vt) in product(W.items(), repeat=2):
            if static and (s[1] or t[1]):
                continue
            a, c = sub(p, s), sub(p, t)
            points = [a, sub(a, t), c]
            normal = -crossed(points, 0)[0][0][1]

            def anomalous(eps):
                return crossed(points, eps)[0][1][0]

            diff1 = anomalous(1) - anomalous(-1)
            diff2 = anomalous(2) - anomalous(-2)
            linear = (8 * diff1 - diff2) / 12
            total += vs * vt * qp * yp * ((0 if mutation == 'omit-normal' else normal * yp) - linear)
    return total


def folded_momentum_mean():
    # With |Q| already integrated in b, omit the same global factor 2νπN0
    # from the folded form. Compare to the unfolded form with half weight.
    line = range(-3, 3)
    k = (0, 1, 2)
    coupling = (F(2), F(3, 4), F(1, 3), F(1, 5), F(1, 9), F(1, 13))
    b = (F(4, 7), F(3, 8), F(7, 10))
    prob = (F(1, 2), F(1, 3), F(1, 6))

    def idx(n):
        return max(n, -1 - n)

    def y(n, x):
        return (1 if n >= 0 else -1) * F((x + 2) * (idx(n) + x + 1), x + 1)

    def q(n, x):
        return (1 if n >= 0 else -1) * b[idx(n)] * prob[x]

    free = sum((2 * n + 1) * q(n, x) * y(n, x)**2 for n, x in product(line, k)) / 2
    me = sum(coupling[abs(n - m)] * q(n, x) * q(m, v) * (y(n, x) - y(m, v))**2
             for n, m, x, v in product(line, line, k, k)) / 4
    mean = [sum(prob[x] * y(n, x) for x in k) for n in range(3)]
    second = [sum(prob[x] * y(n, x)**2 for x in k) for n in range(3)]
    diag = [sum((coupling[abs(n - m)] - coupling[n + m + 1]) * b[m] for m in range(3))
            for n in range(3)]
    penalty = sum(b[n] * diag[n] * second[n] for n in range(3))
    pair = sum(b[n] * b[m] * (coupling[abs(n - m)] + coupling[n + m + 1]) * mean[n] * mean[m]
               for n, m in product(range(3), repeat=2))
    assert free + me == free + penalty - pair
    mutant = sum(b[n] * diag[n] * mean[n]**2 for n in range(3))
    assert free + me != free + mutant - pair


if __name__ == '__main__':
    value = signed()
    assert value == nambu()
    assert value != signed('square')
    assert value != signed('drop-real')
    assert value != nambu('middle-sign')
    assert value != nambu('omit-normal')
    assert signed(static=True) == nambu(static=True) != 0
    folded_momentum_mean()
    print('Exact Nambu interpolation = signed form (full and static); folded ME mean/variance; six mutants rejected.')
