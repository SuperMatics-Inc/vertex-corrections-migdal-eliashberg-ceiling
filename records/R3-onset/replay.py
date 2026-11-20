"""Exact finite checks for the signed-density and weak-coupling onset proofs."""

from fractions import Fraction as F
from itertools import product
from math import comb, factorial


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def closed_density(k, a, b, c, prefactor=48, include_triangle=True):
    length = max(0, min(a + k, b + c) - max(abs(a - k), abs(b - c)))
    fraction = length / (2 * k) if include_triangle else 1
    return fraction / (prefactor * a * b * c)


def sliced_density(k, a, b, c):
    # W = 1. Independently integrate the joint law in distance d, using
    # the spherical cosine density and the energy derivative pointwise.
    edges = sorted({abs(a - k), a + k, abs(b - c), b + c})
    integral = F(0)
    for left, right in zip(edges, edges[1:]):
        d = (left + right) / 2
        cosine = (b * b - c * c - d * d) / (2 * c * d)
        if abs(a - k) < d < a + k and -1 < cosine < 1:
            from_a_angle = d / (2 * a * k)
            from_b_energy = F(1, 2) / (6 * b * c * d)
            integral += (right - left) * from_a_angle * from_b_energy
    return integral / 4  # two unit-width radial densities of 1/2 each


def geometry(candidate):
    radii = tuple(map(F, ('1/8', '1/4', '1/2', '3/4', '1')))
    for args in product(radii, repeat=4):
        actual = candidate(*args)
        check(actual == sliced_density(*args), 'joint geometry normalization')
        _, a, b, c = args
        check(0 <= actual <= 1 / (48 * a * b * c), 'joint geometry envelope')


def fold_arithmetic():
    core = [F(1)] + [7 * 2 ** (2 * h - 1) * F(factorial(h), 2)
                     for h in (1, 2, 3)]
    closed = [sum(comb(m, h) * core[h] * 6 ** (m - h) * 11 ** (3 - m)
                  for h in range(m + 1)) for m in range(4)]
    # Enumerate each real variable's core/exterior status; integrate over
    # each possible location of the smallest core coordinate separately.
    independent = []
    for m in range(4):
        total = F(0)
        for core_flags in product((False, True), repeat=m):
            h = sum(core_flags)
            if h == 0:
                cost = 1
            else:
                one_minimum = F(factorial(h - 1), 2)
                cost = sum(2 ** (h - 1) * 14 * 2 ** (h - 1) * one_minimum
                           for flag in core_flags if flag)
            total += cost * 6 ** (m - h) * 11 ** (3 - m)
        independent.append(total)
    check(core == [1, 7, 56, 672], 'core fold constants')
    check(closed == independent == [1331, 1573, 1936, 2652], 'signed fold costs')


def odd_cancellation(use_abs=False):
    # The constant coefficient on a symmetric core has zero contribution
    # from any odd real factor, even before the continuum limit.
    grid = (F(-1, 2), F(-1, 4), F(1, 4), F(1, 2))
    real = lambda t: -t / (F(1, 16) + t * t)
    even = lambda t: F(1, 4) / (F(1, 16) + t * t)
    leg = (lambda t: abs(real(t))) if use_abs else real
    signed = sum(leg(a) * even(b) * even(c) for a, b, c in product(grid, repeat=3))
    factored = sum(map(real, grid)) * sum(map(even, grid)) ** 2
    check(signed == factored == 0, 'real factor must be folded before absolute value')


def row_identity(keep_width=True):
    kernels = ((7, 3, 1), (9, 8, 8, 4, 2, 1), (4, 1), (2,))
    for k in kernels:
        def transfer(j):
            return k[abs(j)] if abs(j) < len(k) else 0

        for n in range(12):
            row = 0
            for m in range(n + len(k)):
                if m != n:
                    width = 2 * m + 1 if keep_width else 1
                    row += (transfer(n - m) - transfer(n + m + 1)) * width
            expected = (2 * n + 1) * (2 * sum(k[1:]) + transfer(2 * n + 1))
            check(row == expected, 'weighted folded normal row')


def young_guard(r=F(1, 9)):
    h = F(1, 3)
    pi_lower = F(31, 10)
    c_upper = (6 / pi_lower) ** 2
    # Strict pi and integral guards in the proof give c < 4. The stronger
    # rational lower bound on pi gives an independent strict check here.
    check(r / (1 - h) + c_upper * r**2 / (1 - h)**3 < h, 'Young ball image')
    check(r / (1 - h)**2 + 3 * c_upper * r**2 / (1 - h)**4 < 1,
          'Young pairing/normal Schur margin')


def triple_holder():
    # Square-root factors have integer data; cubing removes the final root.
    fixtures = (([1, 2, 3], [2, 4, 1], [3, 1, 2]), ([1] * 3,) * 3)
    for f, g, h in fixtures:
        n = len(f)
        powered = (F(sum(f[a]**2 * g[c]**2 * h[(a+c) % n]**2
                         for a, c in product(range(n), repeat=2)), n*n)) ** 3
        separate = F(1)
        for factors in (f, g, h):
            separate *= F(sum(t**3 for t in factors), n) ** 2
        check(powered <= separate, 'translation-normalized triple Holder')
        if f == g == h:
            check(powered == separate, 'constant Holder extremizer')


def mixed_guard(keep_factor=True):
    delta = F(1, 17)
    for spend in (F(1, 2), F(3, 8)):
        mean = F(2)
        other = delta * mean / spend  # equality case of the completed square
        cost = delta**2 / spend if keep_factor else delta**2
        check(2 * delta * mean * other == spend * other**2 + cost * mean**2,
              'variance absorption must pay full squared cost')


def angular_guard():
    # An even nonconstant convolution on the four-point circle checks the
    # normalization/centering of the angular argument without a Fourier solver.
    a = F(1, 16)
    v = (1 + a, F(1), 1 - a, F(1))
    norm = lambda h: sum(t*t for t in h) / 4
    convolve = lambda h: tuple(sum(v[(i - j) % 4] * h[j] for j in range(4)) / 4
                               for i in range(4))
    eta = sum(abs(t - 1) for t in v) / 4
    check(convolve((1,) * 4) == (1,) * 4, 'angular constant sector')
    for h in ((1, 0, -1, 0), (0, 1, 0, -1), (-2, 1, -4, 5)):
        got = convolve(h)
        check(sum(h) == sum(got) == 0 and norm(got) <= eta**2 * norm(h),
              'angular zero-mean sector')
    check(1 - F(1, 16) * 2 == F(7, 8), 'angular ME floor')
    check(F(7, 8) - F(3, 8) - F(1, 8) == F(3, 8), 'angular crossed floor')


def constants():
    check(40 * 12**3 < 2**17 and F(2**34, 2**40) == F(1, 64), 'L^2 guard')
    check(8 * 63 * F(4, 9) < 256 and 63 * F(4, 9) == 28, 'line conversions')
    check(F(27, 64) < 1 and 2 * 9 * 258 < 5000, 'dressed-width guard for B >= 1')
    check(2 * 5000 + 4 < 11000, 'ME comparison guard for B >= 1')
    check(10 * (28 + F(27, 64) + 42) < 1000, 'projected crossed cost')
    check(1000 + 11000 + 5 < 13000, 'coercivity guard')
    check(F(13000, 2**80) < F(1, 4) and 4 * 13000 == 52000, 'ceiling inflation guard')
    check(F(9, 2) * F(1, 26) / (1 - F(9, 26**2)) < F(9, 50), 'weak ceiling guard')


def main():
    geometry(closed_density)
    fold_arithmetic()
    odd_cancellation()
    row_identity()
    young_guard()
    triple_holder()
    mixed_guard()
    angular_guard()
    constants()
    print('Exact geometry (625 fixtures), signed folding, normal rows, Young, angular and constants: accepted')
    mutants = {
        'wrong angular density normalization': lambda: geometry(
            lambda *args: closed_density(*args, prefactor=24)),
        'omit triangle admissibility': lambda: geometry(
            lambda *args: closed_density(*args, include_triangle=False)),
        'absolute real line before folding': lambda: odd_cancellation(use_abs=True),
        'omit width in folded normal sum': lambda: row_identity(keep_width=False),
        'unsupported Young threshold 1/8': lambda: young_guard(r=F(1, 8)),
        'underpay mixed variance absorption': lambda: mixed_guard(keep_factor=False),
    }
    for name, mutant in mutants.items():
        try:
            mutant()
        except AssertionError as exc:
            print(f'Mutant rejected: {name} ({exc})')
        else:
            raise AssertionError(f'mutant survived: {name}')


if __name__ == '__main__':
    main()
