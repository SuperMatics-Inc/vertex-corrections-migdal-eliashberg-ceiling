"""Independent exact guard checks for G1-ONSET; not a continuum certificate."""

from fractions import Fraction as F
from itertools import product
from math import factorial


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


def young_guards():
    # Split the integral into |t|<=1/2, 1/2<|t|<=1 and |t|>1.
    # The second band has total length one and integrand <9/10.
    require(F(4, 5) ** 3 < F(9, 10) ** 4, 'middle-band integral cap')
    j_upper = 1 + F(9, 10) + 4
    c_upper = (j_upper / 3) ** 2  # pi > 3
    h, r = F(1, 3), F(1, 9)
    q = r / (1 - h) ** 2
    # One more resolvent factor yields the crossed Lipschitz/Schur row.
    image = (1 - h) * (q + c_upper * q**2)
    row = q + 3 * c_upper * q**2
    require(q == F(1, 4) and c_upper < 4, 'Young threshold scaling')
    require((1 - h) * (q + 4 * q**2) == h, 'endpoint image at c=4')
    require(q + 12 * q**2 == 1, 'endpoint row at c=4')
    require(image < h and row < 1, 'true strict guards at r=1/9')
    t = F(1, 26)  # t = sqrt(A*lambda); t/(1-9*t^2) is increasing below 1/3.
    crossover = F(9, 2) * t / ((1 - 3 * t) * (1 + 3 * t))
    require(crossover == F(117, 667) and crossover < F(9, 50), 'crossover')
    print(f'Young r={r}: image <= {image}, Schur <= {row}; crossover {crossover}')


def signed_and_rows():
    folds = []
    for n_real in range(4):
        total = 0
        for core in product((False, True), repeat=n_real):
            h = sum(core)
            # Each minimum coordinate has integral (h-1)!/2. Pair the sign
            # choices at that coordinate; the other radii cost two each.
            cost = h * F(factorial(h - 1), 2) * 14 * 4 ** (h - 1) if h else 1
            total += cost * 6 ** (n_real - h) * 11 ** (3 - n_real)
        folds.append(total)
    require(folds == [1331, 1573, 1936, 2652], 'core/exterior decomposition')
    require(max(folds) / F(48) < 63, 'signed integral scale')
    # Expand by one positive transfer index j at a time, including the
    # exceptional j=2*n+1 whose Hankel destination is the omitted diagonal.
    for n, j in product(range(20), range(1, 60)):
        w = 2 * n + 1
        coeffs = {}
        for m, sign in ((n + j, 1), (n - j, 1), (j - n - 1, -1)):
            if m >= 0 and m != n:
                coeffs[m] = coeffs.get(m, 0) + sign
        require(sum(coeffs.values()) == 2 * (j <= n) + (j == w), 'row count')
        require(sum((2 * m + 1) * c for m, c in coeffs.items())
                == w * (2 + (j == w)), 'weighted row coefficient')
    print(f'Signed fold totals: {", ".join(map(str, folds))}; per-transfer row checks passed')


def spherical_guards():
    # Suppress powers of B>=1 only after assigning each term its power.
    x_cap = F(1, 2**80)  # x <= x_cap * B^-8
    l_cap = 2**17 * F(1, 2**60)  # L <= l_cap * B^-4
    l_square_over_x = 2**34 * F(1, 2**40)
    require(40 * 12**3 < 2**17 and l_square_over_x == F(1, 64), 'L guards')
    require(F(4, 3) * x_cap + F(32, 3) * l_cap + F(1, 4) < 1, 'split f row')
    require(6 * 17 * l_cap < 1 and x_cap / 3 + 4 * l_cap < F(1, 4), 'split box')
    require(9 * l_cap < F(1, 8) and 10 * l_cap < F(1, 8), 'd and delta')

    monomial = F(4, 9) * 63  # upper coefficient of B^2*x using pi > 3
    replace = 3 * 9 * l_square_over_x  # upper coefficient of B*x
    require(8 * monomial < 256 and 8 * monomial + replace < 257, 'normal core')
    width = 2 * (8 + 1) * (257 + 1)
    require(width <= 5000, 'dressed width coefficient of B^3*x')
    require((1 - F(1, 8)) / (1 + F(1, 8)) ** 2 > F(1, 2), 'Q lower ratio')
    require(F(32, 3) < 14, 'exterior mass coefficient')
    scalar = 2 * 5000 + F(4, 3)  # coefficient of B^4*x
    projected = 10 * (monomial + replace + 14 * 3)  # use L<=1
    require(scalar < 11000 and projected < 1000, 'mean block bounds')
    mixed = []
    for floor, spend in ((F(1), F(1, 2)), (F(7, 8), F(3, 8))):
        require(floor - spend - F(1, 8) == F(3, 8), 'variance margin')
        mixed.append(10**2 * l_square_over_x / spend)
    require(mixed == [F(25, 8), F(25, 6)], 'mixed costs divided by x')
    c0 = 13000
    require(11000 + 1000 + max(mixed) < c0, 'coercivity coefficient of B^4*x')
    require(c0 * x_cap < F(1, 4) and 4 * c0 == 52000, 'inflation coefficient')
    for a in (F(1, 2**80), F(1, 2), F(1)):
        guard = 1 - (1 + a) ** -2 - a / 2
        require(guard == a * (3 - a*a) / (2 * (1 + a)**2) > 0, 'inflation gap')
    print(f'Variance 3/8 in both sectors; mixed costs {mixed[0]}, {mixed[1]}; inflation {4*c0}*B^4')


def overlap_guards():
    require(F(1, 5)**2 * 12 * F(1, 9) < 1, 'x_c strictly below x_h under (8)')
    for cap in (F(1, 100), F(1, 10), F(1, 9) - F(1, 10**12)):
        # Synthetic rho with pi*rho/x_h = a*sqrt(12*gap). Test below the join.
        a, gap = F(99, 100), F(1, 9) - cap
        scaled_x = (1 + a) / 2
        r = cap + gap * (a / scaled_x) ** 2
        require(a < scaled_x < 1 and r < F(1, 9), 'strict overlap below x_h')
    require(F(1, 9) + F(1, 8) < F(1, 3), 'uniform overlap norm')
    require(F(1, 9) + F(1, 9) < F(1, 3), 'Einstein overlap norm')
    print('Overlap strictly below join; both small-x norm bounds fit Young ball')


if __name__ == '__main__':
    young_guards()
    signed_and_rows()
    spherical_guards()
    overlap_guards()
