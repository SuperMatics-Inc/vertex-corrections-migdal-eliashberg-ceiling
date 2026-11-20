"""Exact arithmetic for the forward no-go and overlap contraction guards."""

from fractions import Fraction as F
from itertools import product


def forward_margin(delta, tail, window):
    # Rescale momentum to q_c=1 and Green lines to ν=1.
    r = 2 * delta
    assert 2 * r + 8 * r**2 <= 5 * delta
    u_min, u_max, e_max = 1 - 5 * delta, 1 + 5 * delta, 8 * delta
    assert u_min / (u_max**2 + e_max**2) >= F(1, 2)
    assert e_max / u_min <= 16 * delta
    assert tail <= 1
    static_negative = F(1, 4) * F(1, 2)**4 * window**2
    static_real = 3 * 2 * 2**4 * (16 * delta)**2
    nonstatic = (4 + 3 * 2) * 2**4 * 3 * tail
    return static_negative - static_real - nonstatic


def check_corner_geometry():
    # All extrema of linear coordinates on these boxes occur at vertices.
    support = F(1, 16)
    transfer = F(1, 2), F(3, 4)
    for p, q, ell in product((-support, support), transfer, transfer):
        assert max(p - q, p - ell, p - q - ell) < -support
    for p, q, ell in product((-support, support), (-1, 1), (-1, 1)):
        assert max(abs(p - q), abs(p - ell), abs(p - q - ell)) <= 2 + support
    return (transfer[1] - transfer[0]) / 2


def angular_integral(a, b):
    # I(a,b)/π by partial fractions in u=x², independent of parameter
    # differentiation of the integrated simple-pole product.
    if a == b:
        return (F(3, 8) - F(5, 16)) / a**3
    d = b*b - a*a
    double_a = -a*b*a*a / d**2
    double_b = -a*b*b*b / d**2
    simple_a = (a*b - 2*d*double_a) / d**2
    return simple_a * (1/a - 1/b) + double_a/(2*a**3) + double_b/(2*b**3)


def check_angular_mean():
    cutoff = 10
    n_line = range(-cutoff - 1, cutoff + 1)

    def kernel(j):
        return F(1, 1 + ((1 << 13) * j)**2)

    def trial(n):
        return int(n == 0) - int(n == -1)

    ordered, mutant = F(0), F(0)
    for n, m in product(n_line, repeat=2):
        a, b = F(abs(2*n + 1)), F(abs(2*m + 1))
        integral = angular_integral(a, b)
        assert integral == 1 / (2 * (a + b)**3)
        term = kernel(n - m) * (trial(n) - trial(m))**2 * integral
        ordered += term if (n < 0) == (m < 0) else -term
        mutant += term
    folded = -8 * kernel(1) / 16
    folded += 4 * sum((kernel(j) - kernel(j + 1)) / (2 * (2*j + 2)**3)
                      for j in range(1, cutoff + 1))
    assert ordered == folded and 2 * ordered <= -15 * kernel(1) / 16
    assert mutant > 0 > ordered
    assert -F(15, 16) + F(250, 1000) < -F(1, 4)


def weak_contraction(R):
    # π>3, sup x J³ <= 4096.
    assert F(2 * 4096, 3 * 9) < 304
    L = 304 * R**2
    return 3 * R + L < F(1, 2) and 4 * R + 6 * L < 1


if __name__ == '__main__':
    r = F(1, 16)
    assert 2 * r + 8 * r**2 == F(5, 32)
    assert 4 * r + 48 * r**2 == F(7, 16)
    assert 4 * F(1, 4) + 48 * F(1, 4)**2 > 1  # permissive-contraction mutant

    delta, tail = F(1, 2**16), F(1, 2**24)
    window = check_corner_geometry()
    margin = forward_margin(delta, tail, window)
    # Independent denominator clearing of the stated bound in (2).
    # Terms scaled by 2^32: 2^20 - 24576 - 480*256 - 2^19.
    scaled_slack = (1 << 20) - 24576 - 480 * 256 - (1 << 19)
    assert scaled_slack > 0 and margin - F(1, 8192) == F(scaled_slack, 1 << 32)
    assert forward_margin(delta, 8 * tail, window) <= F(1, 8192)
    assert forward_margin(delta, tail, window / 4) <= 0

    check_angular_mean()

    # Certified gamma=2 lower endpoint from ../GAMMA/RESULT.md, not its rounded value.
    c_inf_lower = F('0.1827262477723295078215179')
    assert c_inf_lower > F(1, 8)
    assert c_inf_lower * (1 << 16) > (1 << 12)
    assert F(1, 2**20) + 4 * F(1, 2**10) < F(1, 128)
    assert weak_contraction(F(1, 128))
    assert not weak_contraction(F(1, 16))

    # Independent rational Jacobian check for dξ_b/dU at test radii:
    # implicit dr_b/dU from the squared-radius law; then dξ_b/dR_b.
    for radius, rc, d, width in [(F(1, 2), F(1, 3), F(1, 5), F(7, 3)),
                                  (F(7, 9), F(7, 8), F(3, 4), F(2))]:
        d_r = (2 * rc * d) / (2 * radius)
        d_energy = 2 * width * 3 * radius**2 * d_r
        assert d_energy == 6 * width * radius * rc * d
        assert F(1, 2) / d_energy / (2 * width)**2 == 1 / (48 * width**3 * radius * rc * d)
    print(f'Exact forward margin {margin} > 1/8192; angular-mean and overlap guards; five mutants rejected.')
