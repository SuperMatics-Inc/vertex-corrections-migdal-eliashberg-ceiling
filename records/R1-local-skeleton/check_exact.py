"""Exact finite algebra replay, not an approximation to the infinite operator."""
from fractions import Fraction as F
from itertools import product


def corner_forms(weights, widths, gap, *, middle=-1, normal=1):
    size = len(weights)
    direct = factored = zero_transfer = 0
    for n, j, l in product(range(size), repeat=3):
        a, b, c = (n - j) % size, (n - j - l) % size, (n - l) % size
        coefficient = weights[j] * weights[l] * widths[n] * widths[a] * widths[b] * widths[c]
        row = coefficient * gap[n] * (normal * gap[n] - (gap[a] + middle * gap[b] + gap[c]))
        direct += row
        factored += coefficient * (gap[n] - gap[a] + gap[b] - gap[c]) ** 2
        if j == 0 or l == 0:
            zero_transfer += row
    return direct, F(factored, 4), zero_transfer


def row_moments():
    transfer = list(map(F, [21, 8, 4, 4, 1, 0]))

    def k(j):
        return transfer[j] if j < len(transfer) else 0

    for n in range(12):
        terms = [(m, k(abs(n - m)) - k(n + m + 1)) for m in range(n + len(transfer)) if m != n]
        w = 2 * n + 1
        assert sum(p for _, p in terms) <= w * k(1)
        assert sum((2 * m + 1) * p for m, p in terms) == w * (2 * sum(transfer[1:]) + k(w))


def main():
    weights = [2, 3, 5, 7, 7, 5, 3]
    widths = [3, -2, 1, -3, 2, 4, -1]
    gap = [6, 5, -1, 2, 0, 1, -3]
    direct, factored, static = corner_forms(weights, widths, gap)
    assert direct == factored != 0 and static == 0
    for mutant in ({"middle": 1}, {"normal": 0}):
        mutated, _, static = corner_forms(weights, widths, gap, **mutant)
        assert mutated != factored and static != 0
    assert direct != 4 * factored
    row_moments()
    assert 4 * F(157, 50) ** 2 * F(9, 50) ** 2 > F(5, 4)
    assert F(1, 24) / F(9, 50) ** 2 < F(4, 3)
    assert F(1, 24) / F(9, 50) < F(1, 4)
    assert F(4, 5) + F(1, 15) < F(7, 8)
    assert 2 * F(16, 3) / 100 < 1
    ratio = F("134.27") / F("129.85")
    print("exact algebra replay: signed form, row moments, rational guards, three mutants")
    print(f"G1 rounded-input ratio = {ratio} = {float(ratio):.12f}; excess = {float(100*(ratio-1)):.6f}%")


if __name__ == "__main__":
    main()
