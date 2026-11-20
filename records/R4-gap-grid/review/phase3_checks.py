"""Exact weak-edge constants and receipt/table cross-check, no verifier imports."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def weak():
    # A tighter rational upper endpoint for 2**(-3/4) than the author uses.
    q = F(119, 200)
    require(q**4 > F(1, 8), 'secant endpoint')
    j = 6 - F(22, 15) * (1 - q)
    require(j < F(406, 75) and j*j < 3 * F(333, 106)**2, 'DOS integral')
    r, a = F(3, 25), F(7, 10)
    image = r/a + 3*r*r/a**3
    row = r/a**2 + 9*r*r/a**4
    c = F('0.1827262477723295078215179')  # GAMMA lower certificate
    at_eta30 = F(7, 60) + 1/(12*c*c*30**2)
    at_adiabatic_edge = F(7, 60)*(1 + F(22, 7)**2/(12*2**6))
    require(image < F('0.29737610') < 1-a, 'weak image')
    require(row < F('0.78467306') < 1, 'weak derivative and pairing')
    require(max(at_eta30, at_adiabatic_edge) < r, 'weak finite windows')
    print('weak image, row, R_eta30 and R_adiabatic rounded upward:')
    for x in (image, row, at_eta30, at_adiabatic_edge):
        units = -(-x.numerator * 10**12 // x.denominator)
        print(f'{units // 10**12}.{units % 10**12:012d}')


def receipts(root):
    # Reconstruct peak/extra lists from the result table, not from grid_receipts.py.
    lam = tuple(map(F, ('1/9', '1/4', '1/2', '1', '2')))
    eta = F(1), F(3, 2)
    x = tuple(F(2)**i for i in range(-4, 5))
    grid = set(product(lam, eta, x))
    peak = {(l, e, u) for e, thresholds in zip(eta, ((1,2,4,8,16),(1,2,4,8,8)))
            for l, t in zip(lam, thresholds) for u in x if u >= t}
    rows2 = (('1/4 1/2', '1', '2', '4', '4 8'),
             ('1/8 1/4 1/2', '1', '2', '2 4', '4'))
    rows3 = (('1/16 1/8', '1/8 1/4 1/2', '1/4 1/2 1', '1/2 1 2', '2'),
             ('1/16', '1/16 1/8 1/4 1/2', '1/16 1/8 1/4 1/2 1', '1/8 1/4 1/2 1', '1 2'))

    def table(rows):
        return {(l, e, F(u)) for e, row in zip(eta, rows)
                for l, cell in zip(lam, row) for u in cell.split()}

    for file, expected_grid, expected_success, peak_only in (
        ('replay-weak', grid, peak, True),
        ('phase2-finite-rounding', grid-peak, table(rows2), False),
    ):
        found, passed = set(), set()
        for line in (root / (file + '.log')).read_text().splitlines():
            if not line[:1].isdigit():
                continue
            l, e, u, mark, *_ = line.split(',')
            if e == 'static':
                continue
            k = tuple(map(F, (l, e, u)))
            require(k not in found, 'duplicate earlier node')
            found.add(k)
            if not mark.startswith('UNRESOLVED'):
                require(mark == 'PEAK' if peak_only else F(0) <= F(mark) <= F('0.50417617'),
                        'earlier point bound')
                passed.add(k)
        require(found == expected_grid and passed == expected_success, 'earlier result table')

    files = ('map128', 'averaged-low', 'mid-growth', 'refined-low', 'refined-high',
             'last-slab', 'refined-weak512')
    found = {}
    for stem in files:
        for line in (root / (stem + '.log')).read_text().splitlines():
            if not line.startswith('CLOSED '):
                continue
            _, key, p, q = line.split()
            k = tuple(map(F, key.split(',')))
            found.setdefault(k, (F(p.removeprefix('pairing<=')), F(q.removeprefix('q<='))))
    require(found.keys() == table(rows3), 'radial result table')
    p, q = tuple(map(max, zip(*found.values())))
    require((p, q) == (F('0.72827964'), F('0.70929592')), 'radial reported maxima')
    combined = peak | table(rows2) | table(rows3)
    require(tuple(sum(k[1] == e for k in combined) for e in eta) == (34,40), 'combined counts')
    print(f'28 exact radial nodes; combined 34/45 and 40/45; maxima {float(p):.8f}, {float(q):.8f}')


if __name__ == '__main__':
    weak()
    receipts(Path('D:/t/G1Gap'))
