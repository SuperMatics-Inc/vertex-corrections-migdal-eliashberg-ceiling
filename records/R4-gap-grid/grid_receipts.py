"""Aggregate completed replay receipts; this checks grid bookkeeping, not operators."""
import argparse
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import re


LAM = tuple(map(Q, ('1/9', '1/4', '1/2', '1', '2')))
ETA = (Q(1), Q(3, 2))
X = tuple(Q(2)**i for i in range(-4, 5))
GRID = set(product(LAM, ETA, X))
RADIAL = 'map128 averaged-low mid-growth refined-low refined-high last-slab refined-weak512'.split()
# Historical receipts say q; new receipts say q_G for the same preconditioned map.
CLOSED = re.compile(r'^CLOSED (\S+) pairing<=(\S+) q(?:_G)?<=(\S+)$')
PUBLISHED = (
    ('1/16 1/8', '1/8 1/4 1/2', '1/4 1/2 1', '1/2 1 2', '2'),
    ('1/16', '1/16 1/8 1/4 1/2', '1/16 1/8 1/4 1/2 1', '1/8 1/4 1/2 1', '1 2'),
)
PUBLISHED_NODES = {(lam, eta, Q(x)) for eta, rows in zip(ETA, PUBLISHED)
                   for lam, xs in zip(LAM, rows) for x in xs.split()}


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def series(grid, lam, eta):
    return ','.join(str(x) for x in X if (lam, eta, x) in grid) or 'none'


def radial_receipts(logs, names, left):
    added = {}
    for name in names:
        seen = set()
        for row in (logs / f'{name}.log').read_text().splitlines():
            match = CLOSED.match(row)
            if match is None:
                require(not row.startswith('CLOSED '), f'invalid enclosure in {name}: {row}')
                continue
            key, pair, normal = match.groups()
            node = tuple(map(Q, key.split(',')))
            require(node in left and node not in seen, f'invalid or duplicate radial node {node}')
            seen.add(node)
            a, q = Q(pair), Q(normal)
            require(0 <= a < 1 and 0 <= q < 1, f'unstable or invalid bound in {name}')
            added.setdefault(node, (a, q))
    return added


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('logdir', type=Path)
    parser.add_argument('--radial', action='append', default=[], metavar='LOG_STEM')
    opts = parser.parse_args()
    logs = opts.logdir
    peak, finite = {}, {}
    for name in ('replay-weak', 'phase2-finite-rounding'):
        seen = set()
        for row in (logs / f'{name}.log').read_text().splitlines():
            if not row[:1].isdigit():
                continue
            l, e, x, mark, *rest = row.split(',')
            if e == 'static':
                continue
            node = tuple(map(Q, (l, e, x)))
            require(node in GRID and node not in seen, f'invalid or duplicate node {node}')
            seen.add(node)
            if mark.startswith('UNRESOLVED'):
                continue
            if name == 'replay-weak':
                require(mark == 'PEAK', f'unexpected box guard {mark}')
                peak[node] = Q(rest[1])
            else:
                finite[node] = Q(mark)
        expected = GRID if name == 'replay-weak' else GRID - peak.keys()
        require(seen == expected, f'incomplete grid in {name}')
    require(len(peak) == 31 and len(finite) == 15, 'phase 1 or 2 count changed')
    left = GRID - peak.keys() - finite.keys()
    added = radial_receipts(logs, RADIAL, left)
    require(all(0 <= a < 1 for a in (*peak.values(), *finite.values())), 'invalid earlier bound')
    require(len(PUBLISHED_NODES) == 28 and added.keys() == PUBLISHED_NODES,
            'historical radial receipts differ from the exact published 28-node set')
    bounds = (max(a for a, _ in added.values()), max(q for _, q in added.values()))
    require(bounds[0] <= Q('0.72827964'), 'historical pairing maximum exceeds 0.72827964')
    require(bounds[1] <= Q('0.70929592'), 'historical q_G maximum exceeds 0.70929592')
    total = peak.keys() | finite.keys() | added.keys()
    for e in ETA:
        closed = sum(t[1] == e for t in total)
        print(f'eta={e} closed={closed}/45 radial={sum(t[1] == e for t in added)}')
        for l in LAM:
            print(f'  lambda={l} radial_x={series(added,l,e)} open_x={series(GRID-total,l,e)}')
    # Bounds use the first qualifying historical receipt per node, not a fresh grid replay.
    units = tuple(-(value.numerator * -10**8 // value.denominator) for value in bounds)
    print(f'phase3_historical_receipt_union: pairing<=0.{units[0]:08d} '
          f'preconditioned_normal_map_G_q<=0.{units[1]:08d}')
    if opts.radial:
        extra = radial_receipts(logs, opts.radial, left)
        new_nodes = extra.keys() - PUBLISHED_NODES
        print(f'optional_receipts: closed={len(extra)}, new={len(new_nodes)} (not in published union)')


if __name__ == '__main__':
    main()
