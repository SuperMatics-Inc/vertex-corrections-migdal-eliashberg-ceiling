"""Reduce one-axis fixed-density progress JSONL into convergence records.

The two shifts are empirical differences, not error bounds. Missing or failed
roots remain explicit; only complete grids receive the shift-based uncertainty.
"""
import json
import sys

MOMENTA = (130, 258, 514)
FREQUENCIES = (8, 16, 32, 64)
LAMBDAS = (.5, 1., 2.)
OMEGAS = (.03, .06, .12)
roots, failures = {}, {}
for path in sys.argv[1:]:
    with open(path, encoding='utf8') as f:
        for line in f:
            row = json.loads(line)
            if row['kind'] == 'root':
                key = (row['n'], row['nw'], row['lam'], round(row['omega'], 8), row['vertex'])
                roots[key] = row
            elif row['kind'] in ('error', 'deadline'):
                spec = row['spec'].split(':')
                if len(spec) >= 5:
                    key = (int(spec[0]), int(spec[1]), float(spec[2]), round(float(spec[3]), 8), bool(int(spec[4])))
                    failures[key] = row.get('error', 'deadline before completion')

for lam in LAMBDAS:
    for omega in OMEGAS:
        out = {'lambda': lam, 'omega_over_bare_fermi_energy': omega/.6}
        for vertex in (False, True):
            label = 'vertex' if vertex else 'ME'
            points = [(n, 32) for n in MOMENTA] + [(514, nw) for nw in FREQUENCIES if nw != 32]
            records = {f'{n}/{nw}': roots.get((n, nw, lam, omega, vertex), {}).get('T') for n, nw in points}
            errors = {f'{n}/{nw}': failures[(n, nw, lam, omega, vertex)] for n, nw in points if (n, nw, lam, omega, vertex) in failures}
            complete = all(v is not None for v in records.values())
            x = {'roots_K': records, 'failures': errors, 'all_axes_complete': complete}
            if complete:
                mshift = abs(records['514/32'] - records['258/32'])
                fshift = abs(records['514/64'] - records['514/32'])
                x.update(best_K=records['514/64'], momentum_shift_K=mshift,
                         frequency_shift_K=fshift, observed_uncertainty_K=mshift+fshift+.005,
                         ceiling_K=roots[(514,64,lam,omega,vertex)]['ceiling'],
                         ratio=roots[(514,64,lam,omega,vertex)]['ratio'])
            out[label] = x
        print(json.dumps(out, allow_nan=False), flush=True)
