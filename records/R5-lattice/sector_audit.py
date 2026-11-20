"""Compare the uniform seed with random Pauli-allowed pairing sectors."""
import numpy as np
import lattice

original = lattice.eigs


def seed_for(shape, sector):
    rows, n = shape
    a = np.random.default_rng(4342).standard_normal(shape)
    rev = (-np.arange(n)) % n
    w = a[::-1, :]
    k = a[:, rev]
    wk = w[:, rev]
    if sector == 'even-even':
        return (a + w + k + wk).ravel()
    return (a - w - k + wk).ravel()


def run(n, nw, lam, om, temperature, sector):
    def seeded(op, **kwargs):
        if sector != 'uniform':
            kwargs['v0'] = seed_for((2 * nw, n), sector)
        return original(op, **kwargs)
    lattice.eigs = seeded
    try:
        solver = lattice.Solver(n, nw, lam, om, True)
        value = solver.evaluate(temperature) + 1
        _, u, h = solver.prev
        pherr = max(np.max(np.abs(u - np.roll(u, n // 2, axis=1))),
                    np.max(np.abs(h + np.roll(h, n // 2, axis=1))))
        print(f'{n=} {nw=} {lam=} {om=} {temperature=} {sector=} eigenvalue={value:.12g} particle_hole_error={pherr:.3g}', flush=True)
    finally:
        lattice.eigs = original


for pars in [(66, 8, 1, .1, 104.), (66, 8, 2, .2, 185.)]:
    for sector in ('uniform', 'even-even', 'odd-odd'):
        run(*pars, sector)
