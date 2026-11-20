"""Check particle-hole covariance of the direct and FFT normal maps."""
import numpy as np
import lattice

n, nw = 6, 2
x, w, v, t = lattice.setup(n, nw, 1., .2, 100.)
u = np.broadcast_to(w[:, None], (2 * nw, n)).copy()
h = np.zeros_like(u)
d = u * u + (x[None, :] + h) ** 2
z = -u / d
q = -(x[None, :] + h) / d
C = lambda a, b, c: lattice.direct(a, b, c, v)
zd = lattice.meanconv(z, v) + t * (-C(z, z, z) + C(q, z, q) + C(q, q, z) + C(z, q, q))
hd = lattice.meanconv(q, v) + t * (-C(q, z, z) - C(z, z, q) - C(z, q, z) + C(q, q, q))
ud, hd = w[:, None] - t * zd, t * hd
uf, hf = lattice.normal(u, h, x, w, v, t)
for name, uu, hh in [('direct', ud, hd), ('fft', uf, hf)]:
    pu = np.max(np.abs(uu - np.roll(uu, n // 2, axis=1)))
    ph = np.max(np.abs(hh + np.roll(hh, n // 2, axis=1)))
    print(f'{name} u_ph_defect={pu:.12g} h_ph_defect={ph:.12g}', flush=True)
print(f'direct_fft={max(np.max(np.abs(ud-uf)), np.max(np.abs(hd-hf))):.12g}', flush=True)
