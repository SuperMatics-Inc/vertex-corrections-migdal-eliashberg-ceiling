"""Compare Matsubara/momentum parity sectors at a density-fixed normal state."""
import sys
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs
import lattice

n, nw = int(sys.argv[1]), int(sys.argv[2])
lam, om, temperature = map(float, sys.argv[3:6])
solver = lattice.Solver(n, nw, lam, om, True, .6)
initial = solver.evaluate(temperature) + 1
w, u, h = solver.prev
delta = solver.history[-1][1][3]
x, _, v, t = lattice.setup(n, nw, lam, om, temperature)
x = x + .4 - delta
rev = (-np.arange(n)) % n
print(f'normal_inversion={max(np.max(abs(u-u[:,rev])), np.max(abs(h-h[:,rev]))):.12g}', flush=True)
print(f'normal_frequency={max(np.max(abs(u+u[::-1,:])), np.max(abs(h-h[::-1,:]))):.12g}', flush=True)
apply = lattice.gap(u, h, x, v, t)
if n % 2:
    raise ValueError("sector basis requires even momentum grid")
rng = np.random.default_rng(7182)
for freq_sign in (1, -1):
    for momentum_sign in (1, -1):
        # Work in a reduced basis: projection in the full space would add
        # spurious zero eigenvalues which can dominate a negative sector.
        ks = np.arange(0 if momentum_sign == 1 else 1,
                       n//2+1 if momentum_sign == 1 else n//2)
        reduced_shape = (nw, len(ks))
        size = nw*len(ks)
        def matvec(z):
            pos = np.zeros((nw,n))
            pos[:,ks] = z.reshape(reduced_shape)
            pos[:,(-ks)%n] = momentum_sign*z.reshape(reduced_shape)
            full = np.concatenate((freq_sign*pos[::-1,:], pos),axis=0)
            y = apply(full)
            projected = (y + freq_sign*y[::-1,:] + momentum_sign*y[:,rev]
                         + freq_sign*momentum_sign*y[::-1,:][:,rev])/4
            return projected[nw:,:][:,ks].ravel()
        op = LinearOperator((size,size), matvec=matvec, dtype=float)
        eig = eigs(op, k=1, which='LR', tol=2e-5, maxiter=200,
                   v0=rng.standard_normal(size), return_eigenvectors=False)[0]
        print(f'frequency_parity={freq_sign:+d} momentum_parity={momentum_sign:+d} eigenvalue={eig.real:.12g} imag={eig.imag:.3g}', flush=True)
print(f'uniform_seed={initial:.12g}', flush=True)
