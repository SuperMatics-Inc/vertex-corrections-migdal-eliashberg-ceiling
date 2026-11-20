"""Independent floating-point sine-transform diagnostic, not a certificate."""
import argparse
import os

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from scipy.fft import dst, idst
from scipy.linalg import eigh
from scipy.special import spherical_jn

CEIL = 0.18272624777232952


class Band:
    def __init__(self, cells, count):
        self.count = count
        self.cells = cells
        self.dr = 1 / cells
        self.length = 4
        self.p = np.arange(1, 4 * cells) * np.pi / self.length
        self.r = np.arange(1, cells + 1) * self.dr
        self.dispersion = 2 * self.r**3 - 1
        self.sinc = np.sinc(self.r[:, None] * self.p[None, :] / np.pi)
        self.weight = 3 * self.r**2 * self.dr
        self.weight[-1] /= 2
        self.full = np.arange(-count, count)
        self.pos = np.arange(count)
        self.negative = self.full < 0
        self.width = 2 * self.pos + 1

    def transform(self, z):
        f = np.zeros((*z.shape[:-1], len(self.p)), dtype=complex)
        f[..., :self.cells] = self.r * z
        f[..., self.cells - 1] /= 2
        return dst(f, type=1) * (1.5 * self.dr / self.p)

    def inverse(self, g):
        return idst(g * (self.p / (1.5 * self.dr)), type=1)[..., :self.cells] / self.r

    def transfer(self, lam, x, ratio):
        # ratio is T/Tceil; infinity means Omega/x=0.
        m = self.full[:, None]
        q = self.full[None, :]
        def g(j):
            square = (2 * np.pi * CEIL * np.sqrt(lam) * ratio)**2
            if np.isinf(square):
                return (2 * lam * x / np.pi) * (j == 0)
            return (2 * lam * x / np.pi) / (1 + square * j**2)
        normal = []
        ends = []
        middle = []
        for n in self.pos:
            b = m + q - n
            a = n + m - q
            valid_b = (b >= -self.count) & (b < self.count)
            valid_a = (a >= -self.count) & (a < self.count)
            normal.append((np.clip(b + self.count, 0, 2*self.count-1),
                           g(n-m)*g(n-q)*valid_b))
            ends.append(g(n-q)*valid_b)
            middle.append((np.clip(a + self.count, 0, 2*self.count-1),
                           g(n-q)*g(q-m)*valid_a))
        return g(self.pos[:, None]-self.full[None, :]), normal, ends, middle

    def reflect(self, z):
        return np.concatenate((z[::-1].conj(), z))

    def solve(self, lam, x, ratio):
        rainbow, normal, ends, middle = self.transfer(lam, x, ratio)
        bare = 1j * x * self.width[:, None] - self.dispersion
        sigma = np.zeros_like(bare)
        for alpha in (0.25, 0.5, 0.75, 1):
            for sweeps in range(360):
                z = 1 / (bare - sigma)
                if np.min((bare - sigma).imag) <= 0:
                    raise ValueError('nonpositive width')
                zz = self.reflect(z)
                zh = self.transform(zz)
                crossed = np.asarray([
                    np.einsum('ij,ik,jk,ijk->k', w, zh, zh, zh[b], optimize=True)
                    for b, w in normal
                ])
                update = alpha * (rainbow @ (zz @ self.weight))[:, None]
                update = update + alpha**2 * self.inverse(crossed)
                err = np.max(abs(update - sigma)) / x
                if err < 2e-11:
                    break
                sigma = 0.6 * update + 0.4 * sigma
            else:
                raise ValueError(f'iteration did not converge: {err:.3g}')
        return sigma, zh, (rainbow, normal, ends, middle)

    def eigenvalues(self, z, zh, transfers, sectors):
        rainbow, normal, ends, middle = transfers
        norm = np.sqrt(self.weight)[None, :] * abs(z)
        out = np.empty((self.count, self.count, len(self.p)))
        mult = self.p**2 * (2 / (3 * self.length))
        for n in self.pos:
            b, _ = normal[n]
            a, w = middle[n]
            end = np.einsum('ij,jk,ijk->ik', ends[n], zh, zh[b], optimize=True).real
            mid = np.einsum('ij,jk,ijk->ik', w, zh, zh[a].conj(), optimize=True).real
            lines = 2 * rainbow[n, :, None] * end + mid
            for m in self.pos:
                pair = (lines[self.count + m] + lines[self.count - 1 - m]) * mult
                out[n, m] = pair
        result = []
        for ell in sectors:
            basis = self.sinc if ell == 0 else spherical_jn(ell, self.r[:,None]*self.p)
            matrix = np.empty((self.count, self.cells, self.count, self.cells))
            for n in self.pos:
                for m in self.pos:
                    kernel = (basis * out[n,m]) @ basis.T
                    if ell == 0:
                        kernel += rainbow[n,self.count+m] + rainbow[n,self.count-1-m]
                    matrix[n,:,m,:] = kernel * norm[n,:,None] * norm[m,None,:]
            matrix = matrix.reshape(self.count*self.cells, -1)
            asym = np.max(abs(matrix-matrix.T))
            if asym > 1e-11:
                raise ValueError(f'pair symmetry violated: {asym:.3g}')
            size = len(matrix)
            result.append(eigh(matrix, eigvals_only=True, subset_by_index=(size-1, size-1))[0])
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells', type=int, default=48)
    parser.add_argument('--count', type=int, default=6)
    parser.add_argument('--xs', nargs='+', type=float, default=[0.04, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16])
    parser.add_argument('--lam', nargs='+', type=float, default=[1/9, .25, .5, 1, 2])
    parser.add_argument('--ratio', nargs='+', type=float, default=[float('inf'), 1, 1.5])
    parser.add_argument('--sectors', nargs='+', type=int, default=[0])
    args = parser.parse_args()
    if not (8 <= args.cells <= 512 and 1 <= args.count <= 32 and args.cells*args.count <= 3200):
        parser.error('grid exceeds the diagnostic resource limit')
    if not all(l >= 0 and l % 2 == 0 for l in args.sectors):
        parser.error('only even singlet angular sectors are supported')
    headers = ','.join(f'max_eig_l{l}' for l in args.sectors)
    print(f'lambda,x,T/Tceil,{headers},min_u_over_omega', flush=True)
    unresolved = 0
    for ratio in args.ratio:
        for lam in args.lam:
            for x in args.xs:
                band = Band(args.cells, 1 if np.isinf(ratio) else args.count)
                try:
                    sigma, zh, data = band.solve(lam, x, ratio)
                    shift = 1j * x * band.width[:, None] - band.dispersion - sigma
                    ev = band.eigenvalues(1/shift, zh, data, args.sectors)
                    values = ','.join(f'{v:.10f}' for v in ev)
                    print(f'{lam:.12g},{x:g},{ratio:g},{values},{np.min(shift.imag/(x*band.width[:,None])):.8f}', flush=True)
                except ValueError as error:
                    print(f'{lam:.12g},{x:g},{ratio:g},UNRESOLVED:{error}', flush=True)
                    unresolved += 1
    if unresolved:
        raise SystemExit(f'{unresolved} diagnostic nodes unresolved')


if __name__ == '__main__':
    main()
