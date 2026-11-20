"""Independent analytic-band ME control at density n=0.6 (energies in eV)."""
import json
import sys
import numpy as np
from scipy.linalg import eigvalsh
from scipy.optimize import brentq
from scipy.signal import fftconvolve

KB = 8.617333262145178e-5
FILL = .6
A0 = 1 - FILL


def band(u, a):
    """Uniform N(E)=1/2: averages of 1/(u²+(E+a)²) and Re G."""
    au = np.abs(u)
    q = (np.arctan((1+a)/au) - np.arctan((-1+a)/au))/(2*au)
    g = -.25*np.log((u*u+(1+a)**2)/(u*u+(-1+a)**2))
    return q, g


class ME:
    def __init__(self, nw, lam, omega):
        self.nw, self.lam, self.omega = nw, lam, omega
        self.history = []

    def evaluate(self, temperature):
        t = KB*temperature
        nw, lam, omega = self.nw, self.lam, self.omega
        w = np.pi*t*(2*np.arange(-nw,nw)+1)
        j = np.arange(1-2*nw,2*nw)
        v = 2*lam*omega**2/(omega**2+(2*np.pi*t*j)**2)
        candidates = [x for x in self.history if x[0] >= temperature]
        if candidates:
            oldtemp, oldw, u, h, dmu = min(candidates, key=lambda z:z[0]-temperature)
            u = u + w-oldw
            h = h.copy()
        else:
            u = w.copy()
            h = np.zeros_like(w)
            dmu = 0.
        for it in range(1200):
            aref = A0-dmu
            nref = t*(np.logaddexp(0., (1-aref)/t)-np.logaddexp(0., (-1-aref)/t))
            gref = band(w, np.full_like(w,aref))[1]
            q, g = band(u, A0-dmu+h)
            un = w + t*fftconvolve(u*q, v, mode='same')
            hn = t*fftconvolve(g, v, mode='same')
            density = nref + 4*t*np.sum((g-gref)[nw:])
            dn = FILL-density
            residual = max(np.max(abs(un-u)), np.max(abs(hn-h)), abs(dn))
            scale = max(np.max(abs(u)), np.max(abs(h)), abs(dmu), omega)
            if residual/scale <= 1e-10:
                break
            u += .2*(un-u)
            h += .2*(hn-h)
            dmu += .2*dn
        else:
            raise RuntimeError(f'ME continuum normal failed: relative_residual={residual/scale:g}')
        self.history.append((temperature, w, u.copy(), h.copy(), dmu))
        self.history = self.history[-24:]
        qpos = q[nw:]
        j = np.arange(nw)
        d = j[:,None]-j[None,:]
        s = j[:,None]+j[None,:]+1
        kernel = lambda z: 2*lam*omega**2/(omega**2+(2*np.pi*t*z)**2)
        matrix = t*np.sqrt(qpos[:,None]*qpos[None,:])*(kernel(d)+kernel(s))
        eig = eigvalsh(matrix, subset_by_index=[nw-1,nw-1])[0]
        return eig-1

    def solve(self):
        ceiling = .18272624777232952*self.omega*np.sqrt(self.lam)/KB
        hi = ceiling*1.8
        if self.evaluate(hi) >= 0:
            raise RuntimeError('upper bracket unstable')
        for _ in range(30):
            lo = hi/1.25
            if self.evaluate(lo) > 0:
                break
            hi = lo
        else:
            raise RuntimeError('lower bracket not found')
        root = brentq(self.evaluate,lo,hi,xtol=.002,rtol=1e-6)
        return {'nw':self.nw, 'lam':self.lam, 'omega':self.omega, 'T':root, 'ceiling':ceiling}


if __name__ == '__main__':
    nw = int(sys.argv[1])
    for lam in (.5,1.,2.):
        for omega in (.03,.06,.12):
            try:
                print(json.dumps(ME(nw,lam,omega).solve()), flush=True)
            except Exception as exc:
                print(json.dumps({'nw':nw,'lam':lam,'omega':omega,'error':repr(exc)}), flush=True)
