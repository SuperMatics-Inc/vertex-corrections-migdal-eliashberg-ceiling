"""Independent zero-grid ME control for the uniform finite-band DOS."""
import json
import numpy as np
from scipy.linalg import eigvalsh
from scipy.optimize import brentq

KB = 8.617333262145178e-5
C = 0.18272624777232952

def eigenvalue(temp, lam, omega, nw=96):
    t = KB * temp
    j = np.arange(nw)
    w = np.pi * t * (2*j + 1)
    diff = j[:, None] - j[None, :]
    summ = j[:, None] + j[None, :] + 1
    def interaction(z):
        return 2*lam*omega**2/(omega**2+(2*np.pi*t*z)**2)
    vd, vs = interaction(diff), interaction(summ)
    u = w.copy()
    for _ in range(1000):
        q = np.arctan(1/u)/u
        un = w + t*((vd-vs)@(u*q))
        residual = np.max(np.abs(un-u))
        u = 0.65*u + 0.35*un
        if residual < 2e-12:
            break
    else:
        raise RuntimeError('normal solve did not converge')
    q = np.arctan(1/u)/u
    mat = t*np.sqrt(q[:,None]*q[None,:])*(vd+vs)
    return float(eigvalsh(mat, subset_by_index=[nw-1,nw-1])[0])

def main():
    for lam in (0.5,1.0,2.0):
        for omega in (0.05,0.1,0.2):
            ceiling = C*omega*np.sqrt(lam)/KB
            hi = 1.8*ceiling
            if eigenvalue(hi,lam,omega) >= 1:
                raise RuntimeError('upper bracket unstable')
            while True:
                lo = hi/1.5
                if eigenvalue(lo,lam,omega) > 1:
                    break
                hi = lo
                if lo < .001:
                    raise RuntimeError('lower bracket not found')
            root = brentq(lambda temp: eigenvalue(temp,lam,omega)-1,lo,hi,xtol=1e-5)
            print(json.dumps({'lambda':lam,'omega':omega,'nw':96,'T':root,'ceiling':ceiling}),flush=True)

if __name__ == '__main__':
    main()
