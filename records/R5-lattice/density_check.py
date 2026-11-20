"""Free-density acceptance and bare-reference mutant at fixed density."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
import lattice

n, fill, temperature, omega = 66, .6, 34.749, .03
for nw in (8,16,32):
    x,w,_,t = lattice.setup(n,nw,0.,omega,temperature)
    x = x+1-fill
    u = np.broadcast_to(w[:,None],(2*nw,n)).copy()
    h = np.zeros_like(u)
    exact = lambda d: 2*np.mean(expit(-(x-d)/t))
    proper = lambda d: lattice.density(u,h,x,w,t,nw,d)
    def mutant(d):
        ref = -x[None,:]/(w[:,None]**2+x[None,:]**2)
        xx = x-d
        green = -xx[None,:]/(w[:,None]**2+xx[None,:]**2)
        return 2*np.mean(expit(-x/t))+4*t*np.sum(np.mean(green[nw:]-ref[nw:],axis=1))
    d = brentq(lambda z:proper(z)-fill,-.5,.5,xtol=1e-14)
    dm = brentq(lambda z:mutant(z)-fill,-.5,.5,xtol=1e-14)
    err = abs(proper(d)-exact(d))
    merr = abs(mutant(dm)-exact(dm))
    print(json.dumps(dict(nw=nw,mu_shift=d,density=proper(d),exact=exact(d),error=err,
                          mutant_mu_shift=dm,mutant_density=mutant(dm),
                          mutant_exact=exact(dm),mutant_error=merr)),flush=True)
    if err>1e-12 or abs(proper(d)-fill)>1e-12 or merr<=1e-12:
        raise SystemExit('density acceptance or mutant failed')
