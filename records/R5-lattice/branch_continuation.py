"""Fixed-T fixed-density continuation from lambda=0 for the crossed onset."""
import json
import sys
import numpy as np
from scipy.optimize import brentq
import lattice

n, nw = int(sys.argv[1]), int(sys.argv[2])
target, omega = map(float, sys.argv[3:5])
fill = .6

def branch(T):
    # The lambda=0 reference fixes the same finite-Matsubara number equation.
    x, w, _, t = lattice.setup(n, nw, 0., omega, T)
    x = x + 1-fill
    u0 = np.broadcast_to(w[:,None], (2*nw,n)).copy()
    h0 = np.zeros_like(u0)
    def number(delta):
        return lattice.density(u0,h0,x,w,t,nw,delta)-fill
    delta0 = brentq(number,-.5,.5,xtol=1e-12)
    lattice.Solver.cache.clear()
    lattice.Solver.cache[(n,nw,omega,True,fill,0.)] = [(T,(w,u0,h0,delta0))]
    steps = int(np.ceil(target/.05))
    for lam in np.linspace(target/steps,target,steps):
        solver = lattice.Solver(n,nw,float(lam),omega,True,fill)
        val = solver.evaluate(T)
        last = solver.history[-1][1][3]
    lattice.log(kind='lambda_branch',n=n,nw=nw,lam=target,omega=omega,T=T,
                eigenvalue=val+1,mu_shift=float(last),steps=steps)
    return val

ceiling = lattice.C*omega*np.sqrt(target)/lattice.KB
try:
    if len(sys.argv)>5:
        branch(float(sys.argv[5]))
    else:
        hi = 1.8*ceiling
        if branch(hi) >= 0:
            raise RuntimeError('upper bracket unstable')
        for _ in range(30):
            lo = hi/1.25
            if branch(lo) > 0:
                break
            hi = lo
        else:
            raise RuntimeError('lower bracket not found')
        root = brentq(branch,lo,hi,xtol=.005,rtol=1e-5)
        residual = branch(root)
        if abs(residual)>1e-4:
            raise RuntimeError(f'branch root is discontinuous or inaccurate: eigenvalue-1={residual}')
        lattice.log(kind='lambda_branch_root',n=n,nw=nw,lam=target,omega=omega,T=root,
                    ceiling=ceiling,ratio=root/ceiling,eigenvalue_error=residual)
except Exception as exc:
    lattice.log(kind='lambda_branch_failure',n=n,nw=nw,lam=target,omega=omega,error=str(exc))
    raise
