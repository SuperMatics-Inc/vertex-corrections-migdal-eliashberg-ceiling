"""Onsets from independently cold-started normal solves; detect branch jumps."""
import sys
import numpy as np
from scipy.optimize import brentq
import lattice

n,nw=int(sys.argv[1]),int(sys.argv[2])
omegas=(float(sys.argv[4]),) if len(sys.argv)>4 else (.03,.06,.12)
lambdas=(float(sys.argv[3]),) if len(sys.argv)>3 else (.5,1.,2.)
for omega in omegas:
    for lam in lambdas:
        ceiling=lattice.C*omega*np.sqrt(lam)/lattice.KB
        def evaluate(T):
            lattice.Solver.cache.clear()
            return lattice.Solver(n,nw,lam,omega,True,.6).evaluate(T)
        try:
            hi=1.8*ceiling
            if evaluate(hi)>=0: raise RuntimeError('upper bracket unstable')
            for _ in range(30):
                lo=hi/1.25
                if evaluate(lo)>0: break
                hi=lo
            else: raise RuntimeError('lower bracket not found')
            root=brentq(evaluate,lo,hi,xtol=.005,rtol=1e-5)
            err=evaluate(root)
            if abs(err)>1e-4:
                raise RuntimeError(f'root may be branch discontinuity; eigenvalue-1={err:.6g}')
            lattice.log(kind='cold_root',n=n,nw=nw,lam=lam,omega=omega,T=root,
                        ceiling=ceiling,ratio=root/ceiling,eigenvalue_error=err)
        except Exception as exc:
            lattice.log(kind='cold_failure',n=n,nw=nw,lam=lam,omega=omega,error=str(exc))
