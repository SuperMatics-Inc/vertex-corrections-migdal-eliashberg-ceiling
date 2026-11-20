"""Trace the normal state reached at 17.977 K to its pairing crossing."""
import numpy as np
from scipy.optimize import brentq
import lattice

s=lattice.Solver(66,8,.5,.03,True,.6)
try:
    hi=17.977495436676854
    fhi=s.evaluate(hi)
    if fhi>=0: raise RuntimeError(f'initial eigenvalue-1={fhi}')
    for lo in np.arange(hi-.1,16.,-.1):
        flo=s.evaluate(float(lo))
        if flo>0: break
        hi=float(lo)
    else: raise RuntimeError('no sign change')
    root=brentq(s.evaluate,float(lo),hi,xtol=.001,rtol=1e-6)
    residual=s.evaluate(root)
    if abs(residual)>1e-4: raise RuntimeError(f'branch-switch eigenvalue-1={residual}')
    ceiling=lattice.C*.03*np.sqrt(.5)/lattice.KB
    lattice.log(kind='secondary_root',T=root,ceiling=ceiling,ratio=root/ceiling,eigenvalue_error=residual)
except Exception as exc:
    lattice.log(kind='secondary_failure',error=str(exc))
