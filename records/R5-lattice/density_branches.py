"""Compare cold and high-temperature continuation at one fixed T."""
import sys
import json
import numpy as np
import lattice

if len(sys.argv)==6:
    n,nw=int(sys.argv[1]),int(sys.argv[2])
    lam,omega,T=map(float,sys.argv[3:6])
else:
    n,nw,lam,omega,T=(66,8,.5,.03,float(sys.argv[1]) if len(sys.argv)>1 else 34.749)
states=[]
for label in ('cold','warm'):
    lattice.Solver.cache.clear()
    solver=lattice.Solver(n,nw,lam,omega,True,.6)
    try:
        if label=='warm': solver.evaluate(T+20)
        value=solver.evaluate(T)+1
        w,u,h=solver.prev
        delta=solver.history[-1][1][3]
        states.append((u,h,delta))
        print(json.dumps({'branch':label,'T':T,'mu_shift':delta,'eigenvalue':value}),flush=True)
    except Exception as exc:
        print(json.dumps({'branch':label,'T':T,'failure':str(exc)}),flush=True)
if len(states)==2:
    (u,h,d),(u2,h2,d2)=states
    print(json.dumps({'max_u_difference':float(np.max(abs(u-u2))),
                      'max_h_difference':float(np.max(abs(h-h2))),
                      'mu_difference':abs(d-d2)}),flush=True)
