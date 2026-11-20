"""One-dimensional triangular-band crossed Eliashberg onset, eV units.

Continuum DOS N(E)=1/2 on [-1,1]; N(0)=1/2, A=1, E_F=1.
The 2pi-periodic dispersion 2|wrap(k)|/pi-1 is inversion symmetric;
its hopping is long-ranged, not a nearest-neighbour band. Finite grid has
zero-exterior Matsubara frequencies. This is not a certificate.
"""
import json, os, sys, time
import numpy as np
from scipy.fft import fft, ifft, next_fast_len
from scipy.optimize import brentq, newton_krylov, NoConvergence
from scipy.sparse.linalg import LinearOperator, eigs
from scipy.linalg import eigvalsh
from scipy.special import expit

KB=8.617333262145178e-5
C=0.18272624777232952
OUT='progress.jsonl'
DEADLINE=float(os.environ.get('G1_DEADLINE_EPOCH','inf'))
CURRENT_SPEC=''
def log(**x):
    with open(OUT,'a',encoding='utf8') as f: f.write(json.dumps(x,allow_nan=False)+'\n')
    print(json.dumps(x,allow_nan=False),flush=True)

def shift(a,ell):
    o=np.zeros_like(a)
    if ell>=0:
        if ell<len(a): o[:len(a)-ell]=a[ell:]
    else: o[-ell:]=a[:len(a)+ell]
    return o

def momentum(a):
    f=fft(a,axis=1); return f,f[:,(-np.arange(a.shape[1]))%a.shape[1]]

class FC:
    def __init__(self,v,L):
        self.L=L; self.n=next_fast_len(3*L-2); self.v=fft(v,self.n)[:,None]
    def __call__(self,a):
        return ifft(fft(a,self.n,axis=0)*self.v,axis=0)[self.L-1:2*self.L-1]

def crossed(a,b,c,v):
    L,n=a.shape; af,_=momentum(a); _,bm=momentum(b); cf,_=momentum(c); fc=FC(v,L)
    out=np.zeros_like(af)
    for ell in range(1-L,L): out+=v[ell+L-1]*fc(af*shift(bm,ell))*shift(cf,ell)
    return ifft(out,axis=1).real/n**2

def direct(a,b,c,v):
    L,n=a.shape; out=np.zeros_like(a)
    for m in range(L):
      for k in range(n):
       for mp in range(L):
        for mm in range(L):
         r=m+mm-mp
         if 0<=r<L:
          for x in range(n):
           for y in range(n):
            out[m,k]+=v[m-mp+L-1]*v[mp-mm+L-1]*a[mp,x]*b[mm,y]*c[r,(k+y-x)%n]/n**2
    return out

def audit():
    rng=np.random.default_rng(44); a,b,c=[rng.normal(size=(4,3)) for _ in range(3)]; v=np.array([.1,.2,.5,1.,.5,.2,.1])
    err=float(np.max(abs(crossed(a,b,c,v)-direct(a,b,c,v))))
    log(kind='audit',direct_fft_max_abs_error=err)
    if err>1e-11: raise RuntimeError('crossed audit failed')
    x,w,v,t=setup(3,2,1.,.2,100.)
    u=np.broadcast_to(w[:,None],(4,3)).copy()+rng.normal(scale=.01,size=(4,3))
    h=rng.normal(scale=.01,size=(4,3)); phi=rng.normal(scale=.01,size=(4,3))
    den=u*u+(x[None,:]+h)**2; z=-u/den; h0=-(x[None,:]+h)/den; p=-phi/den
    zn=meanconv(z,v)+t*(-direct(z,z,z,v)+direct(h0,z,h0,v)+direct(h0,h0,z,v)+direct(z,h0,h0,v))
    hn=meanconv(h0,v)+t*(-direct(h0,z,z,v)-direct(z,z,h0,v)-direct(z,h0,z,v)+direct(h0,h0,h0,v))
    pn=meanconv(p,v)+t*(-direct(p,z,z,v)-direct(z,z,p,v)+direct(p,h0,h0,v)+direct(h0,h0,p,v)+direct(z,p,z,v)+direct(h0,p,h0,v))
    un,hh=normal(u,h,x,w,v,t); pp=gap(u,h,x,v,t)(phi)
    errors=[float(np.max(abs(un-(w[:,None]-t*zn)))),float(np.max(abs(hh-t*hn))),float(np.max(abs(pp+t*pn)))]
    log(kind='factor_audit',normal_u_error=errors[0],normal_h_error=errors[1],gap_error=errors[2])
    if max(errors)>1e-10: raise RuntimeError('factored equation audit failed')

def setup(n,nw,lam,om,T):
    k=2*np.pi*np.arange(n)/n; x=2*np.abs((k+np.pi)%(2*np.pi)-np.pi)/np.pi-1
    t=KB*T; L=2*nw; w=np.pi*t*(2*np.arange(-nw,nw)+1); j=np.arange(1-L,L)
    v=2*lam*om**2/(om**2+(2*np.pi*t*j)**2)
    return x,w,v,t

def meanconv(a,v):
    from scipy.signal import fftconvolve
    return fftconvolve(a.mean(axis=1),v,mode='same')[:,None]

def normal(u,h,x,w,v,t):
    den=u*u+(x[None,:]+h)**2; z=-u/den; h=-(x[None,:]+h)/den
    zf,zm=momentum(z); hf,hm=momentum(h); L,n=z.shape; fc=FC(v,L)
    zo=np.zeros_like(zf); ho=np.zeros_like(hf)
    for ell in range(1-L,L):
        mz=shift(zm,ell); mh=shift(hm,ell)
        a=fc(-zf*mz+hf*mh); b=fc(hf*mz+zf*mh)
        dz=shift(zf,ell); dh=shift(hf,ell); ve=v[ell+L-1]
        zo+=ve*(a*dz+b*dh); ho+=ve*(-b*dz+a*dh)
    zn=meanconv(z,v)+t*ifft(zo,axis=1).real/n**2
    hn=meanconv(h,v)+t*ifft(ho,axis=1).real/n**2
    return w[:,None]-t*zn,t*hn

def gap(u,h,x,v,t):
    den=u*u+(x[None,:]+h)**2; z=-u/den; h=-(x[None,:]+h)/den
    zf,zm=momentum(z); hf,hm=momentum(h); L,n=z.shape; fc=FC(v,L)
    # Recompute shifted arrays to avoid O(L² n) persistent memory.
    fixed=[fc(-zf*shift(zm,e)+hf*shift(hm,e)) for e in range(1-L,L)]
    def apply(phi):
        p=-phi/den; pf,pm=momentum(p); out=np.zeros_like(pf)
        for i,e in enumerate(range(1-L,L)):
            mz=shift(zm,e); mh=shift(hm,e); mp=shift(pm,e)
            b=fc(-pf*mz+zf*mp); c=fc(pf*mh+hf*mp)
            out+=v[i]*(b*shift(zf,e)+c*shift(hf,e)+fixed[i]*shift(pf,e))
        return np.broadcast_to(-t*(meanconv(p,v)+t*ifft(out,axis=1).real/n**2),p.shape)
    return apply

def density(u,h,x,w,t,nw,delta):
    """Exact free reference at the *same* shifted chemical potential."""
    xx=x-delta
    ref=-xx[None,:]/(w[:,None]**2+xx[None,:]**2)
    green=-(xx[None,:]+h)/(u*u+(xx[None,:]+h)**2)
    return 2*np.mean(expit(-xx/t))+4*t*np.sum(np.mean(green[nw:]-ref[nw:],axis=1))

class Solver:
    # Successful states are available for continuation between adjacent lambda values.
    cache={}
    def __init__(self,n,nw,lam,om,vertex,fill=1.):
        self.n=n; self.nw=nw; self.lam=lam; self.om=om; self.vertex=vertex; self.fill=fill
        self.prev=None; self.history=[]
    def evaluate(self,T):
        x,w,v,t=setup(self.n,self.nw,self.lam,self.om,T); x=x+(1-self.fill); shape=(len(w),self.n)
        size=np.prod(shape)
        if self.fill==1. and self.n%2: raise ValueError('particle-hole symmetric grid requires even N')
        def project(vec):
            if self.fill!=1.: return vec
            uu=vec[:size].reshape(shape); hh=vec[size:2*size].reshape(shape)
            return np.concatenate((.5*(uu+np.roll(uu,self.n//2,axis=1)).ravel(),
                                   .5*(hh-np.roll(hh,self.n//2,axis=1)).ravel(),[0.]))
        candidates=[(0.,temp,state) for temp,state in self.history if temp>=T]
        if not candidates:
            for (n,nw,om,vertex,fill,lam),states in self.cache.items():
                if (n,nw,om,vertex,fill)!=(self.n,self.nw,self.om,self.vertex,self.fill) or lam==self.lam: continue
                candidates.extend((abs(lam-self.lam),temp,state) for temp,state in states if temp>=T)
        if candidates:
            _,_,(oldw,u,h,delta)=min(candidates,key=lambda a:(a[0],a[1]-T))
            u=u+(w-oldw)[:,None]; h=h.copy()
        else:
            u=np.broadcast_to(w[:,None],shape).copy(); h=np.zeros(shape); delta=0.
        def residual(vec):
            if time.time() >= DEADLINE:
                log(kind='deadline',spec=CURRENT_SPEC)
                raise SystemExit(124)
            uu=vec[:size].reshape(shape); hh=vec[size:2*size].reshape(shape); delta=vec[-1]
            xx=x-delta
            den=uu*uu+(xx[None,:]+hh)**2
            if self.vertex: un,hn=normal(uu,hh,xx,w,v,t)
            else:
                un=np.broadcast_to(w[:,None]+t*meanconv(uu/den,v),shape)
                hn=np.broadcast_to(-t*meanconv((xx[None,:]+hh)/den,v),shape)
            if self.fill==1.:
                dn=delta
            else:
                dn=self.fill-density(uu,hh,x,w,t,self.nw,delta)
            return np.concatenate(((un-uu).ravel(),(hn-hh).ravel(),[dn]))
        def scale(vec):
            return max(float(np.max(np.abs(vec))),self.om,1e-12)
        vec=project(np.concatenate((u.ravel(),h.ravel(),[delta])))
        best=vec.copy(); bestrel=float('inf'); bestabs=float('inf'); dx=[]; dr=[]
        last_vec=last_res=None; previous=float('inf'); method='mixing'
        for it in range(650):
            raw=residual(vec)
            abserr=float(np.max(np.abs(raw))); rel=abserr/scale(vec)
            rr=project(raw)
            if not np.isfinite(rel): raise RuntimeError(f'normal nonfinite; best_relative_residual={bestrel}')
            if rel<bestrel: best=vec.copy(); bestrel=rel; bestabs=abserr
            if rel<=1e-8: break
            if it<250:
                vec=vec+.2*rr
                continue
            method='anderson'
            if rel>1.25*previous:
                dx.clear(); dr.clear(); last_vec=last_res=None
            if last_vec is not None:
                dx.append(vec-last_vec); dr.append(rr-last_res)
                if len(dx)>5: dx.pop(0); dr.pop(0)
            candidate=vec+.2*rr
            if dr:
                dmat=np.column_stack(dr)
                try:
                    gamma=np.linalg.lstsq(dmat,rr,rcond=1e-10)[0]
                    trial=candidate-(np.column_stack(dx)+.2*dmat)@gamma
                    if np.all(np.isfinite(trial)) and scale(trial)<10*scale(vec): candidate=trial
                    else: dx.clear(); dr.clear()
                except np.linalg.LinAlgError:
                    dx.clear(); dr.clear()
            previous=rel; last_vec=vec; last_res=rr; vec=candidate
        else:
            method='newton_krylov'
            try:
                vec=newton_krylov(lambda z:project(residual(project(z))),best,
                                  f_tol=max(1e-15,scale(best)*1e-8),
                                  maxiter=30,inner_maxiter=12,line_search='armijo')
            except NoConvergence as exc:
                vec=np.asarray(exc.args[0])
            except (ValueError,RuntimeError) as exc:
                raise RuntimeError(f'normal failed; relative_residual={bestrel:.6g}, absolute_residual={bestabs:.6g}; newton_krylov={exc}') from exc
            vec=project(vec)
            rr=residual(vec); abserr=float(np.max(np.abs(rr))); rel=abserr/scale(vec)
            if not np.isfinite(rel) or rel>1e-8:
                raise RuntimeError(f'normal failed; relative_residual={rel:.6g}, absolute_residual={abserr:.6g}, best_relative_residual={bestrel:.6g}')
        u=vec[:size].reshape(shape); h=vec[size:2*size].reshape(shape); delta=vec[-1]; res=abserr
        self.prev=(w,u.copy(),h.copy())
        self.history.append((T,(*self.prev,delta))); self.history=self.history[-24:]
        key=(self.n,self.nw,self.om,self.vertex,self.fill,self.lam)
        self.cache[key]=self.history
        if self.vertex:
            apply=gap(u,h,x-delta,v,t); size=np.prod(shape)
            op=LinearOperator((size,size),matvec=lambda y:apply(y.reshape(shape)).ravel(),dtype=float)
            val=eigs(op,k=1,which='LR',tol=2e-5,maxiter=100,v0=np.ones(size),return_eigenvectors=False)[0]
            if abs(val.imag)>1e-7: raise RuntimeError(f'complex eigenvalue {val}')
            val=float(val.real)
        else:
            q=np.mean(1/(u[self.nw:]**2+(x[None,:]-delta+h[self.nw:])**2),axis=1)
            j=np.arange(self.nw); d=j[:,None]-j[None,:]; s=j[:,None]+j[None,:]+1
            f=lambda z: 2*self.lam*self.om**2/(self.om**2+(2*np.pi*t*z)**2)
            val=float(eigvalsh(t*np.sqrt(q[:,None]*q[None,:])*(f(d)+f(s)),subset_by_index=[self.nw-1,self.nw-1])[0])
        log(kind='evaluation',n=self.n,nw=self.nw,lam=self.lam,omega=self.om,fill=self.fill,mu_shift=float(delta),vertex=self.vertex,T=T,eigenvalue=val,residual=float(res),relative_residual=float(rel),iterations=it+1,normal_method=method)
        return val-1
    def solve(self):
        ceiling=C*self.om*np.sqrt(self.lam)/KB
        hi=ceiling*1.8
        fhi=self.evaluate(hi)
        if fhi>=0: raise RuntimeError('upper bound already unstable')
        for _ in range(30):
            lo=hi/1.25; flo=self.evaluate(lo)
            if flo>0: break
            hi=lo
        else: raise RuntimeError('no unstable point found in temperature scan')
        root=brentq(self.evaluate,lo,hi,xtol=.005,rtol=1e-5)
        froot=self.evaluate(root)
        if abs(froot)>1e-4:
            raise RuntimeError(f'root may be branch discontinuity; eigenvalue-1={froot:.6g}')
        log(kind='root',n=self.n,nw=self.nw,lam=self.lam,omega=self.om,fill=self.fill,vertex=self.vertex,T=root,ceiling=ceiling,ratio=root/ceiling)

def main():
    if sys.argv[1]=='audit': audit(); return
    failed=False
    for s in sys.argv[1:]:
        global CURRENT_SPEC
        CURRENT_SPEC=s
        n,nw,lam,om,vertex,*extra=s.split(':')
        fill=float(extra[0]) if extra else 1.
        try: Solver(int(n),int(nw),float(lam),float(om),bool(int(vertex)),fill).solve()
        except Exception as e:
            failed=True
            log(kind='error',spec=s,error=repr(e))
    if failed: sys.exit(1)
if __name__=='__main__': main()
