"""Dyadic momentum-dependent normal centers and continuous Krawczyk bounds."""
import argparse
import os
from fractions import Fraction as Q

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from scipy.fft import dst, idst
from scipy.linalg import solve

from finite import PUBLISHED
from radial_space import Space, down, up
from replay import C_BOUND, ETA, GRID, LAM, VOLUME, box, choose, decimal, iv, require

V = float(VOLUME)


def ends(a):
    return float(np.nextafter(float(a.a),-np.inf)),float(np.nextafter(float(a.b),np.inf))


def fold(j):
    return np.where(j >= 0,j,-1-j)


def dyadic(a):
    return np.ldexp(np.rint(np.ldexp(a.real,36)),-36)+1j*np.ldexp(np.rint(np.ldexp(a.imag,36)),-36)


class Trial:
    """The sine transform is used only for proposing fully momentum-dependent centers."""
    def __init__(self, x, n, cells, numerator, alpha, panel, spatial=None):
        self.n,self.cells = n,cells
        self.spatial = n if spatial is None else spatial
        self.r = (np.arange(cells)+.5)/cells
        self.weight = 3*self.r**2/cells
        self.freq = 2*np.arange(n+2*panel)+1
        self.bare = 1j*self.freq[:,None]-(2*self.r**3-1)/x
        self.q = np.pi/4*np.arange(1,4*cells+1)
        self.y = np.arange(self.spatial)[:,None]-np.arange(-panel,panel+1)[None,:]
        self.b = self.y[:,:,None]+self.y[:,None,:]-np.arange(self.spatial)[:,None,None]
        t = np.arange(-panel,panel+1)
        self.ell = numerator/(1+alpha*t**2)
        # The omitted cubic wings are bounded by the verifier. Retain a longer
        # signed rainbow prefix in the proposal to avoid a gratuitous defect.
        s = np.arange(-max(160,8*n),max(160,8*n)+1)
        self.rain = numerator/(1+alpha*s*s)
        self.ridx = np.arange(n)[:,None]-s
        f = fold(self.ridx)
        self.rfold = np.minimum(f,n-1)
        self.rbare = np.sign(-self.ridx-.5)*1j*x*np.arctan(1/(x*(2*f+1)))
        self.rhigh = f >= n

    def propose(self, theta, old):
        n = self.n
        M = len(self.freq)
        s = np.zeros((M,self.cells),dtype=complex)
        s[:n] = old
        for _ in range(150):
            z = 1/(self.bare-s)
            mean = (z[:n]@self.weight)[self.rfold]
            mean = np.where(self.ridx>=0,mean,mean.conj())
            mean = np.where(self.rhigh,self.rbare,mean)
            hats = dst(z*self.r,type=2,n=len(self.q),axis=1)*((1.5/self.cells)/self.q)
            hats = np.r_[hats[::-1].conj(),hats]
            a = hats[self.y+M]*self.ell[None,:,None]
            c = np.asarray([np.einsum('it,jt,ijt->t',aa,aa,hats[bb+M]) for aa,bb in zip(a,self.b)])
            update = np.repeat((theta*(mean@self.rain))[:,None],self.cells,axis=1)
            update[:self.spatial] += theta**2*idst(c*(self.q/(1.5/self.cells)),type=2,axis=1)[:,:self.cells]/self.r
            delta = (update-s[:n])
            if np.max(abs(delta)/(self.freq[:n,None])) < 2e-11:
                return dyadic(s[:n])
            s[:n] += .65*delta
        raise ArithmeticError('normal proposal did not converge')


class Replay:
    def __init__(self, lam, eta, x, n, cells, panel, space=Space):
        self.n,self.cells,self.j = n,cells,panel
        self.spatial = n
        self.cut = max(4*n,72)
        self.outer = max(2*self.cut,128)
        self.count = self.cut+self.outer+2
        self.w = 2*np.arange(self.count)+1
        self.rows = self.w[:self.cut]
        self.x = float(x)
        self.space = space(self.x,cells)
        self.a0,self.a1 = ends(2*box(lam/x)/iv.pi)
        self.b0,self.b1 = ends((2*iv.pi*C_BOUND)**2*box(lam*eta**2))
        indices = np.arange(-self.outer,self.outer+1)
        self.idx = fold(np.arange(self.cut)[:,None]-indices[None,:])
        self.ell_hi = up(self.a1/(1+self.b0*indices**2))
        self.ell = np.maximum(0.,down(self.a0/(1+self.b1*indices**2)))
        self.ell_err = up(self.ell_hi-self.ell,self.ell_hi)
        self.prefix = np.arange(n)[:,None]-indices[None,:]
        self.local = abs(indices)<=panel
        self.missing = up(2*self.a1/(self.b0*self.outer))
        self.mass = up(sum(self.ell_hi)+self.missing)
        # Evaluate expensive spatial cells only where crossed signs matter most.
        self.signed = min(n,4)

    def g(self,j):
        return self.ell[self.outer+j]

    def ge(self,j):
        return self.ell_err[self.outer+j]

    def sums(self,p,cap):
        """Each profile is <= cap/|2m+1| (weighted where appropriate)."""
        far = self.missing*cap/(2*(self.outer+1)-(2*self.cut-1))
        a = up(p[self.idx[:,self.local]]@self.ell_hi[self.local])
        b = up(p[self.idx[:,~self.local]]@self.ell_hi[~self.local]+far)
        # n>=cut and |2m+1|<cut forces |n-m|>=cut/2.
        high = up(self.mass*cap/self.cut+2*self.ell_hi[self.outer+self.cut//2]*sum(p[self.w<self.cut]))
        return np.r_[a,0.],np.r_[b,high]

    @staticmethod
    def wings(a,b):
        return a[1]*(b[0]+b[1])+a[0]*b[1]

    def envelopes(self,line):
        w = self.w
        pw = np.max(line.peak,axis=1)
        L = (line.norm1,line.norm2,pw)
        sc = up(1/line.width)
        S = [self.sums(p,sc) for p in L]
        maxZ = [max(p) for p in L]
        out = []
        for weighted in (True,False):
            mul = w if weighted else 1.
            dL = (mul*line.p2,mul*pw*line.norm2,mul*pw**2)
            dS = [self.sums(p,sc**2) for p in dL]
            out.append((dL,dS))
        def cross(C, D, F, mx, wing):
            p = self.wings if wing else lambda a,b: (a[0]+a[1])*(b[0]+b[1])
            return np.minimum.reduce((mx[0]*p(C[1],F[1]),mx[1]*p(C[0],F[1]),
                                      mx[1]*p(C[1],F[0]),V*mx[2]*p(C[2],F[2])))
        envs = []
        for wing in (False,True):
            f = up(cross(S,S,S,maxZ,wing))
            derivative = []
            for dL,dS in out:
                end = cross(dS,S,S,maxZ,wing)
                mid = cross(S,dS,S,[max(p) for p in dL],wing)
                derivative.append(up(2*end+mid))
            envs.append((f,*derivative))
        return up(S[0][0]+S[0][1]),tuple(up(d[1][0][0]+d[1][0][1]) for d in out),envs

    def pair_sum(self,line,n,m):
        a = np.zeros(2*self.cells,dtype=complex)
        b = np.zeros_like(a)
        err = np.zeros(2*self.cells)
        l2 = scale = 0.
        for t in range(-self.j,self.j+1):
            h,c,e,s = line.pair(n-t,m-t)
            g,ge = self.g(t),self.ge(t)
            a += g*h
            b += g*c
            err += g*e+ge*(abs(c)+e)
            l2 += g*s+ge*(abs(s)+self.space.norm(h))
            scale += (g+ge)*(1+max(abs(c))+max(abs(h)))
        return a,b,up(err,scale),up(l2,scale)

    def cubic_error(self,line,n,m,h,e):
        f = int(fold(m))
        return up(e*line.norm2[f]+self.space.norm(h)*line.e2[f])

    def residual(self,line,center,theta,rback,fw):
        f = self.idx[:self.n]
        z = line.mean[f]
        z = np.where(self.prefix>=0,z,z.conj())
        means = z@self.ell
        errs = up(line.mean_err[f]@self.ell_hi+abs(z)@self.ell_err
                  +self.missing/(line.width*(2*self.outer-2*self.cut+3)),abs(z)@self.ell_hi)
        v = []
        for n in range(self.n):
            mean,er_mean = means[n],errs[n]
            scale = 1.+abs(z[n])@self.ell_hi
            if n>=self.spatial:
                v.append(up(max(abs(theta*mean-center[n]))+theta*er_mean+theta**2*self.f_abs[n],
                            scale+max(abs(center[n])))/self.rows[n])
                continue
            e = 0.
            out = np.zeros(self.cells,dtype=complex)
            var = out.real.copy()
            for s in range(-self.j,self.j+1):
                m = n-s
                f = int(fold(m))
                ph,_,_,error = self.pair_sum(line,n,m)
                _,cv,rv,_ = self.space.convolve(ph,line.value(m))
                g,ge = self.g(s),self.ge(s)
                err = self.cubic_error(line,n,m,ph,error)
                out += g*cv[:self.cells]
                var += g*rv[:self.cells]+ge*(abs(cv[:self.cells])+rv[:self.cells])
                e += (g+ge)*err
                scale += (g+ge)*(1+max(abs(cv)))
            v.append(up(max(abs(theta*mean+theta**2*out-center[n])+theta**2*var)
                        +theta*er_mean+theta**2*(e+fw[n]),scale+max(abs(center[n])))/self.rows[n])
        high = self.high(theta*rback+theta**2*self.f_abs)
        return np.r_[v,high]

    def high(self,v,weighted=True):
        if not weighted:
            return max(v[self.n:])
        return max(np.r_[v[self.n:-1]/self.rows[self.n:],v[-1]/(2*self.cut+1)])

    def signed_rows(self,line,theta,rb,rb_pair,xb,tail):
        raw_normal = up(xb[0][1]/np.r_[self.rows,2*self.cut+1]*theta**2)
        raw_pair = up(theta*rb_pair+theta**2*xb[0][2])
        for n in range(self.signed):
            normal,pair = {},{}
            def add(dest,m,c,e,g,ge):
                if m not in dest:
                    dest[m] = (np.zeros(2*self.cells,dtype=c.dtype),np.zeros(2*self.cells),0.)
                val,err,scale = dest[m]
                val += g*c
                err += g*e+ge*(abs(c)+e)
                dest[m] = (val,err,scale+(g+ge)*max(abs(c)+e))
            for s in range(-self.j,self.j+1):
                m = n-s
                _,c,e,_ = self.pair_sum(line,n,m)
                a = 2*theta**2*self.g(s)
                ae = 2*theta**2*self.ge(s)
                add(normal,m,c,e,a,ae)
                add(pair,m,c.real,e,a,ae)
                for t in range(-self.j,self.j+1):
                    k = n-s-t
                    _,c,e,_ = line.pair(n-s,n-t)
                    g = theta**2*self.g(s)*self.g(t)
                    ge = up(theta**2*((self.g(s)+self.ge(s))*(self.g(t)+self.ge(t)))-g)
                    add(normal,k,c,e,g,ge)
                    _,c,e,_ = line.pair(n-s,-1-(n-t))
                    add(pair,k,c.real,e,g,ge)
            values = []
            for is_normal,terms in ((True,normal),(False,pair)):
                row = np.zeros(self.cells)
                scale = 1.
                for m,(c,e,sc) in terms.items():
                    f = int(fold(m))
                    if not is_normal and abs(n-m)<=self.outer:
                        c = c + theta*self.g(n-m)
                        e = e + theta*self.ge(n-m)
                    peak = up(line.peak[f]**2)
                    K = up(abs(c)+e,sc+max(abs(c)))
                    a = (np.arange(2*self.cells)+1)/self.cells*K
                    b = self.space.right*peak
                    _,cv,rv,_ = self.space.convolve(a,b)
                    coef = (self.w[f]/self.w[n]) if is_normal else 1.
                    row += coef*(cv[:self.cells].real+rv[:self.cells])
                    scale += coef*max(abs(cv)+rv)
                missing = 0.
                if not is_normal:
                    # Remove only explicitly accounted rainbow terms from the positive sum.
                    # Sum the complement directly, so this bound does not subtract intervals.
                    for s in range(-self.outer,self.outer+1):
                        if n-s not in terms:
                            missing += self.ell_hi[self.outer+s]*line.p2[int(fold(n-s))]
                    missing += self.missing/line.width**2/(2*(self.outer+1)-(2*self.cut-1))**2
                rem = theta**2*(tail[1][n]/self.w[n] if is_normal else tail[2][n])
                values.append(up(max(row)+rem+theta*missing,scale+abs(missing)))
            raw_normal[n] = min(raw_normal[n],values[0])
            raw_pair[n] = min(raw_pair[n],values[1])
        return np.r_[raw_normal[:self.n],self.high(raw_normal,False)],max(raw_pair)

    def inverse(self,line,theta,mutant=None):
        """Return bounds on A and on I-A(I-R), not a trusted computed inverse."""
        n = self.n
        k = np.arange(n)
        W = self.w[:n]
        f = (k[:,None]-k[None,:])
        b = (k[:,None]+k[None,:]+1)
        g = theta*self.ell[self.outer+f]*W[None,:]/W[:,None]
        h = theta*self.ell[self.outer+b]*W[None,:]/W[:,None]
        J = np.zeros((2*n,2*n))
        J[::2,::2] = g+h
        J[1::2,1::2] = g-h
        E = up(theta*(self.ell_err[self.outer+f]+self.ell_err[self.outer+b])*W[None,:]/W[:,None])
        J_err = up(2*E)
        D = np.zeros_like(J)
        d = line.sqavg[:n]
        D[2*k,2*k] = D[2*k+1,2*k+1] = d.real
        D[2*k,2*k+1] = -d.imag
        D[2*k+1,2*k] = d.imag
        # Block matrices act on complexes with the Euclidean norm in each slot.
        def blocks(X):
            return X.reshape(n,2,n,2).transpose(0,2,1,3)
        def norm(X):
            B = abs(blocks(X))
            frobenius = up(np.sqrt(np.sum(B**2,axis=(2,3))))
            columns = up(np.max(np.sum(B,axis=2),axis=2))
            rows = up(np.max(np.sum(B,axis=3),axis=2))
            return np.minimum(frobenius,up(np.sqrt(columns*rows)))
        M = np.round(solve((np.eye(2*n)-D@J).T,J.T).T*2**36)*2**-36
        if mutant == 'inverse':
            M[0,0] += 32
        bM = norm(M)
        def apply(b):
            return up(np.r_[b[:-1]+bM@(line.p2[:n]*b[:-1]),b[-1]])
        d_err = line.sqavg_err[:n]
        j_bound = up(np.maximum(abs(g+h),abs(g-h))+J_err)
        dj_err = up(d_err[:,None]*j_bound+up(abs(d))[:,None]*J_err)
        defect = up(norm(J-M+M@D@J)+J_err+bM@dj_err,1+bM@((1+abs(d))[:,None]*j_bound))
        # Independently bound the actual rounded solve; zeroing or changing it cannot pass this test by fiat.
        krawczyk = up(np.r_[defect@line.p2[:n],0.])
        return apply,krawczyk

    def anchor(self,center,theta):
        center_all = np.zeros((self.count,self.cells),dtype=complex)
        center_all[:self.n] = center
        line = self.space.state(center_all,self.w)
        rb,(rN,rP),(core,wing) = self.envelopes(line)
        self.f_abs = core[0]
        res = self.residual(line,center,theta,rb,wing[0])
        return line,rb,rN,rP,core,wing,res

    def finish(self,record,right,theta0,theta1,delta,mutant=None):
        line,rb,rN,rP,core,wing,res = record
        # All ratios below use the weighted W_n=|2n+1| norm.
        w = np.r_[self.rows,2*self.cut+1]
        summary = lambda arr: np.r_[arr[:self.n],max(arr[self.n:])]
        xD = summary(core[1]/w)
        rD = summary(rN/w)
        if hasattr(line,'growth'):
            growth_rainbow,growth_crossed,_ = line.growth
            rV,xV = summary(growth_rainbow/w),summary(growth_crossed/w)
        else:
            rV,xV = rD,xD
        t = theta1-theta0
        sxD,_ = self.signed_rows(line,1.,rb,rP,(core,wing),wing)
        pre,inv = self.inverse(line,theta0,mutant)
        # Omitted rainbow part of DF_c after removing the finite rank J V.
        # Head rows: only lines m>=N; tail rows: all lines.
        head = []
        profile = line.p2*self.w
        for i in range(self.n):
            row = self.idx[i]
            exists = row>=self.n
            head.append(up(profile[row[exists]]@self.ell_hi[exists]
                           + self.missing/(line.width**2*(2*(self.outer+1)-(2*self.cut-1))))/w[i])
        hole = up(np.r_[head,rD[-1]]*theta0)
        midpoint = up(pre(theta1**2*sxD+hole+t*rD)+inv)
        if max(inv) >= 1:
            return None,'inverse defect >= 1'
        d = up(delta/line.width)
        if d >= .20:
            return None,'anchor shift > 0.20'
        rF,xF = summary(rb/w),summary(core[0]/w)
        # Affine centers: both endpoint residuals plus a second derivative allowance.
        curv = (2*t*d+2*theta1*d*d)/(1-d)**3*rF
        curv += (2*t*t+12*theta1*t*d+12*theta1**2*d*d)/(1-d)**5*xF
        eta = up(pre(np.maximum(res,right[-1])+up(curv/8)))
        coarse = max(eta)
        floor = np.inf
        sparse = (1.06,1.2,1.4,1.65,2.,2.5,3.,3.8,4.8,6.,8.,12.)
        # Expanded candidates are proposals only; every radius faces the same strict test.
        for factor in (*sparse, *(1+i/16 for i in range(1,321))):
            rho = up(factor*coarse)
            e = up((delta+rho)/line.width)
            if e >= .30:
                continue
            # Keep signed insertion coefficients when the perturbed line is
            # the inserted one. Only variation of the two inner lines uses xD.
            k = up(e*(2-e)/(1-e)**2)
            growth = k*(theta1*rV+theta1**2*(sxD+xV/(1-e)**2))
            q = up(midpoint+pre(up(growth)))
            floor = min(floor,max(q))
            if max(q)<1 and max(eta+rho*q)<down(rho):
                z = up(max(pre(right[-1]))/(1-max(q)))
                return (rho,max(q),z),''
        return None,f'normal center defect={coarse:.6g}, minimum guard q={floor:.6g}'

    def pairing(self,record,rho):
        line,rb,rN,rP,core,wing,res = record
        _,p = self.signed_rows(line,1.,rb,rP,(core,wing),wing)
        e = up(rho/line.width)
        if e >= .3:
            return np.inf
        k = up(e*(2-e)/(1-e)**2)
        absolute_growth = line.growth[2] if hasattr(line,'growth') else core[2]
        return up(p*(1+k)+k*max(absolute_growth)/(1-e)**2)


def open_nodes():
    nodes = [(lam,eta,x) for x in GRID for lam in reversed(LAM) for eta in ETA[:2]
             if choose(lam,eta,x) is None and x not in PUBLISHED[ETA.index(eta)][LAM.index(lam)]]
    require(len(nodes) == 44,'phase 2 open grid changed')
    return nodes


def run(primary, trial, steps, mutant=None):
    center = np.zeros((primary.n,primary.cells),dtype=complex)
    records = [primary.anchor(center,0.)]
    witness = [center]
    qmax = 0.
    for i in range(steps):
        lo,hi = i/steps,(i+1)/steps
        new = trial.propose(hi,center)
        r = primary.anchor(new,hi)
        delta = up(np.max(abs(new-center)/primary.w[:primary.n,None]))
        ans,why = primary.finish(records[-1],r,lo,hi,delta,mutant)
        if ans is None:
            return None,f'primary slab={i+1}/{steps}: {why}'
        _,q,z = ans
        qmax = max(qmax,q)
        center = new
        records[-1] = r
        witness.append(center)
    pair = primary.pairing(records[-1],z)
    if not pair < 1:
        return None,f'primary q={qmax:.6g}; final pairing={pair:.6g}'
    return (pair,qmax,witness),' '


def main(proto=Replay, verify=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells',default=128,type=int)
    p.add_argument('--head',default=16,type=int)
    p.add_argument('--panel',default=8,type=int)
    p.add_argument('--steps',default=24,type=int)
    p.add_argument('--lam',type=Q,nargs='+',default=LAM)
    p.add_argument('--eta',type=Q,nargs='+',default=ETA[:2])
    p.add_argument('--xs',type=Q,nargs='+',default=GRID)
    p.add_argument('--node',action='append',metavar='LAMBDA,ETA,X',
                   help='select an exact node explicitly; repeat for a set')
    p.add_argument('--probe',action='store_true',help='do not issue a certificate')
    args = p.parse_args()
    require(32<=args.cells<=1024 and args.cells&(args.cells-1)==0 and 4<=args.head<=64
            and 4<=args.panel<=20 and 2<=args.steps<=256,'dimensions outside arithmetic bounds')
    grid = open_nodes()
    if args.node:
        selected = [tuple(Q(x) for x in node.split(',')) for node in args.node]
        require(len(set(selected))==len(selected) and all(node in grid for node in selected),
                'duplicate or non-open exact node')
        nodes = [node for node in grid if node in selected]
    else:
        nodes = [(l,e,x) for l,e,x in grid if l in args.lam and e in args.eta and x in args.xs]
    require(nodes,'no open nodes requested')
    bad = 0
    for lam,eta,x in nodes:
        print(f'ATTEMPT {lam},{eta},{x}',flush=True)
        r = proto(lam,eta,x,args.head,args.cells,args.panel)
        trial = Trial(float(x),r.n,r.cells,r.a0,r.b1,r.j,r.spatial)
        try:
            z,why = run(r,trial,args.steps)
        except (ArithmeticError, ValueError) as exc:
            z,why = None,f'primary arithmetic enclosure: {exc}'
        if z is None:
            bad += 1
            print('UNRESOLVED',why,flush=True)
            continue
        a,q,witness = z
        if args.probe:
            print(f'PRIMARY pairing<={decimal(a,8)} q<={decimal(q,8)}',flush=True)
            continue
        if verify is None:
            from radial_check import verify
        try:
            out,why = verify(lam,eta,x,r,witness)
        except (ArithmeticError, ValueError) as exc:
            out,why = None,f'second arithmetic enclosure: {exc}'
        if out is None:
            bad += 1
            print('UNRESOLVED',why,flush=True)
            continue
        print(f'CLOSED {lam},{eta},{x} pairing<={decimal(max(a,out[0]),8)} '
              f'q_G<={decimal(max(q,out[1]),8)}',flush=True)
    if bad:
        raise SystemExit(f'{bad}/{len(nodes)} requested Krawczyk nodes unresolved')
    if not args.probe:
        print(f'RADIAL PASS closed={len(nodes)}',flush=True)


if __name__ == '__main__':
    main()
