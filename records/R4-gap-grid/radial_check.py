"""Second radial replay: half-line integer products and analytic DOS moments."""

from radial import Replay, fold, up, down
import numpy as np
from radial_space import Space
from replay import require


class HalfLine(Space):
    reverse_cross = False

    def convolve(self,a,b):
        """Exact integer folded correlations, with interval input quantization."""
        m,n = len(a),len(b)
        umax = max(1.,max(abs(a)),max(abs(b)))
        power = min(24,int(np.floor(np.log2(2.**25/umax))))
        require(power >= 12,'independent fixed point scale')
        unit = 2.**-power
        def encode(x):
            return (np.rint(x.real/unit).astype(np.int64),np.rint(x.imag/unit).astype(np.int64))
        ar,ai = encode(a)
        br,bi = encode(b)
        aa = max(int(max(abs(ar))),int(max(abs(ai))))
        bb = max(int(max(abs(br))),int(max(abs(bi))))
        require(4*min(m,n)*aa*bb < 2**62,'independent exact convolution overflow')
        def correlate(x,y):
            xy = np.zeros(m+n+1,dtype=np.int64)
            xy[1:-1] = np.convolve(x,y)
            sign = 1 if self.reverse_cross else -1
            xy[:m] += sign*np.correlate(x,y,'full')[n-1:]
            xy[:n] += sign*np.correlate(y,x,'full')[m-1:]
            return xy
        re = correlate(ar,br)-correlate(ai,bi)
        im = correlate(ar,bi)+correlate(ai,br)
        sd = -(3*self.h/2)*unit**2*(re.astype(float)+1j*im.astype(float))
        rg = np.arange(m+n)+.5
        # Half-cell integration of the exact affine derivative; total omitted
        # derivative from quantization evaluated before any correlation cancellation.
        p = sd[:-1]
        q = sd[1:]
        tail = np.r_[0j,np.cumsum((p+q)/2)[:-1]]
        central = (tail+(3*p+q)/8)/rg
        hd = self.h*(tail+(2*p+q)/6) if self.mean_height else self.h*rg*central
        # Any average of d is at most a left-interval comparison plus the error.
        baseline = np.r_[sd[0],tail[1:]/(rg[1:]-.5)]
        d0 = np.maximum(abs(p-baseline),abs(q-baseline))
        sc = up(3*self.h*min(m,n)*max(abs(a))*max(abs(b)))
        quant = up(3*self.h*min(m,n)*unit*(max(abs(a))+max(abs(b))+unit),sc)
        e = up(d0/(rg+.5)+abs(baseline-central)+4*quant)
        der = up(np.maximum(abs(p),abs(q))+quant,sc)
        el2 = up(self.h/np.sqrt(12)*self.norm(der)+(m+n)*self.h*np.sqrt(3*(m+n)*self.h)*4*quant)
        return hd,central,e,el2

    def state(self,sigma,widths):
        result = super().state(sigma,widths)
        # Exact flat-DOS primitive inside each radial step. Integrate the even
        # geometric series with a certified remainder, not radial midpoint quadrature.
        v = 1j*widths[:,None]-sigma
        y = (2*np.r_[0.,self.right]**3-1)/self.x
        half = np.diff(y)/2
        mass = self.x*half
        d = v-(y[1:]+y[:-1])/2
        mid = 1/d
        q = (half*mid)**2
        qu = up(abs(q))
        require(np.max(qu) < .30,'independent DOS series radius')
        m = 15
        poly = np.zeros_like(q)
        for k in range(m-1,-1,-1):
            poly = 1/(2*k+1)+q*poly
        val1 = mass*mid*poly
        b1 = up(mass*up(abs(mid))*qu**m/(1-qu),mass*up(abs(mid))/(1-qu))
        val2 = mass/((d-half)*(d+half))
        b2 = up(0.,mass*result.peak**2*up(1+abs(v)+abs(y[1:])+abs(y[:-1])))
        # complex absolute error before the signed frequency summation.
        result.mean = np.sum(val1,axis=1)
        result.mean_err = up(np.sum(b1,axis=1),np.sum(abs(val1),axis=1))
        result.sqavg = np.sum(val2,axis=1)
        result.sqavg_err = up(np.sum(b2,axis=1),np.sum(abs(val2),axis=1))
        return result


class Independent(Replay):
    omit_wing = False

    def __init__(self,lam,eta,x,n,cells,panel):
        super().__init__(lam,eta,x,n,cells,panel,HalfLine)
        self.m = np.arange(self.count)
        left = np.arange(self.cut)[:,None]-self.m
        right = np.arange(self.cut)[:,None]+self.m+1
        ia = np.minimum(abs(left),self.outer)
        ib = np.minimum(right,self.outer)
        self.lower = (self.ell[self.outer+ia]*(abs(left)<=self.outer),
                      self.ell[self.outer+ib]*(right<=self.outer))
        self.raise_to = (self.ell_hi[self.outer+ia]*(abs(left)<=self.outer),
                         self.ell_hi[self.outer+ib]*(right<=self.outer))
        self.central = (abs(left)<=panel,right<=panel)

    def transfer_wing(self,cap):
        return 0. if self.omit_wing else up(self.missing*cap/(2*self.outer-2*self.cut+3))

    def sums(self,p,cap):
        # Group by electron line rather than a contiguous transfer sum. Verify all
        # transfer weights outside the main panel, including both reflected tails.
        mass = self.raise_to[0]+self.raise_to[1]
        center = self.raise_to[0]*self.central[0]+self.raise_to[1]*self.central[1]
        rest = self.raise_to[0]*(~self.central[0])+self.raise_to[1]*(~self.central[1])
        close = up(center@p)
        far = up(rest@p+self.transfer_wing(cap))
        # A different electron split for the no-center high rows.
        width = self.cut+self.cut//4
        sep = self.cut-width//2
        high = up(self.mass*cap/width+2*self.ell_hi[self.outer+sep]*sum(p[self.w<width]))
        return np.r_[close,0.],np.r_[far,high]

    def residual(self,line,center,theta,rback,fw):
        d = np.r_[self.rows,2*self.cut+1]
        pos,neg = (a[:self.n] for a in self.lower)
        posh,negh = (a[:self.n] for a in self.raise_to)
        r = (pos+neg)@line.mean.real+1j*((pos-neg)@line.mean.imag)
        scale_r = (posh+negh)@abs(line.mean)
        er = up((posh+negh)@line.mean_err+(up(posh-pos)+up(negh-neg))@abs(line.mean)
                +self.missing/(line.width*(2*self.outer-2*self.cut+3)),scale_r)
        # Rebuild the crossed term from the ordered double-transfer prefix. Group
        # pair heights by the line hit by the outer convolution, not by derivative.
        head = np.empty(self.n)
        for n in range(self.n):
            if n>=self.spatial:
                head[n] = up(max(abs(theta*r[n]-center[n]))+theta*er[n]+theta**2*self.f_abs[n],
                             scale_r[n]+max(abs(center[n])))/d[n]
                continue
            v = np.zeros(self.cells,dtype=complex)
            e = np.zeros(self.cells)
            scale = 1.
            for m in range(n-self.j,n+self.j+1):
                items = [line.pair(k,m+k-n) for k in range(n-self.j,n+self.j+1)]
                ell = np.array([self.g(n-k) for k in range(n-self.j,n+self.j+1)])
                rel = np.array([self.ge(n-k) for k in range(n-self.j,n+self.j+1)])
                h = np.stack([it[0] for it in items])
                lp = np.array([it[3] for it in items])
                p = ell@h
                lp = up((ell+rel)@lp+rel@self.space.norm(h), (ell+rel)@self.space.norm(h))
                _,s,rv,_ = self.space.convolve(line.value(m),p)
                g,ge = self.g(n-m),self.ge(n-m)
                val = s[:self.cells]
                err = rv[:self.cells]+self.cubic_error(line,n,m,p,lp)
                v += g*val
                e += (g+ge)*err+ge*abs(val)
                scale += (g+ge)*max(abs(val)+err)
            head[n] = up(max(abs(theta*r[n]+theta**2*v-center[n])+theta**2*e)
                         +theta*er[n]+theta**2*fw[n],scale+scale_r[n]+max(abs(center[n])))/d[n]
        return np.r_[head,self.high(theta*rback+theta**2*self.f_abs)]

    def signed_rows(self,line,theta,rb,rb_pair,xb,tail):
        # Group by output electron index, not by the omitted transfer.
        W = np.r_[self.rows,2*self.cut+1]
        normal = up(theta**2*xb[0][1]/W)
        pair = up(theta*rb_pair+theta**2*xb[0][2])
        grid = self.space.h*np.arange(1,2*self.cells+1)
        for n in range(self.signed):
            ks = range(n-self.j,n+self.j+1)
            acc_n = np.zeros(self.cells)
            acc_p = np.zeros(self.cells)
            sc_n = sc_p = 0.
            for m in range(n-2*self.j,n+2*self.j+1):
                N = np.zeros(2*self.cells,complex)
                P = np.full(2*self.cells,theta*self.g(n-m))
                en = np.zeros(self.cells*2)
                ep = np.full(2*self.cells,theta*self.ge(n-m))
                scale = 1.+theta*self.ell_hi[self.outer+n-m]
                if m in ks:
                    for k in ks:
                        _,v,e,_ = line.pair(k,m+k-n)
                        a = 2*theta**2*self.g(n-m)*self.g(n-k)
                        ah = up(2*theta**2*self.ell_hi[self.outer+n-m]*self.ell_hi[self.outer+n-k])
                        N += a*v
                        P += a*v.real
                        en += ah*e+up(ah-a)*abs(v)
                        ep += ah*e+up(ah-a)*abs(v.real)
                        scale += ah*max(abs(v)+e)
                for k in ks:
                    l = m+n-k
                    if l not in ks:
                        continue
                    _,cn,dn,_ = line.pair(k,l)
                    _,cp,dp,_ = line.pair(k,-1-l)
                    a = theta**2*self.g(n-k)*self.g(n-l)
                    ah = up(theta**2*self.ell_hi[self.outer+n-k]*self.ell_hi[self.outer+n-l])
                    N += a*cn
                    P += a*cp.real
                    en += ah*dn+up(ah-a)*abs(cn)
                    ep += ah*dp+up(ah-a)*abs(cp.real)
                    scale += ah*max(max(abs(cn)+dn),max(abs(cp)+dp))
                f = int(fold(m))
                h = self.space.right*up(line.peak[f]**2)
                for K,E,acc,is_normal in ((N,en,acc_n,True),(P,ep,acc_p,False)):
                    _,s,e,_ = self.space.convolve(grid*up(abs(K)+E,scale),h)
                    row = s[:self.cells].real+e[:self.cells]
                    if is_normal:
                        factor = self.w[f]/self.w[n]
                        acc += factor*row
                        sc_n += factor*max(abs(s[:self.cells])+e[:self.cells])
                    else:
                        acc += row
                        sc_p += max(abs(s[:self.cells])+e[:self.cells])
            # Remaining rainbow rows are summed as scalar areas in the folded basis.
            pos = (abs(n-self.m)>2*self.j)
            neg = (n+self.m+1>2*self.j)
            rest = (pos*self.raise_to[0][n]+neg*self.raise_to[1][n])@line.p2
            rest += self.transfer_wing(line.width**-2)/(2*self.outer-2*self.cut+3)
            normal[n] = min(normal[n],up(max(acc_n)+theta**2*tail[1][n]/W[n],sc_n))
            pair[n] = min(pair[n],up(max(acc_p)+theta*rest+theta**2*tail[2][n],sc_p+rest))
        return np.r_[normal[:self.n],max(normal[self.n:])],max(pair)

    def inverse(self,line,theta,mutant=None):
        n = self.n
        dim = 2*n
        W = self.w[:n]
        g,h = (theta*z[:n,:n]*W[None,:]/W[:,None] for z in self.lower)
        a,b = (theta*z[:n,:n]*W[None,:]/W[:,None] for z in self.raise_to)
        J = np.zeros((dim,dim))
        UJ = J.copy()
        J[:n,:n],J[n:,n:] = g+h,g-h
        UJ[:n,:n] = UJ[n:,n:] = up(up(a-g)+up(b-h),a+b)
        s = line.sqavg[:n]
        X = np.diag(np.r_[s.real,s.real])
        X[n:,:n],X[:n,n:] = np.diag(s.imag),np.diag(-s.imag)
        err = line.sqavg_err[:n]
        UX = np.diag(np.r_[err,err])
        UX[n:,:n] = UX[:n,n:] = np.diag(err)
        K = X@J
        inv = np.linalg.inv(np.eye(dim)-K)
        L = np.rint((J@inv)*2**34)*2**-34
        if mutant=='inverse':
            L[n,n] -= 32
        # Entrywise interval product bounds, before absolute row reduction.
        # Center matrices are signed and only L is declared an exact dyadic matrix.
        EK = up(abs(X)@UJ+UX@(abs(J)+UJ),abs(X)@abs(J))
        prod = L@K
        ep = up(abs(L)@EK,abs(L)@(abs(K)+EK))
        remainder = up(abs(J-L+prod)+UJ+ep,abs(J)+abs(L)+abs(L)@(abs(K)+EK))
        # Vector of complex row norms: commuting real blocks are not assumed.
        def complex_bounds(M):
            a,b,c,d = abs(M[:n,:n]),abs(M[:n,n:]),abs(M[n:,:n]),abs(M[n:,n:])
            bound1 = up(np.maximum(a+c,b+d))
            boundinf = up(np.maximum(a+b,c+d))
            return np.minimum(up(np.sqrt(a*a+b*b+c*c+d*d)),up(np.sqrt(bound1*boundinf)))
        bL = complex_bounds(L)
        def apply(vec):
            return up(np.r_[vec[:-1]+bL@(line.p2[:n]*vec[:-1]),vec[-1]])
        return apply,up(np.r_[complex_bounds(remainder)@line.p2[:n],0.])


def verify(lam,eta,x,primary,witness,proto=Independent):
    n = primary.n
    b = proto(lam,eta,x,n,2*primary.cells,primary.j+3)
    require(len(witness)>1 and not np.any(witness[0]),'independent zero coupling witness')
    old = np.repeat(witness[0],2,axis=1)
    record = b.anchor(old,0.)
    maxq = 0.
    slabs = len(witness)-1
    for i,c in enumerate(witness[1:]):
        require(c.shape==(n,primary.cells) and np.all(np.isfinite(c)),'independent radial center')
        next_center = np.repeat(c,2,axis=1)
        delta = up(np.max(abs(next_center-old)/b.w[:n,None]))
        r = b.anchor(next_center,(i+1)/slabs)
        outcome,reason = b.finish(record,r,i/slabs,(i+1)/slabs,delta)
        if outcome is None:
            return None,f'second slab {i+1}/{slabs}: {reason}'
        maxq = max(maxq,outcome[1])
        old,record = next_center,r
    pair = b.pairing(record,outcome[2])
    if not pair < 1:
        return None,f'second pairing={pair:.6g}'
    return (pair,maxq),''


def mutants():
    """Independent analytic witnesses for the actual second implementation."""
    from fractions import Fraction as Q

    m = HalfLine(1.,32)
    h = np.ones(32)
    expected = Q(3)-Q(9,4)*Q(1,64)
    def contains():
        _,v,e,_ = m.convolve(h,h)
        return abs(Q(float(v[0].real))-expected) <= Q(float(e[0]))
    require(contains(),'healthy odd-convolution fixture')
    m.reverse_cross = True
    require(not contains(),'accepted wrong half-line cross sign')

    r = Independent(Q(2),Q(1),Q(1,16),4,32,7)
    # One positive discarded term remains with all represented profiles zero.
    s = r.outer+1
    charge = Q(r.a0)/((1+Q(r.b1)*s*s)*(2*s-1))
    p = np.zeros(r.count)
    require(Q(float(r.sums(p,1.)[1][0])) > charge,'healthy omitted-transfer fixture')
    r.omit_wing = True
    require(Q(float(r.sums(p,1.)[1][0])) < charge,'accepted omitted infinite wing')
    r.omit_wing = False

    r = Independent(Q(1,2),Q(1),Q(1),4,32,7)
    c = np.zeros((r.n,r.cells),complex)
    left,right = r.anchor(c,0.),r.anchor(c,1/64)
    good,why = r.finish(left,right,0.,1/64,0.)
    require(good is not None,f'healthy inverse defect fixture: {why}')
    bad,why = r.finish(left,right,0.,1/64,0.,'inverse')
    require(bad is None and why=='inverse defect >= 1','accepted corrupted inverse')
    print('RADIAL FIXTURES PASS: 3/3 mutants rejected')


if __name__ == '__main__':
    mutants()
