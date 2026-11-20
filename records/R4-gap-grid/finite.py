"""Sign-aware dressed-rainbow tube and full-tail crossed Schur bounds."""
import argparse
import os
from fractions import Fraction as Q

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np

from finite_check import Check
from replay import C_BOUND, ETA, GRID, LAM, VOLUME, box, choose, decimal, iv, require

# Absolute + relative pads for nonnegative array reductions and elementary expressions.
PAD = 2.0**-35
V = float(VOLUME)
# Newly closed rows, in the fixed LAM order. The first peak-box rows are not repeated.
PUBLISHED = (
    ((Q(1,4),Q(1,2)),(1,),(2,),(4,),(4,8)),
    ((Q(1,8),Q(1,4),Q(1,2)),(1,),(2,),(2,4),(4,)),
)


def up(a):
    return np.nextafter(a + PAD*(1+abs(a)), np.inf)


def down(a):
    return np.nextafter(a - PAD*(1+abs(a)), -np.inf)


def endpoints(a):
    return float(np.nextafter(float(a.a),-np.inf)), float(np.nextafter(float(a.b),np.inf))


class Tube:
    def __init__(self, lam, eta, x):
        self.x = float(x)
        self.xlo,self.xhi = endpoints(box(x))
        self.alo,self.ahi = endpoints(2*box(lam/x)/iv.pi)
        self.blo,self.bhi = endpoints((2*iv.pi*C_BOUND)**2*box(lam*eta**2))
        self.n = 96
        self.j = 2*self.n
        self.max_m = self.n+self.j
        self.m = np.arange(self.max_m)
        self.w = 2*self.m+1
        self.rows = np.arange(self.n)
        self.s = np.arange(-self.j,self.j+1)
        mm = self.rows[:,None]-self.s
        self.idx = np.where(mm>=0,mm,-mm-1)
        self.g = up(self.ahi/(1+self.blo*self.s**2))
        self.R = up(sum(self.g)+2*self.ahi/(self.blo*self.j))
        # Evaluate the exact nonnegative reflection difference without subtracting kernels.
        f = (self.rows[:,None]-self.m)**2
        b = (self.rows[:,None]+self.m+1)**2
        self.fold_hi = up(self.ahi*self.bhi*(b-f)/((1+self.blo*f)*(1+self.blo*b)))
        self.fold_lo = np.maximum(0,down(self.alo*self.blo*(b-f)/((1+self.bhi*f)*(1+self.bhi*b))))
        self.weight = np.where(self.w<17,1.,4.)
        self.row_weight = np.r_[self.weight[:self.n],4.]
        self.rainbow_tail = up(2*self.ahi/(self.blo*(self.max_m-self.n)*(2*self.max_m+1)))
        self.low = 8
        self.local = abs(self.s)<=7
        tr = self.s[self.local]
        n = self.rows[:self.low,None,None]
        fold = lambda m: np.where(m>=0,m,-m-1)
        self.panel = (fold(n-tr[None,:,None]),fold(n-tr[None,:,None]-tr[None,None,:]),
                      fold(n-tr[None,None,:]))
        k = self.g[self.local]
        self.panel_mass = k[:,None]*k[None,:]

    def mean(self, v):
        return self.x*np.arctan(1/(self.x*v))


    def centers(self, steps):
        d = np.zeros(self.max_m)
        yield d.copy()
        for i in range(1,steps+1):
            t = i/steps
            for _ in range(30):
                d[:self.n] = t*(self.fold_hi@self.mean(self.w+d))
            yield d.copy()

    def sums(self, profiles, t, radius):
        # Tail of omitted transfer in each head row: n < N, |s| > J.
        h = 4*radius
        tail_peak = 1/(2*(self.j+1)-(2*self.n-1)-h)
        tail_g = 2*self.ahi/(self.blo*self.j)
        # Beyond n>=N, split on low electron |2m+1|<N, not on the transfer.
        high_peak = 1/(self.n-h)
        low_g = self.ahi/(1+self.blo*(self.n//2)**2)
        total = []
        for p,power,scale in profiles:
            mid = p[self.idx]@self.g
            head = up(mid + tail_g*scale*tail_peak**power)
            low = 2*np.sum(p[self.w<self.n])
            tail = up(self.R*scale*high_peak**power + low_g*low)
            total.append(up(t*np.r_[head,tail]))
        return total

    def line(self,v,radius):
        v = down(v)
        floor = down(v-radius*self.weight)
        if np.min(floor) <= 0:
            return None
        peak = up(1/floor)
        # (sqrt(v²+y²)-rho)² >= (v-rho)² + (1-rho/v)y².
        u,beta = floor,floor/v
        root = np.sqrt(beta)
        arg = root/(self.xlo*u)
        a = up(np.arctan(arg))
        dos2 = up(np.minimum(1/u**2,self.xhi*a/(u*root)))
        primitive = up(np.arcsinh(arg))
        l1 = up(np.minimum(np.sqrt(dos2),self.xhi*primitive/root))
        l2 = up(np.sqrt(dos2))
        # Integral and L2 norm of the squared Lorentz majorant.
        l4sq = up(np.sqrt(np.minimum(u**-4, 1/(2*u*u*(u*u+beta/self.xhi**2))
                                     +self.xhi*a/(2*u**3*root))))
        return l1,l2,peak,dos2,l4sq

    def resolve(self, t, radius, line, derivative=None, cap=1):
        """Resolve both low transfers explicitly and bound only the remaining wings."""
        def split(p,power,scale):
            cols = self.idx[:self.low]
            h = t*(p[cols[:,self.local]] @ self.g[self.local])
            bound = 2*(self.j+1)-(2*self.low-1)-4*radius
            tail = 2*self.ahi/(self.blo*self.j)*scale/bound**power
            r = t*(p[cols[:,~self.local]] @ self.g[~self.local]+tail)
            return up(h),up(r)

        def join(a,b):
            ha,ra = a
            hb,rb = b
            return ra*(hb+rb)+ha*rb

        def direct(a,b,c):
            A,B,C = (np.asarray(p)[:,i] for p,i in zip((a,b,c),self.panel))
            bound = np.minimum.reduce((A[0]*B[1]*C[1],A[1]*B[0]*C[1],A[1]*B[1]*C[0],
                                       V*A[2]*B[2]*C[2]))
            return up(t*t*np.sum(self.panel_mass*bound,axis=(1,2)))

        s1,s2,sp = (split(p,1,1) for p in line)
        u1,u2,upk = (max(p) for p in line)
        side1,side2,sidep = join(s1,s2),join(s2,s2),join(sp,sp)
        if derivative is None:
            tail = np.minimum.reduce((u1*side2,u2*side1,V*upk*sidep))
            return up(direct(line,line,line)+tail)
        a1,a2,ap = (split(p,2,cap) for p in derivative)
        left = np.minimum.reduce((u1*join(a2,s2),u2*join(a1,s2),u2*join(a2,s1),V*upk*join(ap,sp)))
        middle = np.minimum.reduce((max(derivative[0])*side2,max(derivative[1])*side1,
                                    V*max(derivative[2])*sidep))
        return up(2*(direct(derivative,line,line)+left)+direct(line,derivative,line)+middle)

    def cross(self, t, v, radius, normal, parts=False, data=None, resolve=True):
        l1,l2,peak,dos2,l4sq = self.line(v,radius) if data is None else data
        s1,s2,sp = self.sums(((l1,1,1),(l2,1,1),(peak,1,1)),t,radius)
        u1,u2,upeak = (max(p) for p in (l1,l2,peak))
        if normal is None:
            q = up(np.minimum.reduce((u1*s2**2,u2*s1*s2,V*upeak*sp**2)))
            if resolve:
                q[:self.low] = np.minimum(q[:self.low],self.resolve(t,radius,(l1,l2,peak)))
            return q
        w = self.weight if normal else 1.
        cap = 4. if normal else 1.
        d1,d2,dp = dos2*w,l4sq*w,peak**2*w
        a1,a2,ap = self.sums(((d1,2,cap),(d2,2,cap),(dp,2,cap)),t,radius)
        ends = np.minimum.reduce((u1*a2*s2,u2*a1*s2,u2*a2*s1,V*upeak*sp*ap))
        middle = np.minimum.reduce((max(d1)*s2**2,max(d2)*s1*s2,V*max(dp)*sp**2))
        cross = up(2*ends+middle)
        if resolve:
            panel = self.resolve(t,radius,(l1,l2,peak),(d1,d2,dp),cap)
            cross[:self.low] = np.minimum(cross[:self.low],panel)
        if normal and parts:
            return up(cross/self.row_weight),up(a1/self.row_weight)
        return up((a1+cross)/self.row_weight if normal else a1+cross)

    def inverse(self,t,v,central):
        # |∫ Z_c²| = x²/(1+(xv)²), not ∫|Z_c|².
        j = up(np.minimum(1/v**2,self.xhi**2/(1+(self.xhi*v)**2)))
        a,d = self.sums(((central[3]*self.weight,2,4),(j*self.weight,2,4)),t,0)
        ell = up(max(a/self.row_weight))
        ell2 = up(ell*max(d/self.row_weight))
        gap = down(1-ell2)
        if gap > 0:
            return up((1+ell)/gap)
        return np.inf

    def normal(self,t,v,radius,inv,central):
        ball = self.line(v,radius)
        crossed,rainbow = self.cross(t,v,radius,True,True,ball)
        plain = up(max(crossed+rainbow))
        if not np.isfinite(inv):
            return plain,np.inf
        importance = up(radius*self.weight**2*(ball[2]+central[2])*ball[1]*central[1])
        diff, = self.sums(((importance,3,2*radius*4**2),),t,radius)
        return plain,up(inv*max(crossed+diff/self.row_weight))

    def residual(self, lo, hi, dl, dh):
        v0 = self.w+np.minimum(dl,dh)
        endpoint = []
        for t,d in ((lo,dl),(hi,dh)):
            v = self.w+d
            means_lo = np.maximum(0,down(self.xlo*np.arctan(1/(self.xhi*v))))
            means_hi = up(self.xhi*np.arctan(1/(self.xlo*v)))
            low_t,high_t = max(0.,np.nextafter(t,-np.inf)),np.nextafter(t,np.inf)
            rl = down(down(low_t*(self.fold_lo@means_lo))-d[:self.n])
            rh = up(up(high_t*(self.fold_hi@means_hi+self.rainbow_tail))-d[:self.n])
            endpoint.append((rl,rh))
        lo,hi = max(0.,np.nextafter(lo,-np.inf)),np.nextafter(hi,np.inf)
        # Affine interpolation cancels the large first-order rainbow center change.
        dd = up(abs(dh-dl))
        bound = 2*(hi-lo)*np.minimum(self.xhi**2,v0**-2)*dd
        bound += 2*hi*np.minimum(self.xhi**3,v0**-3)*dd**2
        modulus = up((self.fold_hi@bound)/8)
        low = down(np.minimum(endpoint[0][0],endpoint[1][0])-modulus)
        high = up(np.maximum(endpoint[0][1],endpoint[1][1])+modulus)
        signed = np.maximum(abs(low),abs(high))
        # For the zero tail center, the absolute rainbow estimate suffices.
        pt = self.line(v0,0)[0]
        tail = self.sums(((pt,1,1),),hi,0)[0][-1]
        r = up(np.r_[signed,tail]+self.cross(hi,v0,0,None))
        return v0,up(max(r/np.r_[self.weight[:self.n],4.]))

    def verify(self, steps):
        centers = self.centers(steps)
        dl = next(centers)
        witness = []
        max_q = 0.
        for i,dh in enumerate(centers):
            lo,hi = i/steps,(i+1)/steps
            v, res = self.residual(lo,hi,dl,dh)
            hi = np.nextafter(hi,np.inf)
            central = self.line(v,0)
            inv = self.inverse(hi,v,central)
            for factor in (1.1,1.25,1.45,1.7,2,2.5,3,3.6,4.3,5.2,6.3,7.7,9.3,11.5):
                radius = up(factor*res)
                plain = pre = np.inf
                applied = False
                if min(down(down(v)-radius*self.weight)) <= 1/8:
                    break
                plain,pre = self.normal(hi,v,radius,inv,central)
                if res < down(radius*(1-plain)):
                    max_q = max(max_q,plain)
                    break
                if np.isfinite(inv) and up(inv*res) < down(radius*(1-pre)):
                    max_q = max(max_q,pre)
                    applied = True
                    break
                if plain >= 1 and pre >= 1:
                    break
            else:
                return None, f'slab {i+1}/{steps}: residual={res:.6g}'
            if min(down(down(v)-radius*self.weight)) <= 1/8 or not (
                    res < down(radius*(1-plain))
                    or np.isfinite(inv) and up(inv*res) < down(radius*(1-pre))):
                return None, f'slab {i+1}/{steps}: residual={res:.6g}, radius={radius:.6g}'
            witness.append((dl[:self.n].copy(),dh[:self.n].copy(),radius,applied))
            dl = dh
        value = max(self.cross(1.,self.w+dh,radius,False))
        if not 0 <= value < 1:
            return None, f'normal q={max_q:.6g}, pairing={value:.6g}'
        return (value,max_q,witness), ''


def fixtures():
    T = Tube(Q(1,4),Q(3,2),Q(1,2))
    a = 2*box(Q(1,2))/iv.pi
    b = (2*iv.pi*C_BOUND)**2*box(Q(9,16))
    signed = a*b/(1+b)
    positive = a+a/(1+b)
    require(box(T.fold_lo[0,0]) < signed < box(T.fold_hi[0,0]),'signed reflection fixture')
    require(positive > box(T.fold_hi[0,0]),'wrong reflection sign survived')

    w = 1/T.w
    y = T.sums(((w,1,1),),1.,0)[0][T.n-1]
    n = T.n-1
    unfolded = sum(a/((1+b*j*j)*abs(2*(n-j)+1)) for j in range(-2*T.j,2*T.j+1))
    require(unfolded < box(y),'transfer-tail fixture')
    missing = w[T.idx[n]] @ T.g
    require(unfolded > box(missing),'missing tail survived')

    # Independent derivative of the spatially constant static cube at unit amplitude.
    t = Q(1,2**12)
    derivative = ((1+t)**3-(1-t)**3)/(2*t)-t*t
    g = 0.5
    T.sums = lambda profiles, *_: [np.full(T.n+1,g) for p in profiles]
    T.weight = np.ones(T.max_m)
    unit = np.ones(T.max_m)
    crossed,_ = T.cross(1.,T.w,0,True,True,(unit,)*5,resolve=False)
    goal = derivative*VOLUME*Q(g)**2
    require(box(goal) < box(crossed[0]) < box(goal+Q(1,2**22)),'all three insertions fixture')
    require(box(2*crossed[0]/3) < box(goal),'missing insertion survived')
    print('FINITE MUTANTS rejected=3 (reflection, cutoff tail, middle insertion)',flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--steps',type=int,default=48)
    p.add_argument('--lam',type=Q,nargs='+',default=LAM)
    p.add_argument('--eta',type=Q,nargs='+',default=ETA[:2])
    p.add_argument('--xs',type=Q,nargs='+',default=GRID)
    p.add_argument('--published',action='store_true')
    args = p.parse_args()
    require(1 <= args.steps <= 256,'steps out of arithmetic bounds')
    require(set(args.lam)<=set(LAM) and set(args.eta)<=set(ETA[:2]) and set(args.xs)<=set(GRID)
            and all(len(set(a))==len(a) for a in (args.lam,args.eta,args.xs)), 'finite grid inputs')
    fixtures()
    counts = dict.fromkeys(args.eta,0)
    total = 0
    closed = set()
    worst = 0.
    print('lambda,eta,x,pairing_upper,normal_upper,preconditioned_slabs',flush=True)
    for lam in args.lam:
        for eta in args.eta:
            for x in args.xs:
                if choose(lam,eta,x) is not None:
                    continue
                total += 1
                result, reason = Tube(lam,eta,x).verify(args.steps)
                if result is not None:
                    ev,q,witness = result
                    second, reason = Check(lam,eta,x,len(witness[0][0])).validate(witness)
                    result = None if second is None else (max(ev,second[0]),max(q,second[1]),second[2])
                if result is None:
                    print(f'{lam},{eta},{x},UNRESOLVED {reason}',flush=True)
                    continue
                counts[eta] += 1
                closed.add((lam,eta,x))
                worst = max(worst,result[0])
                print(f'{lam},{eta},{x},{decimal(result[0],8)},{decimal(result[1],8)},{result[2]}',flush=True)
    print(f'NEW closed={sum(counts.values())}/{total} slices={list(counts.values())}',flush=True)
    if args.published:
        expected = {(lam,eta,x) for eta,rows in zip(ETA[:2],PUBLISHED)
                    for lam,xs in zip(LAM,rows) for x in xs}
        require(args.steps == 48 and total == 59 and closed == expected
                and counts == {Q(1):7,Q(3,2):8} and box(worst) < box(Q('0.50417617')),
                'published finite-grid claim not reproduced')
        print(f'FINITE PASS worst pairing upper<={decimal(worst,8)}; unresolved={total-len(closed)}',
              flush=True)


if __name__ == '__main__':
    main()
