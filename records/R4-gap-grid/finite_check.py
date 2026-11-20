"""Independent folded-frequency, exact integer-dot verifier of dressed tubes."""
from fractions import Fraction as Q

import numpy as np

from replay import C_BOUND, box, iv, require

UNIT = 2.0**-33
DEN = 2**22


def enclosure(value, sign=1):
    return np.nextafter(value + sign*UNIT*(1+abs(value)), sign*np.inf)


def scalar(interval, sign):
    return float(np.nextafter(float(interval.a if sign<0 else interval.b),sign*np.inf))


def atan(z):
    # Two half-angle reductions shorten the domain to tan(pi/16).
    reflected = z > 1
    t = np.minimum(z,1/z)
    for _ in range(2):
        t = t/(1+np.sqrt(1+t*t))
    power = t.copy()
    result = power.copy()
    for j in range(3,28,2):
        power = -power*t*t
        result += power/j
    return np.where(reflected,np.pi/2-4*result,4*result)


def quantize(a, upper):
    a = np.maximum(a,0)
    return (np.ceil(a*DEN) if upper else np.floor(a*DEN)).astype(np.int64)


def dot(k, p, upper=True):
    v = quantize(p,upper)
    require(k.shape[1]*int(np.max(k))*int(np.max(v)) < 2**63,'exact integer dot overflow')
    value = (k@v).astype(float)/(DEN*DEN)
    return np.nextafter(value, np.inf if upper else -np.inf)

def inverse_moment(v, x, weight):
    v_lower = enclosure(v,-1)
    return enclosure(weight/(v_lower*v_lower+1/x**2))



class Check:
    def __init__(self, lam, eta, x, n):
        self.b = (2*iv.pi*C_BOUND)**2*box(lam*eta*eta)
        self.a = 2*box(lam/x)/iv.pi
        self.x = float(x) # Certified grid values are exactly binary.
        self.n = n
        self.m = 3*n+9
        self.w = 2*np.arange(self.m)+1
        self.weight = np.where(self.w<17,1.,4.)
        self.wn = np.r_[self.weight[:n],4.]
        h = np.arange(n)[:,None]
        e = np.arange(self.m)
        alo,ahi = scalar(self.a,-1),scalar(self.a,1)
        blo,bhi = scalar(self.b,-1),scalar(self.b,1)
        kp,km = (h-e)**2,(h+e+1)**2
        up = lambda a,k: enclosure(a/(1+blo*k))
        down = lambda a,k: enclosure(a/(1+bhi*k),-1)
        self.plus = quantize(enclosure(up(ahi,kp)+up(ahi,km)),True)
        self.fold_hi = quantize(enclosure(up(ahi,kp)-down(alo,km)),True)
        self.fold_lo = quantize(enclosure(down(alo,kp)-up(ahi,km),-1),False)
        self.loss = scalar(2*self.a/(self.b*(self.m-self.n)),1)
        self.mass = scalar(sum(self.a/(1+self.b*j*j) for j in range(-self.m,self.m+1))
                           +2*self.a/(self.b*self.m),1)
        # Different high-row split from the primary: low m with w_m < (5/4)N.
        self.cut = n+n//4
        self.dist = n-(self.cut-2)//2
        self.leak = scalar(self.a/(1+self.b*self.dist*self.dist),1)
        # Group the direct panel by external row and flattened transfer pair.
        s = np.arange(-9,10)
        left,right = (a.ravel() for a in np.meshgrid(s,s,indexing='ij'))
        sign_fold = lambda k: np.maximum(k,-k-1)
        q = np.arange(8)[:,None]
        self.nodes = tuple(sign_fold(q-u) for u in (left,left+right,right))
        self.panel_mass = up(ahi,left*left)*up(ahi,right*right)
        self.local = sign_fold(q-s)
        self.local_lo = np.maximum(0,down(alo,s*s))
        self.local_hi = up(ahi,s*s)

    def sum(self,p,power,weight,t,radius):
        worst = 4*radius
        head = dot(self.plus,p)+self.loss*weight/(2*self.m+1-worst)**power
        hi = self.mass*weight/(self.cut-worst)**power
        lo = 2*self.leak*np.sum(p[self.w<self.cut])
        return enclosure(t*np.r_[head,hi+lo])

    def profile(self,v,r):
        v = enclosure(v,-1)
        u = enclosure(v-r*self.weight,-1)
        require(min(u)>1/8,'insufficient independent widths')
        b = np.sqrt(u/v)
        y = b/(self.x*u)
        integr = enclosure(atan(y))
        l1 = enclosure(self.x*np.log(y+np.sqrt(1+y*y))/b)
        f2 = enclosure(np.minimum(1/(u*u),self.x*integr/(u*b)))
        l2 = enclosure(np.sqrt(f2))
        l1 = np.minimum(l1,l2)
        f4 = enclosure(np.sqrt(np.minimum(u**-4,self.x*integr/(2*u**3*b)
                                          +1/(2*u*u*(u*u+(b/self.x)**2)))))
        return np.stack((l1,l2,enclosure(1/u),f2,f4))

    def crossed_panel(self,t,base,d=None):
        def spatial(profiles):
            a,b,c = (p[:,node] for p,node in zip(profiles,self.nodes))
            u = a[0]*b[1]*c[1]
            u = np.minimum(u,a[1]*b[0]*c[1])
            u = np.minimum(u,a[1]*b[1]*c[0])
            u = np.minimum(u,(15/32)*a[2]*b[2]*c[2])
            return enclosure(t*t*np.sum(u*self.panel_mass,axis=1))

        def sides(p,sums):
            h = t*np.einsum('ijk,k->ij',p[:,self.local],self.local_hi)
            lo = t*np.einsum('ijk,k->ij',p[:,self.local],self.local_lo)
            return enclosure(h),enclosure(np.maximum(0,sums[:,:8]-enclosure(lo,-1)))

        def w(h,i,j):
            a,b = h
            return a[i]*b[j]+b[i]*(a[j]+b[j])

        left = sides(*base)
        supremum = np.max(base[0],axis=1)
        outer = np.stack((w(left,1,1),w(left,0,1),(15/32)*w(left,2,2)))
        if d is None:
            wing = np.min(supremum[:,None]*outer,axis=0)
            return enclosure(spatial((base[0],)*3)+wing)
        h,r = sides(*d)
        s,y = left
        mix = r*(s+y)+h*y
        ends = np.minimum.reduce((supremum[0]*mix[1],supremum[1]*(r[0]*(s[1]+y[1])+h[0]*y[1]),
                                  supremum[1]*(r[1]*(s[0]+y[0])+h[1]*y[0]),
                                  (15/32)*supremum[2]*mix[2]))
        middle = np.min(np.max(d[0],axis=1)[:,None]*outer,axis=0)
        return enclosure(2*(spatial((d[0],base[0],base[0]))+ends)
                         +spatial((base[0],d[0],base[0]))+middle)

    def product(self,t,p,r,kind):
        _,_,peak,dsq,dfour = p
        s1,s2,sp = (self.sum(col,1,1,t,r) for col in p[:3])
        l1,l2,linf = np.max(p[:3],axis=1)
        v = 15/32
        if kind=='image':
            q = enclosure(np.minimum(np.minimum(l1*s2*s2,l2*s1*s2),v*linf*sp*sp))
            q[:8] = np.minimum(q[:8],self.crossed_panel(t,(p[:3],np.stack((s1,s2,sp)))))
            return q
        w, cap = (self.weight,4) if kind=='normal' else (1,1)
        d,e,f = dsq*w,dfour*w,peak*peak*w
        a,b,c = (self.sum(col,2,cap,t,r) for col in (d,e,f))
        left = np.minimum(np.minimum(a*l2*s2,b*l1*s2),np.minimum(b*l2*s1,v*c*linf*sp))
        center = np.minimum(np.minimum(max(d)*s2*s2,max(e)*s1*s2),v*max(f)*sp*sp)
        q = enclosure(2*left+center)
        q[:8] = np.minimum(q[:8],self.crossed_panel(t,(p[:3],np.stack((s1,s2,sp))),
                                                  (np.stack((d,e,f)),np.stack((a,b,c)))))
        return a,q

    def rainbow_error(self,t0,t1,d0,d1,v):
        def defect(t,d):
            mean = self.x*atan(1/(self.x*(self.w+d)))
            a,b = max(0,np.nextafter(t,-np.inf)),np.nextafter(t,np.inf)
            low = enclosure(enclosure(a*dot(self.fold_lo,enclosure(mean,-1),False),-1)-d[:self.n],-1)
            high = enclosure(enclosure(b*(dot(self.fold_hi,enclosure(mean))+self.loss/(2*self.m+1)))-d[:self.n])
            return low,high
        lo,hi = max(0,np.nextafter(t0,-np.inf)),np.nextafter(t1,np.inf)
        # Independently bound the affine-interpolation error with exact DOS derivatives.
        delta = abs(d1-d0)
        derivative = 1/(v*v+1/self.x**2)
        curvature = 2*derivative*derivative*np.minimum(v+delta,1/np.sqrt(derivative))
        modulus = enclosure(delta*(2*(hi-lo)*derivative+hi*delta*curvature))
        error = enclosure(dot(self.fold_hi,modulus)/8)
        e0,e1 = defect(t0,d0),defect(t1,d1)
        lower = enclosure(np.minimum(e0[0],e1[0])-error,-1)
        upper = enclosure(np.maximum(e0[1],e1[1])+error)
        return np.maximum(abs(lower),abs(upper))

    def validate(self, witness):
        max_q = 0.
        used = 0
        slabs = len(witness)
        require(slabs>0 and not np.any(witness[0][0]),'initial normal center')
        prior = witness[0][0]
        for i,(d0,d1,r,pre) in enumerate(witness):
            require(np.array_equal(prior,d0) and d1.shape==(self.n,)
                    and np.all(np.isfinite(d1)) and min(d1)>=0 and np.isfinite(r) and r>0,
                    'inhomogeneous endpoints or invalid radius')
            prior = d1
            t0,t1 = i/slabs,(i+1)/slabs
            hi = np.nextafter(t1,np.inf)
            d0,d1 = (np.pad(d,(0,self.m-self.n)) for d in (d0,d1))
            v = self.w+np.minimum(d0,d1)
            central = self.profile(v,0)
            p = self.profile(v,r)
            err = self.rainbow_error(t0,t1,d0,d1,v)
            tail = self.sum(central[0],1,1,hi,0)[-1]
            cubic = self.product(hi,central,0,'image')
            resid = max(enclosure((np.r_[err,tail]+cubic)/self.wn))
            a,cross = self.product(hi,p,r,'normal')
            if not pre:
                q = max(enclosure((a+cross)/self.wn))
            else:
                norm = max(enclosure(self.sum(central[3]*self.weight,2,4,hi,0)/self.wn))
                j = inverse_moment(v,self.x,self.weight)
                jj = max(enclosure(self.sum(j,2,4,hi,0)/self.wn))
                gap = enclosure(1-norm*jj,-1)
                if gap<=0:
                    return None, f'independent inverse at slab {i+1}'
                inverse = enclosure((1+norm)/gap)
                err = r*self.weight**2*(p[2]+central[2])*p[1]*central[1]
                excess = self.sum(enclosure(err),3,32*r,hi,r)
                resid = enclosure(inverse*resid)
                q = max(enclosure(inverse*(cross+excess)/self.wn))
                used += 1
            if not 0 <= q < 1 or not resid < enclosure(r*(1-q),-1):
                return None, f'independent slab {i+1}/{slabs}: q={q:.6g}, residual={resid:.6g}'
            max_q = max(max_q,q)
        a,cross = self.product(1.,self.profile(self.w+d1,r),r,'pair')
        ev = max(enclosure(a+cross))
        if not 0 <= ev < 1:
            return None, f'independent pairing q={ev:.6g}'
        return (ev,max_q,used), ''


if __name__ == '__main__':
    # Exact replay of the reviewer's inward-quantization witness.
    v = 1+1.9998346544058268
    x = .125
    exact = 1/(Q(v)**2+1/Q(x)**2)
    require(Q(int(quantize(1/(v*v+1/x**2),True)),DEN)<exact,'rounding witness lost')
    require(Q(int(quantize(inverse_moment(v,x,1),True)),DEN)>=exact,'inverse moment quantization')
    print('FINITE CHECK inverse-moment rational rounding fixture PASS')
