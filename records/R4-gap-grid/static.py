"""Static-transfer radial solver and continuum enclosure on dyadic grid points."""
import argparse
import os
from fractions import Fraction as Q

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from scipy.signal import fftconvolve
from replay import VOLUME, guards, require, box, iv, decimal, GRID, LAM


# Arithmetic pads are specified and bounded in RIGOR.md, not convergence tolerances.
PAD = 2.0**-26
SIGN_PAD = 2.0**-42


def shifted_correlation(a,b):
    """Independent folding of the positive and negative half-lines."""
    out = np.zeros(len(a)+len(b),dtype=np.result_type(a,b))
    out[1:] = np.convolve(a,b)
    out[:len(a)] -= np.correlate(a,b.conj(),mode='full')[len(b)-1:]
    out[:len(b)] -= np.correlate(b,a.conj(),mode='full')[len(a)-1:]
    return out


INDEPENDENT = False

def up(v):
    return np.nextafter(v + PAD * (1 + abs(v)), np.inf)



def height(values):
    return np.r_[-values[::-1], values]


def convol(a, b, h, rigorous):
    """For constant odd H_a,H_b, return C at left knots, C/H at midpoints, and range radius."""
    i = len(a) + len(b) - 1
    if rigorous and INDEPENDENT:
        conv = shifted_correlation(a,b)
    else:
        op = np.convolve if rigorous else fftconvolve
        conv = op(height(a), height(b))[i:]
    d = -(1.5*h) * np.r_[conv,0.]
    primitive = np.r_[0., np.cumsum(h/2 * (d[1:] + d[:-1]))]
    index = np.arange(len(d)-1)
    hl = primitive[:-1]
    left = np.r_[d[0], hl[1:] / (h*index[1:])]
    mid_a = hl + (h/8) * (3*d[:-1] + d[1:])
    mid = mid_a / (h*(index+.5))
    local = np.maximum(abs(d[:-1]-left), abs(d[1:]-left)) / (index+1)
    dd = np.maximum(abs(d[:-1]), abs(d[1:]))
    if not rigorous:
        return mid, mid_a
    scale = (1 + np.max(abs(a))) * (1 + np.max(abs(b)))
    limit = h*(len(a)+len(b))
    err = PAD*scale
    return left, mid, mid_a, up(local+2*err), up(dd+err), err, (limit+1)*err


class Static:
    def __init__(self, x, cells):
        self.cells = cells
        self.x = float(x)
        self.h = 1/cells
        self.r = (np.arange(cells)+.5)*self.h
        self.rr = self.r+self.h/2
        self.d = (np.arange(2*cells)+1)*self.h
        self.slope = 6*self.rr**2
        self.de = self.slope*self.h/2
        self.bare = 1+1j*self.x-2*self.r**3
        self.measure = 3*self.h

    def update(self, sigma, g):
        z = self.r/(self.bare-sigma)
        _, c = convol(z, z, self.h, False)
        xc, _ = convol(c, z, self.h, False)
        return g*self.measure*np.dot(self.r,z) + g*g*xc[:self.cells]

    def solve(self, sigma, g):
        for _ in range(100):
            step = self.update(sigma, g)-sigma
            if np.max(abs(step)) < self.x*2e-10:
                return sigma
            sigma = sigma + .72*step
        raise ArithmeticError('static proposal failed')

    def moments(self, sigma, radius):
        denominator = self.bare - sigma
        smallest = up(np.maximum(0, abs(denominator.real)-self.de))
        smallest -= 2*PAD*(1+smallest)
        min_norm = np.hypot(np.maximum(0,smallest), denominator.imag)
        min_norm -= 2*PAD*(1+min_norm)
        if np.min(denominator.imag) <= self.x/2+radius+PAD or np.min(min_norm) <= radius:
            return None
        peak = up(1/(min_norm-radius))
        z = denominator.conj() / (denominator.real**2 + denominator.imag**2)
        z_abs = up(abs(z))
        hz = self.r*z
        slope = up(z_abs + self.rr*z_abs*peak*self.slope)
        constant = up(self.rr*z_abs*peak*radius + 4*PAD*(1+z_abs))
        err = up(self.h/np.sqrt(12)*self.norm(slope)+self.norm(constant))
        err1 = up(self.measure*np.dot(self.rr,self.h/4*slope+constant))
        return hz, peak, err, err1

    def norm(self, hz):
        return up(np.sqrt(self.measure*np.dot(abs(hz),abs(hz))))

    def integral(self, hz):
        mid = self.measure*np.dot(self.r,hz)
        return mid, up(PAD*(1+self.norm(hz)))

    def residual(self, sigma, g_lo, g_hi):
        hz, _, e, e1 = self.moments(sigma, 0)
        a = self.norm(hz)
        mz, me = self.integral(hz)
        _, _, hpair, _, slope, _, herr = convol(hz,hz,self.h,True)
        split = up(self.h/np.sqrt(12)*self.norm(slope)+np.sqrt(6)*herr)
        left, _, _, variation, _, er, _ = convol(hpair,hz,self.h,True)
        extra = up(3*e1*(a+e)**2 + a*split+er)
        outside = left[:self.cells]
        var = variation[:self.cells]
        g = (g_lo+g_hi)/2
        halfspan = np.nextafter(max(g-g_lo,g_hi-g),np.inf)
        triple_scale = np.max(abs(outside)+var)
        signed_scale = 1 + g_hi*abs(mz) + g_hi*g_hi*triple_scale + np.max(abs(sigma))
        signed_pad = SIGN_PAD*signed_scale
        bare_residual = np.max(abs(g*mz+g*g*outside-sigma)+g*g*var)
        # Include whole coupling interval, not just the midpoint.
        return up(bare_residual+signed_pad+g_hi*(me+e1)+g_hi*g_hi*extra
                  +halfspan*(abs(mz)+2*g_hi*triple_scale))

    def schur(self, sigma, radius, g_lo, g_hi, pairing):
        if not np.isfinite(radius):
            return np.inf
        data = self.moments(sigma,radius)
        if data is None:
            return np.inf
        hz, peaks, e, _ = data
        if pairing:
            rr, _, _, rv, _, re, _ = convol(hz.real,hz.real,self.h,True)
            qq, _, _, qv, _, qe, _ = convol(hz.imag,hz.imag,self.h,True)
            corr = 3*rr-qq
            pre_cancel = 3*abs(rr)+abs(qq)
            dev = 3*(rv+re+e*(2*self.norm(hz.real)+e)) + qv+qe+e*(2*self.norm(hz.imag)+e)
        else:
            zz, _, _, zv, _, ze, _ = convol(hz,hz,self.h,True)
            corr = 3*zz
            pre_cancel = 3*abs(zz)
            dev = 3*(zv+ze+e*(2*self.norm(hz)+e))
        signed_pad = SIGN_PAD*(1+g_hi*pre_cancel)
        kernel = up(g_hi*(np.maximum(abs(1+g_lo*corr),abs(1+g_hi*corr))+g_hi*dev+signed_pad))
        left, _, _, var, _, err, _ = convol(up(self.d*kernel),up(self.rr*peaks**2),self.h,True)
        return up(np.max(left[:self.cells]+var[:self.cells])+err)


def coupling(g):
    return float(np.nextafter(float(g.a), -np.inf)), float(np.nextafter(float(g.b), np.inf))


def certify(x, lam, cells, steps, budget, lip):
    model = Static(x,cells)
    g_exact = 2*box(lam*x)/iv.pi
    prev = np.zeros(cells, dtype=complex)
    inner = 0.
    max_q = 0.
    for i in range(steps):
        gl = 0. if i == 0 else coupling(g_exact*i/steps)[0]
        gh = coupling(g_exact*(i+1)/steps)[1]
        sigma = model.solve(prev, (gl+gh)/2)
        res = up(model.residual(sigma, gl, gh)+budget)
        junction = up(np.max(abs(sigma-prev))+inner)
        r0 = max(junction, res)
        for fac in (1.01, 1.20, 1.45, 1.8, 2.2, 2.7, 3.3, 4.1, 5.0, 6.3):
            radius = r0*fac
            q = up(model.schur(sigma,radius,gl,gh,False)+lip)
            if q < 1 and res < radius*(1-q)*(1-PAD):
                break
        else:
            return None, f'slab={i}/{steps}; residual/x={res/float(x):.4g}; overlap/x={junction/float(x):.4g}; q={q:.4g}'
        # At the right endpoint sharpen the inherited distance bound for the next slab.
        end_l,end_h = coupling(g_exact*(i+1)/steps)
        inner = up((model.residual(sigma,end_l,end_h)+budget)/(1-q))
        max_q = max(max_q,q)
        prev = sigma
    pairing = up(model.schur(sigma,inner,end_l,end_h,True)+lip)
    if not np.isfinite(pairing) or pairing >= 1:
        return None, f'pairing row={pairing:.6g}'
    return (pairing,max_q), ''


def angular_fixture(a,b,h,r):
    if not r:
        return 3*h*sum(a[i]*b[i] for i in range(min(len(a),len(b))))
    def primitive(t):
        return sum(val*max(Q(0),min(t,h*(i+1))-h*i) for i,val in enumerate(b))
    def angle(t):
        return primitive(r+t)-primitive(abs(r-t))
    crossing = [r, *[v for j in range(len(b)+1) for v in
                    (j*h-r, r-j*h, r+j*h)]]
    integral = Q(0)
    for i,val in enumerate(a):
        lo,hi = i*h,(i+1)*h
        knots = sorted({lo,hi,*[v for v in crossing if lo<v<hi]})
        integral += val*sum((v-u)*(angle(u)+angle(v))/2 for u,v in zip(knots,knots[1:]))
    return 3*integral/(2*r)


def spatial_fixtures():
    a,b = [2,-1,3,-2], [1,-4,1,2,3]
    h = Q(1,8)
    left,mid,_,var,_,err,_ = convol(np.array(a,dtype=float),np.array(b,dtype=float),float(h),True)
    for i in range(len(left)):
        for q in (Q(0), Q(1,4), Q(1,2), Q(3,4), Q(1)):
            exact = angular_fixture(a,b,h,(i+q)*h)
            require(abs(float(exact)-left[i]) <= var[i]+err, 'angular cell range')
            if q == 0 or q == Q(1,2):
                center = left[i] if q == 0 else mid[i]
                require(abs(float(exact)-center) <= err, 'angular knot or midpoint')
    require(abs(-left[0]-float(angular_fixture(a,b,h,0))) > err,'odd-sign mutant survived')
    require(abs(left[2]-float(angular_fixture(a,b,h,Q(5,16)))) > err,'missing midpoint mutant survived')
    require(abs(mid[3]*2/3-float(angular_fixture(a,b,h,Q(7,16)))) > err,'wrong radial factor survived')
    require(abs(left[6]-float(angular_fixture(a,b,h,Q(7,8)))) > err,'shifted support mutant survived')
    print('SPATIAL exact angular fixtures=45; mutants rejected=4',flush=True)


def dynamic_budget(lam,x):
    # Static and all finite Einstein ratios eta >= 2^18.
    perturb = x/2**17
    a = 2*lam/(Q('3.14159')*x)
    delta = 1/(12*Q('0.18272624777232950')**2 * lam*2**36)
    remainder = a*delta*2 + VOLUME*a*a*(2*delta+delta*delta)*8
    lip = a*delta*4 + 3*VOLUME*a*a*(2*delta+delta*delta)*16
    require(remainder < perturb/x, 'nonstatic perturbation budget')
    return perturb, lip


def main():
    global INDEPENDENT
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xs', type=Q, nargs='+', default=[Q('1/4'),Q('1/2'),Q(1),Q(2)])
    p.add_argument('--lam', type=Q, nargs='+', default=[Q(2)])
    p.add_argument('--cells', type=int, default=1024)
    p.add_argument('--steps', type=int, default=40)
    p.add_argument('--independent', action='store_true')
    args = p.parse_args()
    INDEPENDENT = args.independent
    spatial_fixtures()
    require(64 <= args.cells <= 4096 and args.cells & (args.cells-1) == 0, 'dyadic cell limit')
    require(args.steps > 0 and len(set(args.lam)) == len(args.lam)
            and len(set(args.xs)) == len(args.xs), 'positive steps and unique grid coordinates')
    full = {}
    print('lambda,x,w,static_upper_or_status,derivative_upper',flush=True)
    for lam in args.lam:
        for x in args.xs:
            require(x in (Q(2)**i for i in range(-4,5)) and 0 < lam <= 2,'certified input ranges')
            budget, lip = dynamic_budget(lam,x)
            lip_f = up(float(lip))
            upper = 0.
            for w in range(1,101,2):
                a = 2*(lam/x)/(Q('3.14159')*w*w)
                h = Q('1/4')
                f,b = guards(1/(1-h),1/(1-h)**2,a,h,VOLUME)
                if f + budget/(x*w) < h and b + lip < 1:
                    upper = max(upper,up(float(b)+lip_f))
                    print(f'{lam},{x},>={w},{decimal(upper,8)},HOMOTOPY',flush=True)
                    break
                cells = max(args.cells,4096 if x <= Q(1,8) else 1024)
                result, reason = certify(x*w,lam/w,cells,args.steps,float(budget),lip_f)
                if result is None:
                    print(f'{lam},{x},{w},UNRESOLVED {reason}',flush=True)
                    raise ArithmeticError(f'uncertified static grid point: lambda={lam}, x={x}, w={w}')
                ev,q = result
                upper = max(upper,ev)
                print(f'{lam},{x},{w},{decimal(ev,8)},{decimal(q,8)}',flush=True)
            else:
                raise ArithmeticError('unclosed static Matsubara tail')
            full[lam,x] = upper
            print(f'{lam},{x},FULL,{decimal(upper,8)}',flush=True)
    require(len(full) == len(args.lam)*len(args.xs), 'not every requested point reached FULL')
    peak = max(full.values())
    if set(args.lam) == set(LAM) and set(args.xs) == set(GRID):
        require(len(full) == 45 and Q(peak) < Q('0.70431955'), 'published static grid target')
    print(f'STATIC PASS FULL={len(full)} max<={decimal(peak,8)}',flush=True)


if __name__ == '__main__':
    main()
