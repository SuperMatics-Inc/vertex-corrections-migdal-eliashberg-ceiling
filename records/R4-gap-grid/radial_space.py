"""Continuous radial convolution of bounded step heights and interval remainders."""
import numpy as np

from replay import require

# A decimal pad is not used to account for cancellation. The calls enclosing signed
# quantities take a separate absolute scale from before the cancellation.
PAD = 2.**-30


def up(a, scale=None):
    if scale is None:
        scale = abs(a)
    return np.nextafter(a + PAD*(1+scale), np.inf)


def down(a, scale=None):
    if scale is None:
        scale = abs(a)
    return np.nextafter(a - PAD*(1+scale), -np.inf)


class Space:
    """H_f(r)=r*f(r) on a dyadic grid; step approximation is applied to H_f."""
    mean_height = False

    def __init__(self, x, cells):
        self.x = x
        self.cells = cells
        self.h = 1/cells
        self.r = (np.arange(cells)+.5)*self.h
        self.right = self.r+self.h/2
        self.xi = (2*self.r**3-1)/x
        self.xi_dev = 3*self.right**2*self.h/x
        self.xi_d = 6*self.right**2/x
        self.xi_dd = 12*self.right/x

    def norm(self, a):
        return up(np.sqrt(3*self.h*np.sum(abs(a)**2,axis=-1)))

    def convolve(self, a, b):
        """Piecewise quadratic primitive of the correlation of full odd arrays."""
        aa,bb = np.r_[-a[::-1],a],np.r_[-b[::-1],b]
        D = -1.5*self.h*np.r_[np.convolve(aa,bb)[len(a)+len(b)-1:],0.]
        return self.primitive(D,up(6*self.h*min(len(a),len(b))*max(abs(a))*max(abs(b))))

    def primitive(self, d, scale):
        """Returns representative heights, midpoint values and a cell-uniform error."""
        err = up(0.,scale)
        k = np.arange(len(d)-1)
        h = self.h
        a = np.r_[0.,np.cumsum((d[:-1]+d[1:])*(h/2))][:-1]
        mid = a + h*(3*d[:-1]+d[1:])/8
        val = mid/(h*(k+.5))
        # variation about the left value, then translate to the midpoint representative
        v = np.r_[d[0],a[1:]/(h*k[1:])]
        rv = up((np.maximum(abs(d[:-1]-v),abs(d[1:]-v))/(k+1)) + abs(v-val) + 4*err, scale)
        variation = up(np.maximum(abs(d[:-1]),abs(d[1:]))+err,scale)
        e2 = up((h/np.sqrt(12))*self.norm(variation)+(len(d)-1)*h*np.sqrt(3*(len(d)-1)*h)*4*err,scale)
        height = a + h*(2*d[:-1]+d[1:])/6 if self.mean_height else mid
        return height, val, rv, e2

    def state(self, sigma, widths):
        """Enclose the center line, unsquared and squared averages, and step errors."""
        w = widths[:,None]
        den = 1j*w-self.xi-sigma
        width = down(den.imag, w+abs(sigma))
        require(np.min(width/w) > .25,'nonpositive radial center widths')
        amp = up(abs(den),w+abs(self.xi)+abs(sigma))
        hor = np.maximum(0.,down(abs(den.real)-up(self.xi_dev),abs(self.xi)+abs(sigma)))
        peak = up(1/down(np.hypot(hor,width)))
        z = 1/den
        zm = up(abs(z),amp)
        h0 = self.r*z
        slope = up(zm+self.right*zm*peak*self.xi_d)
        # Cover midpoint complex division, center addition and dispersion arithmetic.
        arithmetic = up(0.,(1+amp)*peak**2) + up(0.,abs(h0))
        e2 = up(self.h/np.sqrt(12)*self.norm(slope)+self.norm(arithmetic))
        e1 = up(3*self.h*np.sum(self.right*(self.h/4*slope+arithmetic),axis=-1))
        p2 = up((peak**2) @ (np.diff(np.r_[0.,self.right]**3)))
        # Midpoint integral of r^2 Z^p. Z'=Z^2 xi', Z''=2Z^3 xi'^2+Z^2 xi''.
        rs = self.right
        p = peak
        def integrate(degree):
            f = z**degree
            fp = degree*p**(degree+1)*self.xi_d
            fpp = degree*(degree+1)*p**(degree+2)*self.xi_d**2+degree*p**(degree+1)*self.xi_dd
            rem = 3*self.h**3/24*np.sum(2*p**degree+4*rs*fp+rs**2*fpp,axis=-1)
            mean = 3*self.h*np.sum(self.r**2*f,axis=-1)
            scale = 3*self.h*np.sum(rs**2*(p**degree+degree*p**(degree-1)*arithmetic),axis=-1)
            return mean,up(rem+3*self.h*np.sum(rs**2*degree*p**(degree-1)*arithmetic,axis=-1),scale)
        return Lines(self,h0,peak,e1,e2,p2,integrate(1),integrate(2),min(np.min(width/w),1.))


class Lines:
    def __init__(self, mesh, h, peak, e1, e2, p2, mean, sqavg, width):
        self.mesh = mesh
        self.h = h
        self.peak = peak
        self.p2 = p2
        self.a2 = mesh.norm(h)
        self.a1 = up(3*mesh.h*np.sum(mesh.right*abs(h),axis=-1))
        self.norm1 = up(np.minimum(np.sqrt(p2),self.a1+e1))
        self.norm2 = up(np.minimum(np.sqrt(p2),self.a2+e2))
        self.e1,self.e2 = e1,e2
        self.mean,self.mean_err = mean
        self.sqavg,self.sqavg_err = sqavg
        self.width = width
        self.cache = {}

    @staticmethod
    def fold(i):
        return i if i >= 0 else -1-i

    def value(self, index):
        f = self.fold(index)
        return self.h[f] if index >= 0 else self.h[f].conj()

    def pair(self,i,j):
        key = tuple(sorted((i,j)))
        if key in self.cache:
            return self.cache[key]
        flip = tuple(sorted((-1-i,-1-j)))
        if flip in self.cache:
            h,c,*error = self.cache[flip]
            return (h.conj(),c.conj(),*error)
        a,b = self.fold(i),self.fold(j)
        h,c,r,e = self.mesh.convolve(self.value(i),self.value(j))
        pair_sup = up(min(self.e2[a]*self.norm2[b]+self.a2[a]*self.e2[b],
                          self.e2[b]*self.norm2[a]+self.a2[b]*self.e2[a]))
        pair_l2 = up(min(self.e1[a]*self.norm2[b]+self.a1[a]*self.e2[b],
                         self.e1[b]*self.norm2[a]+self.a1[b]*self.e2[a],
                         self.e2[a]*self.norm1[b]+self.a2[a]*self.e1[b],
                         self.e2[b]*self.norm1[a]+self.a2[b]*self.e1[a]))
        ans = (h,c,up(r+pair_sup),up(e+pair_l2))
        self.cache[key] = ans
        return ans
