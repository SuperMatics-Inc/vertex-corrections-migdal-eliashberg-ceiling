"""Use cell-mean cancellation in the normal residual and signed radial kernels."""
from typing import NamedTuple

from radial import Replay, fold, main, up
from radial_check import HalfLine, Independent, verify
from radial_space import Lines, Space
import numpy as np


class Moments(NamedTuple):
    slope: float
    integral_error: float  # sum of absolute integrals of (true - represented) over the full cells
    variation: float       # total variation of the odd, zero-extended step representative
    peak: float            # peak of represented heights, NOT peak of the radial function
    error2: float          # radial L2 of the difference from the representative


def variation(h):
    return up(2*(abs(h[...,0])+np.sum(abs(np.diff(h,axis=-1)),axis=-1)+abs(h[...,-1])),
              8*np.sum(abs(h),axis=-1))


def correlation_error(h,a,b):
    return up(h*h/8*(a.slope*b.variation+b.slope*a.variation)
              +a.integral_error*b.peak+b.integral_error*a.peak+(2/3)*a.error2*b.error2)


class MomentLines(Lines):
    def __init__(self, base, den, inexact_convolution, cell_error=None):
        self.__dict__.update(vars(base))
        self.curves = {}
        self.inexact_convolution = inexact_convolution
        mesh = self.mesh
        p = self.peak
        L = up(p+mesh.right*p**2*mesh.xi_d)
        if cell_error is None:
            L2 = up(2*p**2*mesh.xi_d+mesh.right*(2*p**3*mesh.xi_d**2+p**2*mesh.xi_dd))
            rounding = up(0.,(1+abs(den))*p**2)+up(0.,abs(self.h))
            cell_error = up(mesh.h**2/24*L2+rounding)
        integrated = up(2*mesh.h*np.sum(cell_error,axis=-1))
        self.moments = tuple(Moments(*row) for row in zip(up(np.max(L,axis=-1)),integrated,
                             variation(self.h),up(np.max(abs(self.h),axis=-1)),self.e2))

    def pair(self,i,j):
        key = tuple(sorted((i,j)))
        if key in self.cache:
            return self.cache[key]
        ref = tuple(sorted((-1-i,-1-j)))
        if ref in self.cache:
            h,c,*error = self.cache[ref]
            return h.conj(),c.conj(),*error
        a,b = self.fold(i),self.fold(j)
        ma,mb = self.moments[a],self.moments[b]
        mesh = self.mesh
        h,c,r,e = mesh.convolve(self.value(i),self.value(j))
        sup = up(1.5*correlation_error(mesh.h,ma,mb))
        su = min(sup,up(min(self.e2[a]*self.norm2[b]+self.a2[a]*self.e2[b],
                            self.e2[b]*self.norm2[a]+self.a2[b]*self.e2[a])))
        err = min(up(np.sqrt(8)*sup),up(min(self.e1[a]*self.norm2[b]+self.a1[a]*self.e2[b],
                            self.e1[b]*self.norm2[a]+self.a1[b]*self.e2[a],
                            self.e2[a]*self.norm1[b]+self.a2[a]*self.e1[b],
                            self.e2[b]*self.norm1[a]+self.a2[b]*self.e1[a])))
        grad = up(self.a2[a]*self.a2[b])
        curvature = 0. if mesh.mean_height else up(1.5*min(ma.variation*mb.peak,mb.variation*ma.peak))
        sc = up(6*mesh.h*mesh.cells*ma.peak*mb.peak)
        derr = up(0.,sc)
        if self.inexact_convolution:
            # Same quantization scale as the integer verifier, evaluated on its
            # actual operand peaks. Use the unpadded peaks for the branch only.
            mx = max(1.,max(abs(self.h[a])),max(abs(self.h[b])))
            unit = 2.**(-min(24,int(np.floor(np.log2(2.**25/mx)))))
            derr += up(3*mesh.h*mesh.cells*unit*(ma.peak+mb.peak+unit),sc)
        mean = up(4*(mesh.h**2/24*curvature+8*derr))
        self.curves[key] = (grad,mean,e,err)
        ans = h,c,up(r+su),up(e+err)
        self.cache[key] = ans
        return ans

    def curve(self,i,j):
        key = tuple(sorted((i,j)))
        ref = tuple(sorted((-1-i,-1-j)))
        return self.curves[key] if key in self.curves else self.curves[ref]


class MomentSpace(Space):
    def state(self,sigma,widths):
        return MomentLines(super().state(sigma,widths),1j*widths[:,None]-self.xi-sigma,False)


class MomentHalfLine(HalfLine):
    def state(self,sigma,widths):
        return MomentLines(super().state(sigma,widths),1j*widths[:,None]-self.xi-sigma,True)


class CubicMoments:
    def setup_transfers(self):
        self.tight = min(self.outer//2,4*self.j)
        self.tight_rows = min(6,self.n)
        s = np.arange(-self.tight,self.tight+1)
        n = np.arange(self.tight_rows)[:,None,None]
        a = fold(n-s[None,:,None])
        b = fold(n-s[None,:,None]-s[None,None,:])
        c = fold(n-s[None,None,:])
        g = self.ell_hi[self.outer+s]
        self.tight_idx = (a,b,c)
        self.tight_weight = up(g[:,None]*g[None,:])
        self.tight_wing = (abs(s[:,None])>self.j)|(abs(s[None,:])>self.j)

    def envelopes(self,line):
        r,rD,G = super().envelopes(line)
        q = self.tight_rows
        p = np.max(line.peak,axis=1)
        L = (line.norm1,line.norm2,p)
        # The relative resolvent excursion of line m is at most e * rel[m].
        rel = np.minimum(1.,up(line.width*self.w*p))
        scales = up(1/line.width)
        idx = self.idx[:q]
        near = abs(np.arange(-self.outer,self.outer+1))<=self.tight

        def parts(profiles,cap):
            inc = [up(v[idx[:,near]]@self.ell_hi[near]) for v in profiles]
            out = [up(v[idx[:,~near]]@self.ell_hi[~near]
                      +self.missing*cap/(2*self.outer-2*self.cut+3)) for v in profiles]
            return inc,out

        def outer(A,C,profileB):
            def w(i,j):
                return self.wings((A[0][i],A[1][i]),(C[0][j],C[1][j]))
            M = [max(v) for v in profileB]
            return up(np.minimum.reduce((M[0]*w(1,1),M[1]*w(0,1),
                                         M[1]*w(1,0),(15/32)*M[2]*w(2,2))))

        def inner(A,B,C,variation_of=None):
            a,b,c = self.tight_idx
            e = np.minimum.reduce((A[0][a]*B[1][b]*C[1][c],A[1][a]*B[0][b]*C[1][c],
                                   A[1][a]*B[1][b]*C[0][c],(15/32)*A[2][a]*B[2][b]*C[2][c]))
            if variation_of is not None:
                e *= up((rel[variation_of[0]]+rel[variation_of[1]])/2)
            e = up(e*self.tight_weight)
            return up(np.sum(e,axis=(1,2))),up(np.sum(e*self.tight_wing,axis=(1,2)))

        def tighten(a,b,c,ap,cp,factor):
            i,j = inner(a,b,c)
            o = outer(ap,cp,b)
            return [up(factor*(h+o)) for h in (i,j)]

        profile = parts(L,scales)
        C = [tighten(L,L,L,profile,profile,1.)]
        d_rainbow = self.sums(up(self.w*line.p2*rel),up(scales**2))
        growth = [up(d_rainbow[0]+d_rainbow[1])]
        for scale in (self.w,1.):
            D = (scale*line.p2,scale*p*line.norm2,scale*p*p)
            enriched = parts(D,up(scales**2))
            end = tighten(D,L,L,enriched,profile,2.)
            mid = tighten(L,D,L,profile,profile,1.)
            C.append([up(a+b) for a,b in zip(end,mid)])
            a,b,c = self.tight_idx
            iv_end,_ = inner(D,L,L,(b,c))
            iv_mid,_ = inner(L,D,L,(a,c))
            outer_end = outer(enriched,profile,L)
            outer_mid = outer(profile,profile,D)
            absolute = (G[0][1] if np.ndim(scale) else G[0][2]).copy()
            absolute[:q] = np.minimum(absolute[:q],up(2*(iv_end+outer_end)+iv_mid+outer_mid))
            growth.append(absolute)
        for order in range(2):
            for slot in range(3):
                G[order][slot][:q] = np.minimum(G[order][slot][:q],C[slot][order])
        line.growth = tuple(growth)
        return r,rD,G

    def cubic_error(self,line,n,m,h,e):
        old = super().cubic_error(line,n,m,h,e)
        slope = me = ms = er = 0.
        scale = 0.
        for t in range(-self.j,self.j+1):
            g,ge = self.g(t),self.ge(t)
            hi = self.ell_hi[self.outer+t]
            a,b,c,d = line.curve(n-t,m-t)
            s,t2 = int(fold(n-t)),int(fold(m-t))
            slope += g*a
            me += g*b
            ms += g*c
            # The exact positive-weight sum minus the sum of convolved
            # nominal step lines. Include coefficient mismatch in this term.
            er += g*d+ge*min(line.norm1[s]*line.norm2[t2],line.norm2[s]*line.norm1[t2])
            scale += hi*(1+a+b+c+d)
        p = Moments(up(slope,scale),up(me,scale),variation(h),up(max(abs(h))),up(ms,scale))
        f = int(fold(m))
        return min(old,up(up(er,scale)*line.norm2[f]+1.5*correlation_error(self.space.h,p,line.moments[f]),scale))


class MomentReplay(CubicMoments,Replay):
    def __init__(self,lam,eta,x,n,cells,panel):
        super().__init__(lam,eta,x,n,cells,panel,MomentSpace)
        self.setup_transfers()
        self.spatial = min(n,16)


class MomentIndependent(CubicMoments,Independent):
    def __init__(self,lam,eta,x,n,cells,panel):
        super().__init__(lam,eta,x,n,cells,panel)
        self.space = MomentHalfLine(self.x,cells)
        self.setup_transfers()
        self.spatial = min(n,16)


def second(lam,eta,x,primary,witness):
    return verify(lam,eta,x,primary,witness,proto=MomentIndependent)


if __name__ == '__main__':
    main(MomentReplay,second)
