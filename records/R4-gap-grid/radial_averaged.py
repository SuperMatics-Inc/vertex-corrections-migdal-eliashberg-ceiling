"""Tighter Krawczyk residuals using cell averages of r times the spatial line."""
from radial import main, up
from radial_check import HalfLine, verify
from radial_moments import MomentIndependent, MomentLines, MomentReplay
from radial_space import Lines, Space
import numpy as np


class AveragedHeights:
    mean_height = True
    inexact_convolution = False

    def state(self, sigma, widths):
        base = super().state(sigma, widths)
        w = widths[:,None]
        v = 1j*w-sigma
        h,r = self.h,self.r
        lo,hi = r-h/2,self.right
        y0,y1 = (2*lo**3-1)/self.x,(2*hi**3-1)/self.x
        d0,d,d1 = v-y0,v-self.xi,v-y1
        height = (lo/d0+4*r/d+hi/d1)/6
        p,a,b,c = base.peak,self.xi_d,self.xi_dd,12/self.x
        L = up(p+hi*p*p*a)
        Z3 = up(6*p**4*a**3+6*p**3*a*b+p*p*c)
        Z4 = up(24*p**5*a**4+36*p**4*a*a*b+6*p**3*b*b+8*p**3*a*c)
        scale = (1+w+abs(sigma)+abs(self.xi)+abs(y0)+abs(y1))*p*p+abs(height)
        cell_error = up(h**4/720*(4*Z3+hi*Z4),scale)
        # Complex variance identity around the true average, and average of
        # |r-t| on one cell, before adding the enclosure of that average.
        e2 = up(h/np.sqrt(12)*self.norm(L)+self.norm(cell_error))
        e1 = up(3*h*np.sum(hi*(h*L/3+cell_error),axis=-1))
        proxy = Lines(self,height,p,e1,e2,base.p2,(base.mean,base.mean_err),
                      (base.sqavg,base.sqavg_err),base.width)
        return MomentLines(proxy,d,self.inexact_convolution,cell_error)


class AveragedDirect(AveragedHeights,Space):
    pass


class AveragedHalfLine(AveragedHeights,HalfLine):
    inexact_convolution = True


class AveragedReplay(MomentReplay):
    def __init__(self,lam,eta,x,n,cells,panel):
        super().__init__(lam,eta,x,n,cells,panel)
        self.space = AveragedDirect(self.x,cells)


class AveragedIndependent(MomentIndependent):
    def __init__(self,lam,eta,x,n,cells,panel):
        super().__init__(lam,eta,x,n,cells,panel)
        self.space = AveragedHalfLine(self.x,cells)


def second(lam,eta,x,primary,witness):
    return verify(lam,eta,x,primary,witness,proto=AveragedIndependent)


if __name__ == '__main__':
    main(AveragedReplay,second)
