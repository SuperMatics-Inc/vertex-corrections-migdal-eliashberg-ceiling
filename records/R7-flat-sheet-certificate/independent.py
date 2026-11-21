#!/usr/bin/env python3
"""independent.py -- second, independently written check of the VERTEX-UNIV 1D
flat-sheet, vertex-corrected counterexample certificate (cert.py).

Independence from cert.py
-------------------------
* arithmetic backend: MPFR with *explicit directed rounding* (gmpy2.RoundDown /
  RoundUp, 160-bit), wrapped in a small real/complex interval layer (RI/CI) written
  here for every model quantity.  The one exception is the signed dense matrix
  product of the contraction enclosure (`krawczyk_modulus`), which is carried out in
  `flint` ball arithmetic (`acb_mat`, 220 bits) after an exact, outward-padded
  interval-to-ball embedding of each RI/CI entry; its result is re-rounded upward
  into an mpfr endpoint.  cert.py uses FLINT's `arb`/`acb` (python-flint) throughout
  instead.  The RI/CI backend is validated against exact rational
  (`fractions.Fraction`) arithmetic by `--selftest`.
* algebraic arrangement: the correlation sums are kept as shifted tables in the
  *momentum* index,
      V1[p][q][s] = sum_m z(p,m) z(q,(m+s) mod Nk)      (crossed / rainbow-pair),
      V2[p][q][s] = sum_m z(p,m) z(q,(s-m) mod Nk),
      Q1[p][q][s] = sum_m z(p,m) conj(z(q,(s-m) mod Nk)),
  and the Jacobian is obtained from a four-term expansion in this arrangement,
  instead of the (K(n,a)K(n,c) table-product) form used in cert.py.
* contraction enclosure: the *signed* product E = I - A (I - DPhi(box)) is formed
  explicitly in python-flint `acb_mat` arithmetic (a triangle-inequality bound
  || |A| W ||_inf would be >= 1 and useless), with the row sup-norm taken by
  overlap-safe `arb.max` reductions; the entrywise boxes come from this file's own
  gmpy2 directed-rounding intervals, not from cert.py's arb tables.
* fixed point: Banach on T(S) = S - A(Phi(S) - S) (a contraction because k < 1),
  so existence and uniqueness of Sigma* in the box follow without any contraction of
  Phi itself (the phase-blind ||DPhi(Sigma_0)||_inf is about 2.96 >> 1).
  A = inv(I - DPhi(Sigma_0)) is a float matrix lifted to exact doubles; its
  nonsingularity is *proved* by k < 1 (if A x = 0 with x != 0 then
  x = (I - A(I - DPhi)) x, so ||x|| <= k ||x|| < ||x||, a contradiction).
* spectrum: with S(Sigma) = (B(Sigma) + B(Sigma)^T)/2 the symmetric part of the
  linearised gap operator, lambda_max(S(Sigma)) >= <v, S v>/<v,v> for every Sigma in
  the box and every real v != 0 (Rayleigh quotient of a real symmetric matrix), and
  the quotient is bounded below by the *lower endpoint* of the same expression
  evaluated on the entrywise enclosure of B over the box.  No antisymmetry penalty is
  subtracted (the v1 version's "- ||B-B^T||_inf/2" correction was not a valid
  inequality); the enclosure itself is symmetrised, so the bound covers every matrix
  of the entrywise ball with that symmetric part.

MODEL (identical to cert.py; exact rational / exact-decimal constants)
----------------------------------------------------------------------
M(Nk, Nw; lam, m, eta, mu0) -- strictly 1D two-sheet band
    xi_j = 2 |2j+1-Nk| / Nk - 1 - mu0 ,  j = 0..Nk-1            (mu0 = -2/5)
    w_n = 2n + 1 , n = -Nw..Nw-1 , M = 2 Nw                     (Matsubara)
    tau = eta C_inf sqrt(lam) m , nu = pi tau                    (C_inf exact decimal)
    K(d) = m^2 / (m^2 + (2 pi tau d)^2) ,  d = n - a             (phonon kernel)
    gRB = 2 lam tau ,  CX = (2 lam tau / Nk)^2                   (couplings)
    z(n,j) = 1 / (i w_n nu - xi_j - Sigma(n,j))
    Phi(Sigma)(n,j) = gRB sum_a K(n,a) <z(a)> + CX sum_{a,c} K(n,a) K(n,c) x
                      x [z(a,j) z(c,.) ... ]        (rainbow + first-crossed skeleton)
The normal-state map is exactly the one implemented in
VERTEX-UNIV/verify_1d_counterexample.py:solve_finite_m_lattice (audited artifact),
transcribed in the V1/V2 arrangement above; the transcription is checked numerically
against that artifact ([C]) and against central differences of Phi ([B]).

Usage:  python independent.py [--Nk 24 --Nw 12 --lam 8 --m 1/100 ...]
        python independent.py --selftest     (backend validation, exact rationals)
Exit 0 iff every certified inequality holds.
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from fractions import Fraction

import gmpy2
import numpy as np

PREC = 160                     # bits of every directed-rounding operation
C_INF = "0.18272624777232952"  # exact decimal (same literal as the audited artifact)
THETA = "5e-3"                 # certified margin: lambda_max >= 1 + THETA
RHO = "1e-12"                  # sup-norm radius of the fixed-point box

# Guarded mutants: each one must make this checker exit 1.  The names and the nature
# of each fault mirror the mutants of cert.py, so that every fault is required to be
# caught by *both* checkers.
MUTANTS = (
    "no_vertex",        # drop the vertex correction from the model itself
    "m_threshold",      # move to m = 2/100 (below the crossing)
    "rho_zero",         # verify a point instead of a box
    "res_scale",        # inflate the residual enclosure by 2^20
    "rb_scale",         # inconsistent model: rainbow kernel scaled by 6/5 in the box only
    "ray_thresh",       # demand lambda_max >= 1.01 instead of 1.005
    "diag_drop",        # under-report the diagonal of B by 2e-2 (Rayleigh bound fails)
    "drop_lb_t3",       # drop the vertex term from the linearised operator only
    "point_not_ball",   # collapse the box to a point while rho stays 1e-12
    "float_leak",       # replace the exact interval 1/Nk by the binary64 lift
)

RoundDown, RoundUp = gmpy2.RoundDown, gmpy2.RoundUp


# ==================================================================================
#  directed-rounding real interval [a, b]
# ==================================================================================
def _dn():
    return gmpy2.context(precision=PREC, round=RoundDown)


def _up():
    return gmpy2.context(precision=PREC, round=RoundUp)


def _neg_exact(x):
    """exact negation of an mpfr endpoint (gmpy2's unary minus rounds to the
    *current* context precision, so it must be evaluated at full precision)."""
    with gmpy2.context(precision=PREC, round=gmpy2.RoundToNearest):
        return -x


def _f(x):
    """float view of an mpfr endpoint (for printing only)."""
    return float(x)


def _fmt(d, sig, up):
    """sig significant digits of the float d, rounded down (up False) or up (up True)."""
    from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
    if d == 0 or not math.isfinite(d):
        return repr(d)
    e = math.floor(math.log10(abs(d)))
    places = max(0, sig - 1 - e)
    return str(Decimal(d).quantize(Decimal(1).scaleb(-places),
                                   rounding=ROUND_CEILING if up else ROUND_FLOOR))


def _dlo(x, sig=12):
    """Decimal string that is <= the exact mpfr endpoint x."""
    d = _f(x)
    if not math.isfinite(d):
        return repr(d)
    return _fmt(d - 8.0 * abs(math.ulp(d)), sig, False)


def _dup(x, sig=12):
    """Decimal string that is >= the exact mpfr endpoint x."""
    d = _f(x)
    if not math.isfinite(d):
        return repr(d)
    return _fmt(d + 8.0 * abs(math.ulp(d)), sig, True)


class RI:
    """Closed real interval [a, b] with directed-rounding mpfr endpoints."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=None):
        self.a = a
        self.b = a if b is None else b
        if self.a > self.b:                                   # pragma: no cover
            raise ValueError("empty interval")

    # ---- constructors -------------------------------------------------------
    @staticmethod
    def exact(x):
        """Exact lift of an int / float / mpz / mpq / mpfr value."""
        if isinstance(x, RI):
            return x
        if isinstance(x, float):
            v = gmpy2.mpfr(x)                                 # doubles are exact in mpfr
            return RI(v, v)
        if isinstance(x, bool):
            return RI(gmpy2.mpfr(int(x)), gmpy2.mpfr(int(x)))
        if isinstance(x, int):
            # gmpy2.mpfr(int) is exact only below the working precision; integers
            # beyond 2^53 are lifted through mpq with directed rounding instead, so
            # that no big integer can be silently rounded.
            if abs(x) < 2 ** 53:
                v = gmpy2.mpfr(x)
                return RI(v, v)
            with _dn():
                lo = gmpy2.mpfr(gmpy2.mpq(x))
            with _up():
                hi = gmpy2.mpfr(gmpy2.mpq(x))
            return RI(lo, hi)
        if isinstance(x, gmpy2.mpfr):
            return RI(x, x)
        if isinstance(x, gmpy2.mpq):
            with _dn():
                lo = gmpy2.mpfr(x)
            with _up():
                hi = gmpy2.mpfr(x)
            return RI(lo, hi)
        raise TypeError("RI.exact: unsupported %r" % (type(x),))

    @staticmethod
    def decimal(s):
        """Enclosure of an exact decimal literal such as '1e-12' or '0.182726'."""
        with _dn():
            lo = gmpy2.mpfr(gmpy2.mpq(str(s)))
        with _up():
            hi = gmpy2.mpfr(gmpy2.mpq(str(s)))
        return RI(lo, hi)

    @staticmethod
    def rational(p, q=1):
        with _dn():
            lo = gmpy2.mpfr(gmpy2.mpq(int(p)) / gmpy2.mpq(int(q)))
        with _up():
            hi = gmpy2.mpfr(gmpy2.mpq(int(p)) / gmpy2.mpq(int(q)))
        return RI(lo, hi)

    # ---- basic ------------------------------------------------------------
    def __repr__(self):
        return "[%.17g, %.17g]" % (_f(self.a), _f(self.b))

    def __add__(self, o):
        if isinstance(o, CI):
            return NotImplemented
        o = RI.exact(o)
        with _dn():
            lo = self.a + o.a
        with _up():
            hi = self.b + o.b
        return RI(lo, hi)

    __radd__ = __add__

    def __neg__(self):
        return RI(_neg_exact(self.b), _neg_exact(self.a))

    def __sub__(self, o):
        if isinstance(o, CI):
            return NotImplemented
        return self + (-RI.exact(o))

    def __rsub__(self, o):
        return RI.exact(o) - self

    def __mul__(self, o):
        if isinstance(o, CI):
            return NotImplemented
        o = RI.exact(o)
        with _dn():
            p_lo = [self.a * o.a, self.a * o.b, self.b * o.a, self.b * o.b]
        with _up():
            p_hi = [self.a * o.a, self.a * o.b, self.b * o.a, self.b * o.b]
        return RI(min(p_lo), max(p_hi))

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, CI):
            return NotImplemented
        o = RI.exact(o)
        return self * o.recip()

    def __rtruediv__(self, o):
        return RI.exact(o) * self.recip()

    def recip(self):
        if self.a > 0:
            with _dn():
                lo = gmpy2.mpfr(1) / self.b
            with _up():
                hi = gmpy2.mpfr(1) / self.a
            return RI(lo, hi)
        if self.b < 0:
            with _dn():
                lo = gmpy2.mpfr(1) / self.b
            with _up():
                hi = gmpy2.mpfr(1) / self.a
            return RI(lo, hi)
        raise ZeroDivisionError("interval divisor contains 0: %r" % (self,))

    def sqrt(self):
        assert self.a >= 0
        with _dn():
            lo = gmpy2.sqrt(self.a) if self.a > 0 else gmpy2.mpfr(0)
        with _up():
            hi = gmpy2.sqrt(self.b)
        return RI(lo, hi)

    def __pow__(self, k):
        assert isinstance(k, int) and k > 0
        r = self
        for _ in range(k - 1):
            r = r * self
        return r

    # ---- tests (certainty semantics, like mpmath.iv) ----------------------
    def __lt__(self, o):
        o = RI.exact(o)
        return bool(self.b < o.a)

    def __gt__(self, o):
        o = RI.exact(o)
        return bool(self.a > o.b)

    def __le__(self, o):
        o = RI.exact(o)
        return bool(self.b <= o.a)

    def __ge__(self, o):
        o = RI.exact(o)
        return bool(self.a >= o.b)

    def contains(self, x):
        x = RI.exact(x)
        return bool(self.a <= x.a and x.b <= self.b)

    def mid(self):
        with gmpy2.context(precision=PREC, round=gmpy2.RoundToNearest):
            return (self.a + self.b) / 2

    def __float__(self):
        return _f(self.mid())

    def __abs__(self):
        return self.abs_iv()

    def abs_iv(self):
        if self.a >= 0:
            return self
        if self.b <= 0:
            return -self
        return RI(gmpy2.mpfr(0), max(_neg_exact(self.a), self.b))


# ==================================================================================
#  complex interval box
# ==================================================================================
class CI:
    """Axis-aligned complex interval (re x im box) built from RI components."""

    __slots__ = ("re", "im")

    def __init__(self, re, im=None):
        self.re = RI.exact(re) if not isinstance(re, RI) else re
        self.im = RI.exact(0 if im is None else im) if not isinstance(im, RI) else im
        if im is None:
            self.im = RI.exact(0)

    def __repr__(self):
        return "(%r + %rj)" % (self.re, self.im)

    def __add__(self, o):
        o = promote(o)
        return CI(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    def __neg__(self):
        return CI(-self.re, -self.im)

    def __sub__(self, o):
        return self + (-promote(o))

    def __rsub__(self, o):
        return promote(o) - self

    def __mul__(self, o):
        o = promote(o)
        return CI(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = promote(o)
        d = o.re * o.re + o.im * o.im
        return CI((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)

    def __rtruediv__(self, o):
        return promote(o) / self

    def conjugate(self):
        return CI(self.re, -self.im)

    def abs_iv(self):
        """enclosure of |z| (lower endpoint clamped at 0 to absorb product slack)."""
        sq = self.re * self.re + self.im * self.im
        return RI(sq.a if sq.a > 0 else gmpy2.mpfr(0), sq.b).sqrt()


def promote(x):
    if isinstance(x, CI):
        return x
    if isinstance(x, RI):
        return CI(x)
    if isinstance(x, np.complexfloating):
        x = complex(x)
    if isinstance(x, complex):
        return CI(RI.exact(float(x.real)), RI.exact(float(x.imag)))
    if isinstance(x, np.integer):
        x = int(x)
    if isinstance(x, np.floating):
        x = float(x)
    return CI(RI.exact(x))


def conj(x):
    return promote(x).conjugate()


def cmod(x):
    """|x| as an RI in the interval path, and float|complex abs in the float path."""
    if isinstance(x, (CI, RI)):
        return promote(x).abs_iv()
    return abs(x)


def _zero_like(z):
    z0 = z[0][0]
    return promote(0) if isinstance(z0, (CI, RI)) else 0j


def real_iv(v):
    return promote(v).re


def iv_mpf(x):
    """scalar interval from an int / float / decimal string / 1-element list."""
    if isinstance(x, RI):
        return x
    if isinstance(x, (list, tuple)):
        vals = [iv_mpf(v) for v in x]
        lo = min(v.a for v in vals)
        hi = max(v.b for v in vals)
        return RI(lo, hi)
    if isinstance(x, str):
        return RI.decimal(x)
    return RI.exact(x)


def iv_pi():
    with _dn():
        lo = gmpy2.const_pi()
    with _up():
        hi = gmpy2.const_pi()
    return RI(lo, hi)


# ==================================================================================
#  model formulas (identical code driven by float complex or by CI intervals)
# ==================================================================================
def corr_tables(z, M, Nk):
    """V1[p][q][s] = sum_m z(p,m) z(q,m+s) ; V2 = sum_m z(p,m) z(q,s-m) ;
       Q1[p][q][s] = sum_m z(p,m) conj(z(q,s-m))   (indices mod Nk)."""
    zero = _zero_like(z)
    V1, V2, Q1 = [], [], []
    for p in range(M):
        r1, r2, r3 = [], [], []
        for q in range(M):
            a1 = [zero] * Nk
            a2 = [zero] * Nk
            a3 = [zero] * Nk
            for m in range(Nk):
                zm = z[p][m]
                for s in range(Nk):
                    if s == 0:
                        a1[0] = a1[0] + zm * z[q][m]
                        a2[0] = a2[0] + zm * z[q][(-m) % Nk]
                        a3[0] = a3[0] + zm * conj(z[q][(-m) % Nk])
                    else:
                        a1[s] = a1[s] + zm * z[q][(m + s) % Nk]
                        a2[s] = a2[s] + zm * z[q][(s - m) % Nk]
                        a3[s] = a3[s] + zm * conj(z[q][(s - m) % Nk])
            r1.append(a1)
            r2.append(a2)
            r3.append(a3)
        V1.append(r1)
        V2.append(r2)
        Q1.append(r3)
    return V1, V2, Q1


def phi_tab(z, V1, V2, K, gRB, CX, M, Nk):
    """Phi(Sigma) in the shifted-table arrangement.

    Phi(n,j) = gRB sum_a K(n,a) <z(a)> + CX sum_{a,c} K(n,a)K(n,c) sum_x z(b,x)
               V2[a][c][(j+x) mod Nk],   b = a+c-n,
    with V2[a][c][s] = sum_m z(a,m) z(c,(s-m) mod Nk) and <z(a)> = (1/Nk) sum_j z(a,j).
    (The triple momentum convolution of the skeleton is written as a cyclic
    correlation of z(b,.) against the pair table V2, i.e. nothing pre-summed over
    three indices at once.)
    """
    zero = _zero_like(z)
    means = []
    for a in range(M):
        acc = zero
        for j in range(Nk):
            acc = acc + z[a][j]
        means.append(acc / Nk)
    out = [[zero] * Nk for _ in range(M)]
    for n in range(M):
        Kn = K[n]
        for j in range(Nk):
            s1 = zero
            for a in range(M):
                s1 = s1 + Kn[a] * means[a]
            s2 = zero
            for a in range(M):
                for c in range(M):
                    b = a + c - n
                    if 0 <= b < M:
                        w = Kn[a] * Kn[c]
                        acc = zero
                        for x in range(Nk):
                            acc = acc + z[b][x] * V2[a][c][(j + x) % Nk]
                        s2 = s2 + w * acc
            out[n][j] = gRB * s1 + CX * s2
    return out


def jac_tab(z, V1, V2, K, gRB, CX, M, Nk):
    """dPhi/dSigma, per entry, from the two correlation tables.

    dPhi(n,j)/dSigma(a0,j0) = z(a0,j0)^2 [ gRB K(n,a0)/Nk
        + 2 CX K(n,a0) sum_p K(n,p) V1[p+a0-n][p][(j-j0) mod Nk]
        +   CX sum_a K(n,a) K(n,a0-a+n) V2[a][a0-a+n][(j+j0) mod Nk] ]

    (rainbow pair term, the two derivative slots of the summed pair giving the
    symmetric factor 2, and the derivative of the third -- "b" -- momentum pair.)
    """
    zero = _zero_like(z)
    z2 = [[z[a][j] * z[a][j] for j in range(Nk)] for a in range(M)]
    J = [[zero] * (M * Nk) for _ in range(M * Nk)]
    for n in range(M):
        Kn = K[n]
        for a0 in range(M):
            Kna0 = Kn[a0]
            for j0 in range(Nk):
                g = z2[a0][j0]
                col = a0 * Nk + j0
                for j in range(Nk):
                    t_rb = gRB * Kna0 / Nk
                    t_pair = zero
                    for p in range(M):
                        b = p + a0 - n
                        if 0 <= b < M:
                            t_pair = t_pair + Kn[p] * V1[b][p][(j - j0) % Nk]
                    t_third = zero
                    for a in range(M):
                        c = a0 - a + n
                        if 0 <= c < M:
                            t_third = t_third + Kn[a] * Kn[c] * V2[a][c][(j0 + j) % Nk]
                    J[n * Nk + j][col] = g * (t_rb + CX * (2 * Kna0 * t_pair + t_third))
    return J


def gap_tab(z, Q1, V1, K, gRB, CX, M, Nk, drop_t3=False):
    """The model's linearised gap operator (vertex corrected), real symmetric part.

    B(n,j;a0,j0) = |z(n,j)| |z(a0,j0)| [ gRB K(n,a0)/Nk
                     + 2 CX K(n,a0) sum_p K(n,p) V1[p+a0-n][p][(j-j0) mod Nk]
                     +   CX sum_a K(n,a) K(n,a0-a+n) Q1[a][a0-a+n][(j0+j) mod Nk] ]
    with Q1[a][c][s] = sum_m z(a,m) conj(z(c,(s-m) mod Nk)).

    Pair term, from the definition: the momentum dependence is
        sum_p K(n,p) sum_m z(p,m) z(p+a0-n, m + j0 - j)
    and substituting m -> m' = m + j0 - j gives
        sum_{m'} z(p+a0-n,m') z(p, m' + j - j0) = V1[p+a0-n][p][(j-j0) mod Nk],
    which is the shifted-table form used below.  (The v1 transcription used the
    shift (j0-j) with the same index order, which is only correct for z even in the
    momentum index; `--selftest-gap` checks the form against a direct evaluation of
    the definition for a *non-even* Sigma.)
    """
    az = [[cmod(z[a][j]) for j in range(Nk)] for a in range(M)]
    zero = _zero_like(z)
    D = M * Nk
    out = [[zero] * D for _ in range(D)]
    for n in range(M):
        Kn = K[n]
        for j in range(Nk):
            row = out[n * Nk + j]
            for a0 in range(M):
                Kna0 = Kn[a0]
                for j0 in range(Nk):
                    t_pair = zero
                    for p in range(M):
                        b = p + a0 - n
                        if 0 <= b < M:
                            t_pair = t_pair + Kn[p] * V1[b][p][(j - j0) % Nk]
                    t_third = zero
                    if not drop_t3:
                        for a in range(M):
                            c = a0 - a + n
                            if 0 <= c < M:
                                t_third = t_third + Kn[a] * Kn[c] * Q1[a][c][(j0 + j) % Nk]
                    val = Kna0 * (gRB / Nk) + CX * (2 * Kna0 * t_pair + t_third)
                    if isinstance(val, (CI, RI)):
                        val = promote(val).re
                    row[a0 * Nk + j0] = az[n][j] * (az[a0][j0] * val)
    return out, az


# ==================================================================================
#  float stage: constants, propagator, map, branch-selected fixed point
# ==================================================================================
def float_setup(Nk, Nw, lam, m, eta, mu0):
    M = 2 * Nw
    xi = [2.0 * abs(2 * j + 1 - Nk) / Nk - 1.0 - mu0 for j in range(Nk)]
    w = [float(2 * n + 1) for n in range(-Nw, Nw)]
    tau = eta * float(C_INF) * np.sqrt(lam) * m
    nu = np.pi * tau
    K = [[1.0 if n == a else m * m / (m * m + (2 * np.pi * tau * (n - a)) ** 2) for a in range(M)]
         for n in range(M)]
    gRB = 2.0 * lam * tau
    CX = (2.0 * lam * tau / Nk) ** 2
    return dict(M=M, Nk=Nk, Nw=Nw, xi=xi, w=w, tau=tau, nu=nu, K=K, gRB=gRB, CX=CX,
                m=m, lam=lam, eta=eta, mu0=mu0)


def float_z(s, P):
    M, Nk = P["M"], P["Nk"]
    return [[1.0 / (complex(0.0, P["w"][n] * P["nu"]) - P["xi"][j] - s[n][j]) for j in range(Nk)]
            for n in range(M)]


def float_phi(s, P):
    V1, V2, _ = corr_tables(float_z(s, P), P["M"], P["Nk"])
    return phi_tab(float_z(s, P), V1, V2, P["K"], P["gRB"], P["CX"], P["M"], P["Nk"])


def float_me_profile(P):
    """rainbow-only profile: the branch connected to the non-interacting state."""
    M, Nk = P["M"], P["Nk"]
    s = [[0j] * Nk for _ in range(M)]
    for _ in range(4000):
        z = float_z(s, P)
        K, gRB = P["K"], P["gRB"]
        new = [[0j] * Nk for _ in range(M)]
        mv = 0.0
        for n in range(M):
            acc = 0j
            for a in range(M):
                za = sum(z[a][j] for j in range(Nk)) / Nk
                acc += K[n][a] * za
            acc *= gRB
            for j in range(Nk):
                new[n][j] = acc
                mv = max(mv, abs(acc - s[n][j]))
        s = [[0.5 * new[n][j] + 0.5 * s[n][j] for j in range(Nk)] for n in range(M)]
        if mv < 1e-16:
            break
    return s


def float_solve(P, tol=1e-15, itmax=4000):
    """damped fixed-point iteration from the rainbow profile, then Newton steps."""
    M, Nk = P["M"], P["Nk"]
    s = float_me_profile(P)
    for _ in range(itmax):
        sn = float_phi(s, P)
        diff = max(abs(sn[a][j] - s[a][j]) for a in range(M) for j in range(Nk))
        s = [[0.5 * sn[a][j] + 0.5 * s[a][j] for j in range(Nk)] for a in range(M)]
        if diff < tol:
            break
    D = M * Nk
    for _ in range(4):
        V1, V2, _ = corr_tables(float_z(s, P), M, Nk)
        Jf = jac_tab(float_z(s, P), V1, V2, P["K"], P["gRB"], P["CX"], M, Nk)
        Jm = np.array([[complex(x) for x in r] for r in Jf])
        ph = float_phi(s, P)
        r = np.array([ph[a][j] - s[a][j] for a in range(M) for j in range(Nk)])
        A = np.linalg.inv(np.eye(D) - Jm)
        dx = A @ r
        for a in range(M):
            for j in range(Nk):
                s[a][j] += dx[a * Nk + j]
        if np.max(np.abs(dx)) < 1e-18:
            break
    return s


# ==================================================================================
#  interval stage: constants and the Sigma box
# ==================================================================================
def box_around(s0, rho_str, M, Nk):
    """Sup-norm box  {|Re(Sigma - Sigma_0)|, |Im(Sigma - Sigma_0)| <= rho}  as intervals.

    The endpoints are formed by directed-rounding addition/subtraction of the exact
    double Sigma_0, so [v - r, v + r] is guaranteed to contain the exact box.
    """
    r = iv_mpf(rho_str)
    box = []
    for a in range(M):
        row = []
        for j in range(Nk):
            vr = RI.exact(float(s0[a][j].real))
            vi = RI.exact(float(s0[a][j].imag))
            re = RI((vr - r).a, (vr + r).b)
            im = RI((vi - r).a, (vi + r).b)
            row.append(CI(re, im))
        box.append(row)
    return box, r


def krawczyk_modulus(Ainv, DF, D):
    """Rigorous upper bound of || I - A (I - DPhi) ||_inf  (A = exact double matrix).

    The diagonal-dominance alternative |||A||I-DPhi||| is always >= 1 because
    A (I - DPhi) = I, so the signed matrix product has to be formed.  Done in
    python-flint ball arithmetic (`acb_mat`): every input interval is embedded by an
    exact lift of its mpfr endpoints, inflated by a relative 2^-60 so that the ball
    radius cannot be under-rounded by the constructor, and every subsequent
    operation is outward rounded by construction.  Row sums of the resulting ball
    matrix are taken with `arb.max` (a hull, not an undecided comparison), so the
    overlap of two balls can never drop an upper bound.
    """
    import flint
    from flint import acb, acb_mat, arb
    saved = flint.ctx.prec
    flint.ctx.prec = 220
    pad = arb(1) + arb(2) ** -60

    def lift(x):
        q = gmpy2.mpq(x)
        return arb(int(q.numerator)) / arb(int(q.denominator))

    def ball(lo, hi):
        """arb ball containing the real interval [lo, hi]."""
        a, b = lift(lo), lift(hi)
        rad = (b - a) * pad if b != a else arb(0)
        return arb(a, rad)

    try:
        one = acb(arb(1), arb(0))
        E = acb_mat(D, D)                      # E = I - DPhi
        for i in range(D):
            for j in range(D):
                e = DF[i][j]
                c = acb(ball(e.re.a, e.re.b), ball(e.im.a, e.im.b))
                E[i, j] = (one - c) if i == j else (-c)
        A = acb_mat(D, D)
        for i in range(D):
            row = Ainv[i]
            for j in range(D):
                z = row[j]
                A[i, j] = acb(arb(float(z.real)), arb(float(z.imag)))
        P = A * E
        worst = arb(0)
        for i in range(D):
            acc = arb(0)
            for j in range(D):
                entry = (arb(1) - P[i, j]) if i == j else (-P[i, j])
                acc = acc + entry.abs_upper()
            worst = worst.max(acc)          # hull: overlapping balls keep the larger
        # worst.upper() is a rigorous upper bound (arf); float() rounds to nearest,
        # so inflate by a relative 2^-40 to keep the float lift one-sided.
        with gmpy2.context(round=gmpy2.RoundUp):
            return gmpy2.mpfr(float(worst.upper())) * (1 + gmpy2.mpfr(2) ** -40)
    finally:
        flint.ctx.prec = saved


def iv_setup(Nk, Nw, lam, m_num, m_den, eta_num, eta_den, mu_num, mu_den, s0, rho_str):
    M = 2 * Nw
    m = iv_mpf(m_num) / iv_mpf(m_den)                       # m = Omega/W, exact rational
    mu0 = iv_mpf(mu_num) / iv_mpf(mu_den)
    et = iv_mpf(eta_num) / iv_mpf(eta_den)
    cinf = iv_mpf(C_INF)
    lamv = iv_mpf(lam)
    tau = et * cinf * lamv.sqrt() * m
    nu = iv_pi() * tau
    w = [iv_mpf(2 * n + 1) for n in range(-Nw, Nw)]
    xi = [iv_mpf(2 * abs(2 * j + 1 - Nk)) / iv_mpf(Nk) - iv_mpf(1) - mu0 for j in range(Nk)]
    K = []
    for n in range(M):
        row = []
        for a in range(M):
            if n == a:
                row.append(iv_mpf(1))
            else:
                dd = iv_mpf((n - a) ** 2)
                row.append(m * m / (m * m + (iv_mpf(4) * iv_pi() ** 2 * tau * tau) * dd))
        K.append(row)
    gRB = iv_mpf(2) * lamv * tau
    CX = (iv_mpf(2) * lamv * tau / iv_mpf(Nk)) ** 2
    sball, r = box_around(s0, rho_str, M, Nk)
    return dict(M=M, Nk=Nk, w=w, nu=nu, xi=xi, K=K, gRB=gRB, CX=CX, sball=sball, rho=r,
                tau=tau)


def iv_z(Q):
    """z over the box, plus a rigorous lower bound on min |Im denominator|."""
    M, Nk = Q["M"], Q["Nk"]
    z = []
    min_im = None
    for n in range(M):
        row = []
        for j in range(Nk):
            re = -(Q["xi"][j]) - Q["sball"][n][j].re
            im = Q["w"][n] * Q["nu"] - Q["sball"][n][j].im
            row.append(CI(iv_mpf(1)) / CI(re, im))
            # |Im| >= min(|lo|,|hi|) when the interval does not contain 0
            if im.a > 0:
                low = im.a
            elif im.b < 0:
                low = _neg_exact(im.b)
            else:
                low = gmpy2.mpfr(0)
            min_im = low if min_im is None else min(min_im, low)
        z.append(row)
    return z, min_im


# ==================================================================================
#  the check
# ==================================================================================
def load_reference():
    import importlib.util
    import os
    from pathlib import Path
    env = os.environ.get("VERTEX_UNIV")
    root = Path(env) if env else Path(__file__).resolve().parent.parent
    path = root / "verify_1d_counterexample.py"
    if not path.exists():
        raise FileNotFoundError(str(path))
    spec = importlib.util.spec_from_file_location("vertex_ref_indep", str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check(Nk=24, Nw=12, lam=8, m_num=1, m_den=100, eta_num=1, eta_den=1,
          mu_num=-2, mu_den=5, rho=RHO, theta=THETA, verbose=True, skip_ref=False,
          mutate="none"):
    t00 = time.time()
    if mutate not in ("none",) + MUTANTS:
        raise SystemExit("unknown mutation %r" % (mutate,))
    if mutate == "m_threshold":
        m_num, m_den = 2, 100                      # (below the crossing)
    if mutate == "rho_zero":
        rho = "0"
    theta_use = "1e-2" if mutate == "ray_thresh" else theta
    drop_t3 = (mutate == "drop_lb_t3")
    say = (lambda s: print(s, flush=True)) if verbose else (lambda s: None)
    m = m_num / m_den
    mu0 = mu_num / mu_den
    eta = eta_num / eta_den
    P = float_setup(Nk, Nw, lam, m, eta, mu0)
    if mutate == "no_vertex":
        P["CX"] = 0.0
    M, D = P["M"], P["M"] * Nk
    say("independent checker (MPFR directed rounding, %d bits): Nk=%d Nw=%d M=%d D=%d "
        "lam=%d m=%d/%d eta=%d/%d mu0=%d/%d rho=%s theta=%s"
        % (PREC, Nk, Nw, M, D, lam, m_num, m_den, eta_num, eta_den, mu_num, mu_den, rho, theta))
    fails = []

    # ---------------- float candidate and transcription self-tests
    s0 = float_solve(P)
    ph0 = float_phi(s0, P)
    res_f = max(abs(ph0[a][j] - s0[a][j]) for a in range(M) for j in range(Nk))
    say("[A] float fixed point candidate: ||Phi(S0)-S0||_inf = %.3e" % res_f)
    if res_f > 1e-12:
        fails.append("float candidate is not a fixed point (%.1e)" % res_f)

    V1f, V2f, _ = corr_tables(float_z(s0, P), M, Nk)
    Jf = jac_tab(float_z(s0, P), V1f, V2f, P["K"], P["gRB"], P["CX"], M, Nk)
    Jm = np.array([[complex(x) for x in r] for r in Jf])
    rng = np.random.default_rng(2024)
    err_fd = 0.0
    for _ in range(6):
        a0 = int(rng.integers(M))
        j0 = int(rng.integers(Nk))
        h = 1e-7
        sp = [[s0[a][j] for j in range(Nk)] for a in range(M)]
        sm = [[s0[a][j] for j in range(Nk)] for a in range(M)]
        sp[a0][j0] += h
        sm[a0][j0] -= h
        fp, fm = float_phi(sp, P), float_phi(sm, P)
        col = np.array([(fp[a][j] - fm[a][j]) / (2 * h) for a in range(M) for j in range(Nk)])
        err_fd = max(err_fd, float(np.max(np.abs(col - Jm[:, a0 * Nk + j0]))))
    V1g, _, Q1f = corr_tables(float_z(s0, P), M, Nk)
    Bf, _ = gap_tab(float_z(s0, P), Q1f, V1g, P["K"], P["gRB"], P["CX"], M, Nk)
    Bm_f = np.array([[complex(x) for x in r] for r in Bf])
    lb_float = float(np.linalg.eigvalsh(0.5 * (Bm_f + Bm_f.conj().T))[-1])
    say("[B] Jacobian vs central differences: %.2e ; float lambda_max = %.9f (asym %.1e)"
        % (err_fd, lb_float, float(np.max(np.abs(Bm_f - Bm_f.conj().T)))))
    if err_fd > 1e-5:
        fails.append("Jacobian differs from central differences of Phi (%.1e)" % err_fd)

    if not skip_ref:
        try:
            ref = load_reference()
            out = ref.solve_finite_m_lattice(lam=float(lam), m_omega=m, eta=eta, mu0=mu0,
                                             N_k=Nk, N_w=Nw)
            d = abs(out["full"] - lb_float)
            say("[C] audited artifact solve_finite_m_lattice: lambda_full %.9f (delta %.1e) ; "
                "lambda_me %.9f ; min_u %.4f" % (out["full"], d, out["me"], out["min_u"]))
            if d > 1e-5:
                fails.append("model disagrees with the audited artifact (%.1e)" % d)
        except Exception as exc:                       # pragma: no cover
            say("[C] audited artifact cross-check skipped: %s" % exc)

    # ---------------- interval enclosures over the box
    ti = time.time()
    Q = iv_setup(Nk, Nw, lam, m_num, m_den, eta_num, eta_den, mu_num, mu_den, s0, rho)
    if mutate == "no_vertex":
        Q["CX"] = iv_mpf(0)
    if mutate == "rb_scale":
        Q["gRB"] = Q["gRB"] * RI.rational(6, 5)

    # ---------------- guards: the box really is a box, and 1/Nk is exact
    if mutate == "point_not_ball":
        Q["sball"], _ = box_around(s0, "0", M, Nk)      # the v1 fault: a point, not a box
    rho_g = iv_mpf(rho)
    box_ok = True
    for a in range(M):
        for j in range(Nk):
            vr = RI.exact(float(s0[a][j].real))
            vi = RI.exact(float(s0[a][j].imag))
            lo_re, hi_re = (vr - rho_g).a, (vr + rho_g).b
            lo_im, hi_im = (vi - rho_g).a, (vi + rho_g).b
            if not (Q["sball"][a][j].re.a <= lo_re and Q["sball"][a][j].re.b >= hi_re
                    and Q["sball"][a][j].im.a <= lo_im and Q["sball"][a][j].im.b >= hi_im):
                box_ok = False
    invNk_g = RI.exact(1.0 / Nk) if mutate == "float_leak" else iv_mpf(1) / iv_mpf(Nk)
    inv_ok = bool((invNk_g * iv_mpf(Nk)).contains(1))
    say("[D0] guards: box covers Sigma_0 +- rho in every component -> %s ; "
        "1/Nk exact ((1/Nk)*Nk contains 1) -> %s" % (box_ok, inv_ok))
    if not box_ok:
        fails.append("box guard: the enclosure does not cover the Sigma box of radius rho")
    if not inv_ok:
        fails.append("constant guard: 1/Nk is not an exact reciprocal (binary64 leak)")
    z, min_im = iv_z(Q)
    say("[D] min |Im (i w_n nu - xi - Sigma)| over the box >= %.6e" % _f(min_im))
    if not (min_im > 0):
        fails.append("box reaches a propagator pole")
    V1, V2, Q1 = corr_tables(z, M, Nk)
    PHI = phi_tab(z, V1, V2, Q["K"], Q["gRB"], Q["CX"], M, Nk)
    say("[E] tables and Phi(box) enclosed (%.1fs)" % (time.time() - ti))
    # Krawczyk residual at the *point* centre Sigma_0 (widths are then pure rounding,
    # ~2^-160 per operation); the box only enters through DPhi.
    Q0 = iv_setup(Nk, Nw, lam, m_num, m_den, eta_num, eta_den, mu_num, mu_den, s0, "0")
    if mutate == "no_vertex":
        Q0["CX"] = iv_mpf(0)
    if mutate == "rb_scale":
        Q0["gRB"] = Q0["gRB"] * RI.rational(6, 5)
    z0, _ = iv_z(Q0)
    W1, W2, _ = corr_tables(z0, M, Nk)
    PHI0 = phi_tab(z0, W1, W2, Q0["K"], Q0["gRB"], Q0["CX"], M, Nk)
    res = [[PHI0[a][j] - CI(iv_mpf(float(s0[a][j].real)), iv_mpf(float(s0[a][j].imag)))
            for j in range(Nk)] for a in range(M)]
    if mutate == "res_scale":
        big = iv_mpf(2) ** 20
        res = [[r * big for r in row] for row in res]
    DF = jac_tab(z, V1, V2, Q["K"], Q["gRB"], Q["CX"], M, Nk)
    say("[F] DPhi(box) enclosed (%.1fs)" % (time.time() - ti))

    # ---------------- Krawczyk / Banach:  A = inv(I - DF(Sigma_0))  (float doubles)
    Ainv = np.linalg.inv(np.eye(D) - Jm)
    rowA = float(np.max(np.abs(Ainv).sum(axis=1)))
    say("[G] ||A||_inf = %.3e (A's nonsingularity follows from k < 1, not from this float)" % rowA)
    r_hi = gmpy2.mpfr(0)
    for i in range(D):
        acc = promote(0)
        for l in range(D):
            a_il = Ainv[i, l]
            if a_il != 0.0:
                acc = acc + CI(iv_mpf(float(a_il.real)), iv_mpf(float(a_il.imag))) * res[l // Nk][l % Nk]
        b = cmod(acc).b
        if b > r_hi:
            r_hi = b
    # k = || I - A (I - DPhi) ||_inf over the box, in arb (matrix product, signed
    # cancellation preserved: a triangle-inequality bound of |A||I-DPhi| would be >= 1
    # and useless).  A is the exact double matrix from the float stage.
    k_hi = krawczyk_modulus(Ainv, DF, D)
    rho_iv = iv_mpf(rho)
    r_iv = RI(gmpy2.mpfr(0), r_hi)
    k_iv = RI(gmpy2.mpfr(0), k_hi)
    cond_ctr = bool(k_hi < 1)
    cond_map = bool((r_iv + k_iv * rho_iv).b <= rho_iv.a)
    say("[H] Banach/Krawczyk: r = ||A(Phi(S0)-S0)||_inf <= %s ; "
        "k = ||I - A(I-DPhi)||_inf <= %s" % (_dup(r_hi), _dup(k_hi)))
    say("    r + k rho <= %s <= rho = %s -> %s ; k < 1 -> %s ; unique fixed point with "
        "||Sigma* - Sigma_0||_inf <= r/(1-k) <= %s"
        % (_dup((r_iv + k_iv * rho_iv).b), _dlo(rho_iv.a), cond_map, cond_ctr,
           _dup((r_iv / (RI.exact(1) - k_iv)).b) if cond_ctr else "inf"))
    if not cond_ctr:
        fails.append("contraction modulus k >= 1")
    if not cond_map:
        fails.append("Krawczyk/Banach containment r + k rho <= rho fails")

    # ---------------- Rayleigh quotient over the box
    tR = time.time()
    Biv, az = gap_tab(z, Q1, V1, Q["K"], Q["gRB"], Q["CX"], M, Nk, drop_t3=drop_t3)
    if mutate == "diag_drop":
        for i in range(D):
            Biv[i][i] = Biv[i][i] - RI.rational(1, 50)
    w, Vv = np.linalg.eigh(0.5 * (Bm_f + Bm_f.conj().T))
    v = Vv[:, -1].real
    num = promote(0)
    den = RI.exact(0)
    vr = [RI.exact(float(v[i])) for i in range(D)]
    for i in range(D):
        den = den + vr[i] * vr[i]
        acc = promote(0)
        for j in range(D):
            s_ij = (Biv[i][j] + Biv[j][i]) * iv_mpf(1) / iv_mpf(2)
            acc = acc + promote(s_ij) * vr[j]
        num = num + promote(vr[i]) * acc
    q = promote(num) / promote(den)
    assert den.a > 0
    q_lo = q.re.a                        # lower endpoint: <v,Sv>/<v,v> >= q_lo (all signs)
    delta = gmpy2.mpfr(0)                # diagnostic only: ||B - B^T||_inf upper bound
    for i in range(D):
        acc = RI.exact(0)
        for j in range(D):
            if i != j:
                acc = acc + cmod(Biv[i][j] - Biv[j][i])
        if acc.b > delta:
            delta = acc.b
    # Rayleigh lower bound on the symmetric part: no antisymmetry penalty is
    # subtracted (S = (B+B^T)/2 is itself the object enclosed here, and for a real
    # symmetric matrix lambda_max(S) >= <v,Sv>/<v,v> for every real v != 0).
    bound_lo = q_lo
    thr = (iv_mpf(1) + iv_mpf(theta_use)).b
    ok_ray = bool(bound_lo >= thr)
    say("[I] Rayleigh: <v,B v>/<v,v> >= %s  (B enclosed entrywise, symmetrised) ; "
        "lambda_max >= %s  (%.1fs)" % (_dlo(q_lo), _dlo(bound_lo), time.time() - tR))
    say("    diagnostic only: ||B-B^T||_inf <= %s (enclosure width);  "
        "criterion: bound >= 1 + %s  (slack >= %s) -> %s"
        % (_dup(delta, 4), theta_use, _dlo(bound_lo - thr), ok_ray))
    if not ok_ray:
        fails.append("Rayleigh lower bound below 1 + theta")

    say("[J] %s  (%s)   wall %.1fs"
        % ("PASS" if not fails else "FAIL",
           "; ".join(fails) if fails else "all certified inequalities hold",
           time.time() - t00))
    return fails


# ==================================================================================
#  transcription self-test of the shifted-table forms against the definition
#  (run on a small lattice at a *non-even* Sigma, where a wrong momentum shift shows)
# ==================================================================================
def _mid_float(x):
    """midpoint of an RI / CI (or a plain number) as float / complex (test output only)."""
    if isinstance(x, CI):
        return complex(_mid_float(x.re), _mid_float(x.im))
    if isinstance(x, RI):
        return float((x.a + x.b) / 2)
    return x


def sx_ref(F, G, H, K, CX, M, Nk):
    """Crossed-vertex sum, transcribed directly from the audited artifact's
    `cx_op_nonlinear`/`apply_full` structure (FFT in the Matsubara index, mirrored
    momentum index on the H slot).  This is the *definition* used by the
    transcription test; the tables under test are V1/V2/Q1."""
    idxneg = (-np.arange(Nk)) % Nk
    Ff = np.fft.fft(F, axis=1)
    Gf = np.fft.fft(G, axis=1)
    Hf = np.fft.fft(H[:, idxneg], axis=1)
    b = np.arange(M)[:, None, None] + np.arange(M)[None, :, None] - np.arange(M)[None, None, :]
    out = np.zeros((M, Nk), dtype=complex)
    for n in range(M):
        acc = np.zeros(Nk, dtype=complex)
        for a in range(M):
            bs = b[a, :, n]
            ok = (bs >= 0) & (bs < M)
            if not ok.any():
                continue
            cs = np.nonzero(ok)[0]
            term = np.array([K[n][c] for c in cs])[:, None] * (Gf[cs] * Hf[bs[cs]])
            acc += K[n][a] * Ff[a] * term.sum(axis=0)
        out[n] = acc
    return CX * np.fft.ifft(out, axis=1)


def rb_ref(F, K, gRB, M, Nk):
    km = np.asarray(K, dtype=float)
    return gRB * (km @ np.mean(F, axis=1))[:, None] * np.ones(Nk)


def phi_ref(z, K, gRB, CX, M, Nk):
    zz = np.asarray(z, dtype=complex)
    return rb_ref(zz, K, gRB, M, Nk) + sx_ref(zz, zz, zz, K, CX, M, Nk)


def apply_ref(z, v, K, gRB, CX, M, Nk):
    """The model's linearised gap operator (definitional form, artifact structure)."""
    zz = np.asarray(z, dtype=complex)
    vv = np.asarray(v, dtype=complex).reshape(M, Nk)
    az = np.abs(zz)
    h = az * vv
    return az * (rb_ref(h, K, gRB, M, Nk)
                 + 2.0 * sx_ref(h, zz, zz, K, CX, M, Nk).real
                 + sx_ref(zz, np.conj(zz), h, K, CX, M, Nk).real)


def selftest_gap(Nk=5, Nw=3, seed=3):
    """Compare phi_tab / gap_tab with the direct definitions at a NON-even Sigma.

    The v1 transcription error in `gap_tab` (the pair term used the shift (j0-j) with
    the index order V1[b][p]) is invisible for a Sigma that is even in the momentum
    index, so this test uses a Sigma with an explicit odd component in j.
    """
    M = 2 * Nw
    rng = np.random.default_rng(seed)
    K = [[1.0 if n == a else 1.0 / (1.0 + (n - a) ** 2) for a in range(M)] for n in range(M)]
    gRB, CX = 0.05, 1e-3
    xi = [2.0 * abs(2 * j + 1 - Nk) / Nk - 1.0 + 0.4 for j in range(Nk)]
    w = [float(2 * n + 1) for n in range(-Nw, Nw)]
    s = [[(0.011 * rng.normal()) + 0.02 * j + 1j * (0.010 * rng.normal())
          for j in range(Nk)] for n in range(M)]
    odd = max(abs(s[n][j] - s[n][(-j) % Nk]) for n in range(M) for j in range(Nk))
    z = [[1.0 / (1j * w[n] * 0.7 - xi[j] - s[n][j]) for j in range(Nk)] for n in range(M)]
    V1, V2, Q1 = corr_tables(z, M, Nk)
    ph_t = np.asarray(phi_tab(z, V1, V2, K, gRB, CX, M, Nk), dtype=complex)
    ph_d = phi_ref(z, K, gRB, CX, M, Nk)
    e_phi = float(np.max(np.abs(ph_t - ph_d)))
    B_t, _ = gap_tab(z, Q1, V1, K, gRB, CX, M, Nk)
    B_t = np.array([[_mid_float(x) for x in row] for row in B_t], dtype=complex)
    D = M * Nk
    B_r = np.zeros((D, D), dtype=complex)
    for k in range(D):
        e = np.zeros(D, dtype=complex)
        e[k] = 1.0
        B_r[:, k] = apply_ref(z, e, K, gRB, CX, M, Nk).ravel()
    e_B = float(np.max(np.abs(np.asarray(B_t) - B_r)))
    sc = float(np.max(np.abs(B_r)))
    print("selftest-gap (Nk=%d, M=%d, non-even Sigma: max|s(n,j)-s(n,-j)| = %.3e)"
          % (Nk, M, odd))
    print("  phi_tab vs definition:  max abs err %.3e" % e_phi)
    print("  gap_tab vs definition:  max abs err %.3e  (scale %.3e, rel %.3e)"
          % (e_B, sc, e_B / sc))
    ok = (e_phi < 1e-12) and (e_B < 1e-12) and odd > 1e-3
    print("  -> %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


# ==================================================================================
#  backend self-test against exact rational arithmetic
# ==================================================================================
def selftest(ntrial=4000, seed=7):
    """Validate RI/CI against exact rational arithmetic (`fractions.Fraction`).

    Random rationals p, q are lifted to intervals; every enclosure produced by an
    operation must contain the exact rational result, and the width must be small.
    """
    rng = np.random.default_rng(seed)
    t0 = time.time()
    bad = 0
    width_bad = 0

    def qq(fr):
        return gmpy2.mpq(fr.numerator) / gmpy2.mpq(fr.denominator)

    def inside(iv, fr):
        """exact-rational membership test of an endpoint interval (no rounding)."""
        e = qq(fr)
        return bool(gmpy2.mpq(iv.a) <= e and e <= gmpy2.mpq(iv.b))

    def down(fr):
        with gmpy2.context(precision=PREC, round=RoundDown):
            return gmpy2.mpfr(qq(fr))

    def up(fr):
        with gmpy2.context(precision=PREC, round=RoundUp):
            return gmpy2.mpfr(qq(fr))

    def sqrt_exact_lo(fr):
        with gmpy2.context(precision=PREC, round=RoundDown):
            return gmpy2.sqrt(gmpy2.mpfr(qq(fr)))

    for _ in range(ntrial):
        pn, pd = int(rng.integers(-2 ** 40, 2 ** 40)), int(rng.integers(1, 2 ** 40))
        qn, qd = int(rng.integers(-2 ** 40, 2 ** 40)), int(rng.integers(1, 2 ** 40))
        rn, rd = int(rng.integers(-2 ** 40, 2 ** 40)), int(rng.integers(1, 2 ** 40))
        rlist = int(rng.integers(0, 2 ** 40))

        A = RI.rational(pn, pd)
        B = RI.rational(qn, qd)
        C = RI.rational(rn, rd)

        # absolute-width criterion: the enclosure of one operation must be tighter
        # than 2^(-PREC+8) relative to the operand scale (cancellation is allowed, so
        # the test is absolute, not relative to the result).
        scale = gmpy2.mpq(max(abs(pn), abs(pd), abs(qn), abs(qd), abs(rn), abs(rd)))
        lim = gmpy2.mpq(1) / gmpy2.mpq(2) ** (PREC - 8) * (scale + 1)
        for iv_acc, ex_acc in ((A * B, Fraction(pn * qn, pd * qd)),
                               (A + B, Fraction(pn, pd) + Fraction(qn, qd)),
                               (A - B, Fraction(pn, pd) - Fraction(qn, qd)),
                               ((A * B) * C, Fraction(pn * qn * rn, pd * qd * rd)),
                               ((A + B) * C - A, (Fraction(pn, pd) + Fraction(qn, qd))
                                * Fraction(rn, rd) - Fraction(pn, pd))):
            if not inside(iv_acc, ex_acc):
                bad += 1
            if not (gmpy2.mpq(iv_acc.b) - gmpy2.mpq(iv_acc.a)) <= lim:
                width_bad += 1

        if qn != 0:
            exd = Fraction(pn * qd, pd * qn)
            D = A / B
            if not inside(D, exd):
                bad += 1
            if not (gmpy2.mpq(D.b) - gmpy2.mpq(D.a)) <= lim / abs(qq(Fraction(qn, qd))):
                width_bad += 1

        if pn > 0:
            S = RI.rational(pn, pd).sqrt()
            sq_lo = sqrt_exact_lo(Fraction(pn, pd))
            if not (S.a <= sq_lo and gmpy2.mpq(S.a) ** 2 <= qq(Fraction(pn, pd))
                    <= gmpy2.mpq(S.b) ** 2):
                bad += 1

        # complex algebra: (a + b i) * (c + d i) and division
        ca = CI(RI.rational(pn, pd), RI.rational(qn, qd))
        cb = CI(RI.rational(rlist, 1), RI.rational(rd, qd))
        cp = ca * cb
        a_, b_, c_, d_ = (Fraction(pn, pd), Fraction(qn, qd), Fraction(rlist, 1),
                          Fraction(rd, qd))
        if not inside(cp.re, a_ * c_ - b_ * d_) or not inside(cp.im, a_ * d_ + b_ * c_):
            bad += 1
        if abs(a_) + abs(b_) > 0:
            cq = (ca * ca.conjugate()).re
            if not inside(cq, a_ * a_ + b_ * b_):
                bad += 1
            if not (cq.a >= 0):
                bad += 1
        if abs(a_) > Fraction(1, 10 ** 6):
            mod = cmod(ca)
            ex2 = qq(a_ * a_ + b_ * b_)
            if not (mod.a <= sqrt_exact_lo(a_ * a_ + b_ * b_)
                    and gmpy2.mpq(mod.a) ** 2 <= ex2 <= gmpy2.mpq(mod.b) ** 2):
                bad += 1

        # exact lifts must be exact
        E = RI.exact(float(Fraction(pn, pd)))
        if E.b - E.a != 0:
            bad += 1
        dec_fr = Fraction(18272624777232952, 10 ** 17)
        dec = RI.decimal("0.18272624777232952")
        if not (dec.a == down(dec_fr) and dec.b == up(dec_fr) and inside(dec, dec_fr)):
            bad += 1

    print("selftest: %d trials, containment failures %d, width failures %d (%.1fs)"
          % (ntrial, bad, width_bad, time.time() - t0))
    return 0 if (bad == 0 and width_bad == 0) else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="independent second checker (MPFR directed rounding)")
    ap.add_argument("--Nk", type=int, default=24)
    ap.add_argument("--Nw", type=int, default=12)
    ap.add_argument("--lam", type=int, default=8)
    ap.add_argument("--m", default="1/100")
    ap.add_argument("--eta", default="1/1")
    ap.add_argument("--mu0", default="-2/5")
    ap.add_argument("--rho", default=RHO)
    ap.add_argument("--theta", default=THETA)
    ap.add_argument("--no-ref", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--selftest-gap", action="store_true")
    ap.add_argument("--mutate", default="none", choices=("none",) + MUTANTS)
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if a.selftest_gap:
        return selftest_gap()
    mn, md = a.m.split("/")
    en, ed = a.eta.split("/")
    un, ud = a.mu0.split("/")
    fails = check(Nk=a.Nk, Nw=a.Nw, lam=a.lam, m_num=int(mn), m_den=int(md),
                  eta_num=int(en), eta_den=int(ed), mu_num=int(un), mu_den=int(ud),
                  rho=a.rho, theta=a.theta, skip_ref=a.no_ref, mutate=a.mutate)
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
