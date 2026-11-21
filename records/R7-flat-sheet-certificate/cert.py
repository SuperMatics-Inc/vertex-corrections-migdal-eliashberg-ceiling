#!/usr/bin/env python3
"""cert.py -- rigorous interval-arithmetic certificate for the VERTEX-UNIV
1D flat-sheet, vertex-corrected counterexample at the Migdal-Eliashberg ceiling.

MODEL  (finite lattice; exact rational constants; this is the object that is certified)
-------------------------------------------------------------------------------------
M(N_k, N_w; lambda, m, eta):  strictly 1D periodic band with two flat sheets,

  momentum   k_j = -pi + (j + 1/2) 2 pi / N_k,          j = 0 .. N_k-1
  dispersion xi_j = 2 |2j+1-N_k| / N_k - 1 - mu0        ( = 2|k_j|/pi - 1 - mu0 ),  mu0 = -2/5
  Matsubara  w_n = 2n + 1,  n = -N_w .. N_w-1,          M = 2 N_w
  temperature at the Migdal-Eliashberg ceiling:
             tau = eta * C_inf * sqrt(lambda) * m,   nu = pi tau,
             C_inf = 0.18272624777232952,  m = Omega/W = 1/100,  eta = 1
  phonon kernel  K(d) = m^2 / (m^2 + (2 pi tau d)^2)   (Lambda = 2 pi tau, i.e. omega = Omega)
  couplings      gRB = 2 lambda tau,   CX = (2 lambda tau / N_k)^2
  z(n,j)      = 1 / (i w_n nu - xi_j - Sigma(n,j))      (Nambu off-diagonal propagator)

  normal-state map (rainbow + first-crossed skeleton, Schrodi-Oppeneer-Aperis
  arXiv:1911.12872v2 Eq. (4), bare phonon, fixed bare chemical potential, no Coulomb):

  Phi(Sigma)(n,j) = gRB sum_a K(n-a) (1/N_k) sum_m z(a,m)
                  + CX  sum_{a,c} K(n-a) K(n-c) sum_{m,m'} z(a,m) z(c,m') z(b,(m+m'-j) mod N_k),
      b = a + c - n  (terms with b outside 0..M-1 omitted)

  linearised gap operator of the model (az gauge, .real taken, as in the reference
  implementation verify_1d_counterexample.py::solve_finite_m_lattice):

  B(Sigma)_{(n,j),(a_0,j_0)} = az(n,j) Re{ az(a_0,j_0) [
        gRB K(n-a_0) / N_k
      + 2 CX K(n-a_0) sum_{p: 0 <= b' < M} K(n-p) sum_m z(p,m) z(b', m+u) ]
      + CX az(a_0,j_0) sum_{a: 0 <= a_0-a+n < M} K(n-a) K(n-c)
                      sum_m z(a,m) conj z(c,(w-m) mod N_k) },
      az = |z|,  u = (j_0 - j) mod N_k,  b' = p + a_0 - n,  w = (j_0 + j) mod N_k,
      c = a_0 - a + n   (i.e. the third sum runs over a and c with a + c = a_0 + n).

CERTIFIED STATEMENTS (all decided inequalities in python-flint arb/acb ball arithmetic)
---------------------------------------------------------------------------------------
  (N) FIXED POINT.  Phi has exactly one fixed point Sigma* in the sup-norm box
      ||Sigma - Sigma_0||_inf <= rho, proved by the Krawczyk/Banach criterion
      T(Sigma) = Sigma - A(Phi(Sigma)-Sigma) with A = inv(I - DPhi(Sigma_0)):
          || A (Phi(Sigma_0) - Sigma_0) || + || I - A (I - DPhi(box)) || rho <= rho,
          || I - A (I - DPhi(box)) || < 1.
  (K) CROSSING.  For every Sigma in that box, the model's linearised gap operator
      B(Sigma) satisfies
          lambda_max(B(Sigma)) >= 1 + theta          (theta = 5e-3 by default),
      by a rigorous enclosure, over the whole box, of the Rayleigh quotient
      <v, B(Sigma) v> / <v, v> of one fixed real vector v.  B(Sigma) is real
      symmetric for EVERY Sigma (CERT.md section 3, Proposition 1), so no
      antisymmetry penalty is needed or subtracted (the penalty used in cert.v1.py
      was unsound, see REVIEW3).  The enclosed ||B - B^T||_inf is printed only as a
      diagnostic of the enclosure width.

The floating point stage only *selects* Sigma_0 and v.  Exit code 0 iff every
certified inequality of every requested model holds.

Usage (everything below has been executed; see CERT.md for the exit codes):
  uv run --no-project --with mpmath --with numpy --with scipy --with python-flint \
      python cert.py
  ... python cert.py --set ladder --json receipt.json
  ... python cert.py --mutate no_vertex          (must exit 1)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time

import numpy as np

try:
    from flint import acb, acb_mat, arb, ctx
except Exception:  # pragma: no cover
    sys.exit("cert.py requires python-flint:  uv run --no-project --with python-flint ...")

C_INF = "0.18272624777232952"
THETA = 5e-3          # certified margin: lambda_max >= 1 + THETA
RHO = "1e-12"         # sup-norm radius of the fixed-point box
PREC = 200            # arb working precision (bits)

# name -> (lambda, m numerator/denominator, eta, mu0 numerator/denominator, N_k, N_w)
MODELS = {
    "main": (8, (1, 100), 1, (-2, 5), 24, 12),      # certified: D = 576
    "m200": (8, (1, 200), 1, (-2, 5), 24, 12),      # certified: D = 576 (vertex decisive: ME column < 1)
    "m400": (8, (1, 400), 1, (-2, 5), 24, 12),      # NUMERICAL only (rainbow column alone already > 1 here)
    "nw16": (8, (1, 100), 1, (-2, 5), 24, 16),      # NUMERICAL only
    "nk48": (8, (1, 100), 1, (-2, 5), 48, 12),      # NUMERICAL only
    "nk64": (8, (1, 100), 1, (-2, 5), 64, 12),      # NUMERICAL only
    "m20": (8, (2, 100), 1, (-2, 5), 24, 12),       # NUMERICAL only (below the crossing)
    "lam6": (6, (1, 100), 1, (-2, 5), 24, 12),      # NUMERICAL only (crosses as well)
}
# "main" (m = 1/100) and "m200" (m = 1/200) are certified by cert.py: they are the content of
# the exit code.  The ladder entries are printed as NUMERICAL diagnostics (no interval
# certificate, so they cannot fail the exit code).
CERTIFIED_MODELS = ("main", "m200")
SETS = {
    "main": ["main", "m200"],
    "m200": ["m200"],
    "ladder": ["m400", "nk48", "nk64", "nw16", "m20", "lam6"],
    "all": ["main", "m200", "m400", "nk48", "nk64", "nw16", "m20", "lam6"],
}

MUTANTS = (
    "no_vertex",        # drop the vertex correction from the model itself
    "m_threshold",      # move to m = 2/100 (below the crossing)
    "rho_zero",         # verify a point instead of a box
    "res_scale",        # inflate the residual enclosure by 2^20 (r must then exceed rho)
    "rb_scale",         # inconsistent model: gRB scaled by 6/5 (box and PHI), PHI additionally
                        # scaled by 6/5 through box_sgn, so PHI is inflated by 36/25 vs DFb by 6/5
    "ray_thresh",       # demand lambda_max >= 1.01 instead of 1.005
    "diag_drop",        # under-report the diagonal of B by 2e-2 (Rayleigh bound fails)
    "drop_lb_t3",       # drop the vertex term from the linearised operator only
    "point_not_ball",   # certify a point instead of a box (radius 0, rho = 1e-12)
    "float_leak",       # lift the certified couplings from binary64 (radius-0 balls)
)


# ==================================================================================
#  float64 reference of the model (candidate selection and numerical self-tests)
# ==================================================================================
class FloatModel:
    def __init__(self, lam, m, Nk, Nw, eta, mu0, cx_zero=False):
        self.M = 2 * int(Nw)
        self.Nk, self.Nw = int(Nk), int(Nw)
        self.lam, self.m = float(lam), float(m)
        self.mu0, self.eta = float(mu0), float(eta)
        Nk_ = self.Nk
        self.xi = 2.0 * np.abs(2.0 * np.arange(Nk_) + 1.0 - Nk_) / Nk_ - 1.0 - self.mu0
        self.w = (2 * np.arange(-self.Nw, self.Nw) + 1).astype(float)
        self.tau = self.eta * float(C_INF) * math.sqrt(self.lam) * self.m
        self.nu = math.pi * self.tau
        self.gRB = 2.0 * self.lam * self.tau
        self.CX = 0.0 if cx_zero else (2.0 * self.lam * self.tau / Nk_) ** 2
        self.idxneg = (-np.arange(Nk_)) % Nk_
        self._s = None

    def K(self, d):
        return self.m ** 2 / (self.m ** 2 + (2.0 * math.pi * self.tau * d) ** 2)

    def kmat(self):
        M = self.M
        return [[self.K(n - a) for a in range(M)] for n in range(M)]

    def z_at(self, s):
        return 1.0 / (1j * self.w[:, None] * self.nu - self.xi[None, :] - s)

    def rb(self, F):
        km = self.kmat()
        kmean = np.array(km).dot(np.mean(F, axis=-1))
        return self.gRB * kmean[:, None] * np.ones(self.Nk)

    def sx(self, F, G, H):
        """Reference implementation (FFT), identical structure to verify_1d_counterexample.py."""
        M, Nk = self.M, self.Nk
        b = np.arange(M)[:, None, None] + np.arange(M)[None, :, None] - np.arange(M)[None, None, :]
        Ff = np.fft.fft(F, axis=1)
        Gf = np.fft.fft(G, axis=1)
        Hf = np.fft.fft(H[:, self.idxneg], axis=1)
        out_f = np.zeros((M, Nk), dtype=complex)
        for n in range(M):
            acc = np.zeros(Nk, dtype=complex)
            for a in range(M):
                bs = b[a, :, n]
                ok = (bs >= 0) & (bs < M)
                if not ok.any():
                    continue
                cs = np.nonzero(ok)[0]
                term = (np.array([self.K(n - c) for c in cs]))[:, None] * (Gf[cs] * Hf[bs[cs]])
                acc += self.K(n - a) * Ff[a] * term.sum(axis=0)
            out_f[n] = acc
        return self.CX * np.fft.ifft(out_f, axis=1)

    def F(self, s):
        z = self.z_at(s)
        return self.rb(z) + self.sx(z, z, z)

    def apply_full(self, v, s=None):
        """The model's linearised gap operator (no even projection)."""
        s = self._s if s is None else s
        v = np.asarray(v).reshape(self.M, self.Nk)
        z = self.z_at(s)
        az = np.abs(z)
        h = az * v
        return az * (self.rb(h) + 2.0 * self.sx(h, z, z).real + self.sx(z, np.conj(z), h).real)

    def me_profile(self):
        """Rainbow-only (Migdal-Eliashberg) self-energy; the physical branch start."""
        km = np.asarray(self.kmat())
        s = np.zeros((self.M, self.Nk), dtype=complex)
        for _ in range(4000):
            sn = self.gRB * (km @ np.mean(self.z_at(s), axis=1))[:, None] * np.ones(self.Nk)
            d = float(np.max(np.abs(sn - s)))
            s = 0.5 * sn + 0.5 * s
            if d < 1e-16:
                break
        return s

    def solve(self):
        # branch selection: start from the rainbow-only profile, i.e. the branch that is
        # continuously connected to the non-interacting state (this reproduces the audited
        # implementation's branch; see CERT.md section 'branch selection').
        s = self.me_profile()
        for _ in range(2000):
            sn = self.F(s)
            d = float(np.max(np.abs(sn - s)))
            s = 0.5 * sn + 0.5 * s
            if d < 1e-15:
                break
        km = self.kmat()
        D = self.M * self.Nk
        for _ in range(4):
            DF = closed_form(self.M, self.Nk, km, [list(r) for r in self.z_at(s)], self.gRB, self.CX, "jac")
            A = np.linalg.inv(np.eye(D) - DF)
            r = (self.F(s) - s).ravel()
            dx = (A @ r).reshape(self.M, self.Nk)
            s = s + dx
            if float(np.max(np.abs(dx))) < 1e-18:
                break
        self._s = s
        return s


# ==================================================================================
#  closed forms of Phi and of the two dense operators, in generic arithmetic
#  (Python complex for the float stage, flint acb for the certified stage)
# ==================================================================================
def tables(z, M, Nk):
    """Cyclic correlations of z (lists of complex floats or of acb):

        Ro[p][q][u]  = sum_m z(p,m) z(q,(m+u) mod Nk)
        Cv[a][c][s]  = sum_m z(a,m) z(c,(s-m) mod Nk)
        Cvc[a][c][s] = sum_m z(a,m) conj(z(c,(s-m) mod Nk))
    """
    zero = z[0][0] * 0
    zc = [[z[a][j].conjugate() for j in range(Nk)] for a in range(M)]
    Ro, Cv, Cvc = [], [], []
    for p in range(M):
        zp = z[p]
        r1, r2, r3 = [], [], []
        for q in range(M):
            zq, zqc = z[q], zc[q]
            ro = [zero] * Nk
            cv = [zero] * Nk
            cvc = [zero] * Nk
            for m in range(Nk):
                zpm = zp[m]
                for u in range(Nk):
                    ro[u] = ro[u] + zpm * zq[(m + u) % Nk]
                for s in range(Nk):
                    k = (s - m) % Nk
                    cv[s] = cv[s] + zpm * zq[k]
                    cvc[s] = cvc[s] + zpm * zqc[k]
            r1.append(ro)
            r2.append(cv)
            r3.append(cvc)
        Ro.append(r1)
        Cv.append(r2)
        Cvc.append(r3)
    return Ro, Cv, Cvc


def closed_form(M, Nk, Kmat, z, gRB, CX, kind, sgn=None):
    """Dense D x D matrix (D = M Nk), row (n,j), column (a0,j0).

      kind='jac' :  dPhi/dSigma                      (complex-linear Jacobian)
      kind='lin' :  the model's linearised gap operator B

    `z` may hold complex floats or flint acb intervals; a ball argument yields an
    enclosure of the whole family of operators (with `az = |z|` an interval too).
    """
    zero = z[0][0] * 0
    z2 = [[z[a][j] * z[a][j] for j in range(Nk)] for a in range(M)]
    az = [[abs(z[a][j]) for j in range(Nk)] for a in range(M)]
    Ro, Cv, Cvc = tables(z, M, Nk)
    # exact 1/Nk: in the certified (acb) stage a binary64 lift would not contain 1/Nk
    invNk = (arb(1) / arb(Nk)) if isinstance(zero, acb) else 1.0 / Nk
    if sgn is None or sgn.get("fac") is None:
        fac = arb(1) if isinstance(zero, acb) else 1.0
    else:
        fnum, fden = sgn["fac"]
        fac = (arb(fnum) / arb(fden)) if isinstance(zero, acb) else fnum / fden
    drop_t3 = bool(sgn is not None and sgn.get("drop_t3", False))
    D = M * Nk
    if isinstance(zero, acb):
        out = acb_mat(D, D)
    elif kind == "jac":
        out = np.zeros((D, D), dtype=np.complex128)
    else:
        out = np.zeros((D, D), dtype=np.float64)
    for n in range(M):
        Kr = Kmat[n]
        for a0 in range(M):
            Kn = Kr[a0]
            v2 = [(p, p + a0 - n) for p in range(M) if 0 <= p + a0 - n < M]
            v3 = [(a, a0 - a + n) for a in range(M) if 0 <= a0 - a + n < M]
            S2 = [zero] * Nk
            S3 = [zero] * Nk
            S3c = [zero] * Nk
            for (p, b) in v2:
                w, rp = Kr[p], Ro[p][b]
                for u in range(Nk):
                    S2[u] = S2[u] + w * rp[u]
            if not drop_t3:
                for (a, c) in v3:
                    w, cva, cvc = Kr[a] * Kr[c], Cv[a][c], Cvc[a][c]
                    for s in range(Nk):
                        S3[s] = S3[s] + w * cva[s]
                        S3c[s] = S3c[s] + w * cvc[s]
            for j0 in range(Nk):
                for j in range(Nk):
                    t12 = gRB * Kn * invNk * fac + 2.0 * CX * Kn * S2[(j0 - j) % Nk]
                    if kind == "jac":
                        out[n * Nk + j, a0 * Nk + j0] = z2[a0][j0] * (t12 + CX * S3[(j0 + j) % Nk])
                    else:
                        val = az[a0][j0] * t12 + CX * az[a0][j0] * S3c[(j0 + j) % Nk]
                        out[n * Nk + j, a0 * Nk + j0] = az[n][j] * val.real
    return out


def phi_closed(M, Nk, Kmat, z, gRB, CX, sgn=None):
    """Phi in generic arithmetic; returns an M x Nk list of lists."""
    zero = z[0][0] * 0
    if sgn is None or sgn.get("fac") is None:
        fac = arb(1) if isinstance(zero, acb) else 1.0
    else:
        fnum, fden = sgn["fac"]
        fac = (arb(fnum) / arb(fden)) if isinstance(zero, acb) else fnum / fden
    mk = [sum((z[a][j] for j in range(Nk)), zero) / Nk for a in range(M)]
    Ro, _, _ = tables(z, M, Nk)
    out = [[zero] * Nk for _ in range(M)]
    for n in range(M):
        rb = sum((Kmat[n][a] * mk[a] for a in range(M)), zero) * gRB * fac
        for j in range(Nk):
            acc = zero
            for a in range(M):
                for c in range(M):
                    b = a + c - n
                    if not (0 <= b < M):
                        continue
                    s = zero
                    for m in range(Nk):
                        s = s + z[a][m] * Ro[c][b][(m - j) % Nk]
                    acc = acc + Kmat[n][a] * Kmat[n][c] * s
            out[n][j] = rb + CX * acc
    return out


# ==================================================================================
#  rigorous upper/lower bounds for the imaginary part of an acb entry
#  (needed to certify that no propagator denominator touches the real axis)
# ==================================================================================
def im_hi(x):
    """Rigorous arb upper bound on |x| (verifies its own rounding)."""
    u = x.abs_upper()
    try:
        e = math.frexp(float(u))[1]
    except (OverflowError, ValueError):  # pragma: no cover
        e = 0
    for _ in range(80):
        b = arb(2.0) ** e
        if bool(u <= b):
            return b
        e += 1
    raise RuntimeError("im_hi: no enclosure found")  # pragma: no cover


def imag_upper(z):
    return im_hi(z.imag)


def abs_lower(x):
    """Rigorous arb lower bound on |x| (0 when the sign of x is not decided).

    x.lower()/x.upper() are the rounded-down/up endpoints of the enclosure, so
    x > 0  =>  |x| >= x.lower()  and  x < 0  =>  |x| >= -x.upper()."""
    lo = arb(x.lower())      # lower endpoint of the enclosure
    if bool(lo > 0):
        return lo
    hi = arb(-x.upper())     # lower bound on -x, i.e. a bound on |x| when x < 0
    if bool(hi > 0):
        return hi
    return arb(0)


def imag_lower(z):
    """Rigorous arb lower bound on |Im z|."""
    return abs_lower(z.imag)


# ==================================================================================
#  decimal reporting of rigorous bounds (never claim more than the ball proves)
# ==================================================================================
def _fmt(d, sig, up):
    """Format the finite float d with sig significant digits, rounded down (up=False)
    or up (up=True).  The result is a decimal point q with q <= d (or q >= d) exactly."""
    from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
    if d == 0 or not math.isfinite(d):
        return repr(d)
    e = math.floor(math.log10(abs(d)))
    places = max(0, sig - 1 - e)
    q = Decimal(d).quantize(Decimal(1).scaleb(-places),
                            rounding=ROUND_CEILING if up else ROUND_FLOOR)
    return str(q)


def lower_str(x, sig=12):
    """Decimal string q with q <= (lower endpoint of the arb ball x)."""
    d = float(x.lower())
    if not math.isfinite(d):
        return repr(d)
    return _fmt(d - 8.0 * abs(math.ulp(d)), sig, False)


def upper_str(x, sig=12):
    """Decimal string q with q >= (upper endpoint of the arb ball x)."""
    d = float(x.upper())
    if not math.isfinite(d):
        return repr(d)
    return _fmt(d + 8.0 * abs(math.ulp(d)), sig, True)


# ==================================================================================
#  interval model
# ==================================================================================
class ArbModel:
    def __init__(self, lam, m_rat, eta, mu0_rat, Nk, Nw, cx_zero=False, fac_one=True):
        ctx.prec = PREC
        self.M = 2 * int(Nw)
        self.Nk, self.Nw = int(Nk), int(Nw)
        Nk_ = self.Nk
        self.lam = arb(lam)
        self.m = arb(m_rat[0]) / arb(m_rat[1])
        self.eta = arb(eta) if not isinstance(eta, arb) else eta
        self.mu0 = arb(mu0_rat[0]) / arb(mu0_rat[1])
        self.c_inf = arb(C_INF)
        self.tau = self.eta * self.c_inf * self.lam.sqrt() * self.m
        self.nu = arb.pi() * self.tau
        self.w = [arb(2 * n + 1) for n in range(-self.Nw, self.Nw)]
        self.xi = [arb(2 * abs(2 * j + 1 - Nk_)) / arb(Nk_) - arb(1) - self.mu0 for j in range(Nk_)]
        two_pi_tau = arb(2) * arb.pi() * self.tau
        m2 = self.m * self.m
        self.Kmat = []
        for n in range(self.M):
            row = []
            for a in range(self.M):
                if n == a:
                    row.append(arb(1))
                else:
                    d = arb(abs(n - a))
                    row.append(m2 / (m2 + (two_pi_tau * d) * (two_pi_tau * d)))
            self.Kmat.append(row)
        self.gRB = arb(2) * self.lam * self.tau
        self.CX = arb(0) if cx_zero else (arb(2) * self.lam * self.tau / arb(Nk_)) ** 2
        self.fac = fac_one

    def z_ball(self, sball):
        """z over a complex box sball (list of lists of acb); returns (z, min|Im u|)."""
        z = []
        min_imu = None
        for n in range(self.M):
            row = []
            for j in range(self.Nk):
                u = acb(-self.xi[j], self.w[n] * self.nu) - sball[n][j]
                row.append(acb(1) / u)
                li = imag_lower(u)
                # arb-overlap-safe minimum: never drop a smaller lower bound
                min_imu = li if min_imu is None else min_imu.min(li)
            z.append(row)
        return z, min_imu

    def phi(self, z, sgn=None):
        return phi_closed(self.M, self.Nk, self.Kmat, z, self.gRB, self.CX, sgn)

    def dphi(self, z, sgn=None):
        return closed_form(self.M, self.Nk, self.Kmat, z, self.gRB, self.CX, "jac", sgn)

    def lb(self, z, sgn=None):
        return closed_form(self.M, self.Nk, self.Kmat, z, self.gRB, self.CX, "lin", sgn)


# ==================================================================================
#  helpers
# ==================================================================================
def ball_of(s0, rho_s):
    """Complex box  s0(n,j) +- rho   (rho an exact decimal radius) as an acb rectangle.

    Each component is the interval hull of the two endpoints s0 -+ rho, so the ball
    really contains the whole closed box of radius rho (the v1 helper returned the
    *shifted* point s0 + rho with a rounding-size radius, which the box guard below
    now rejects).
    """
    pm = arb(rho_s)
    out = []
    for n in range(s0.shape[0]):
        row = []
        for j in range(s0.shape[1]):
            re = arb(s0[n][j].real)
            im_ = arb(s0[n][j].imag)
            row.append(acb((re - pm).union(re + pm), (im_ - pm).union(im_ + pm)))
        out.append(row)
    return out


def point_of(s0):
    """Degenerate "box": the exact double point Sigma_0, radius 0 (mutant support)."""
    return [[acb(arb(s0[n][j].real), arb(s0[n][j].imag)) for j in range(s0.shape[1])]
            for n in range(s0.shape[0])]


def from_np(Mf):
    """Exact lift of a float complex matrix to an acb matrix (doubles are exact in arb)."""
    D = Mf.shape[0]
    A = acb_mat(D, D)
    for i in range(D):
        for j in range(D):
            A[i, j] = acb(arb(Mf[i, j].real), arb(Mf[i, j].imag))
    return A


def eye(D):
    A = acb_mat(D, D)
    for i in range(D):
        A[i, i] = acb(1)
    return A


def rowsum_inf(A, D):
    """Rigorous arb upper bound on the infinity (max row sum) norm of an acb matrix."""
    hi = arb(0)
    for i in range(D):
        s = arb(0)
        for j in range(D):
            s = s + A[i, j].abs_upper()
        hi = hi.max(s)          # arb-overlap-safe maximum
    return hi


# ==================================================================================
#  the certificate for one model
# ==================================================================================
def certify_one(name, lam, m_rat, eta, mu0_rat, Nk, Nw, mutate="none", rho_s=RHO,
                theta=THETA, prec=PREC, verbose=True, fid=True, interval=True):
    ctx.prec = prec
    rep = {"model": name, "Nk": Nk, "Nw": Nw, "lam": lam, "m": "%d/%d" % m_rat,
           "eta": eta, "mu0": "%d/%d" % mu0_rat, "mutate": mutate, "checks": [], "fail": []}
    t0 = time.time()
    M = 2 * Nw
    D = M * Nk
    rep["M"], rep["D"] = M, D

    def emit(s):
        if verbose:
            print(s, flush=True)

    emit("[model %s] Nk=%d Nw=%d  M=%d D=%d  lam=%s m=%d/%d eta=%s  mutate=%s  prec=%d bits"
         % (name, Nk, Nw, M, D, lam, m_rat[0], m_rat[1], eta, mutate, prec))

    # ------------------------------------------------------------------ mutants
    cx_zero = mutate == "no_vertex"
    m_float = 2.0 / 100 if mutate == "m_threshold" else _r(m_rat)
    rho_use = "0" if mutate == "rho_zero" else rho_s
    theta_use = 1e-2 if mutate == "ray_thresh" else theta
    lb_sgn = {"drop_t3": True} if mutate == "drop_lb_t3" else None
    box_sgn = {"fac": (6, 5)} if mutate == "rb_scale" else None
    if mutate == "m_threshold":
        m_rat = (2, 100)
        rep["m"] = "2/100"

    # --------------------------------------------------------- 1) float stage
    tf = time.time()
    fm = FloatModel(lam, m_float, Nk, Nw, eta, _r(mu0_rat), cx_zero=cx_zero)
    s0 = fm.solve()
    res_f = float(np.max(np.abs(fm.F(s0) - s0)))
    km = fm.kmat()
    zf = [list(r) for r in fm.z_at(s0)]
    emit("[1/6] float candidate Sigma_0:  ||Phi(Sigma_0) - Sigma_0||_inf = %.3e   (%.1fs)"
         % (res_f, time.time() - tf))
    rep["float_residual"] = res_f
    rep["float_stage_s"] = time.time() - tf

    # model fidelity (numerical, not part of the certified decision)
    fid_rep = {}
    if fid:
        tf = time.time()
        ph_f = np.array([[complex(x) for x in r] for r in phi_closed(M, Nk, km, zf, fm.gRB, fm.CX)])
        e1 = float(np.max(np.abs(ph_f - fm.F(s0))))
        DF_f = closed_form(M, Nk, km, zf, fm.gRB, fm.CX, "jac")
        LB_f = closed_form(M, Nk, km, zf, fm.gRB, fm.CX, "lin")
        rng = np.random.default_rng(11)
        e2 = 0.0
        for _ in range(3):
            v = rng.standard_normal((M, Nk))
            e2 = max(e2, float(np.max(np.abs(LB_f @ v.ravel() - fm.apply_full(v).ravel()))))
        e3 = 0.0
        h = 1e-7
        for _ in range(6):
            a0 = int(rng.integers(M))
            j0 = int(rng.integers(Nk))
            S1 = s0.copy()
            S1[a0, j0] += h
            S2 = s0.copy()
            S2[a0, j0] -= h
            d = ((np.array([[complex(x) for x in r] for r in phi_closed(M, Nk, km, [list(q) for q in fm.z_at(S1)], fm.gRB, fm.CX)])
                  - np.array([[complex(x) for x in r] for r in phi_closed(M, Nk, km, [list(q) for q in fm.z_at(S2)], fm.gRB, fm.CX)])) / (2 * h)).ravel()
            e3 = max(e3, float(np.max(np.abs(d - DF_f[:, a0 * Nk + j0]))))
        S = 0.5 * (LB_f + LB_f.T)
        w, V = np.linalg.eigh(S)
        v = np.ascontiguousarray(V[:, -1])
        lam_float = float(w[-1])
        rq_float = float(v @ LB_f @ v / (v @ v))
        fid_rep = {"closed_vs_fft_phi": e1, "closed_vs_apply_full": e2, "jac_vs_finite_difference": e3,
                   "lambda_max_float": lam_float, "rayleigh_float": rq_float,
                   "antisym_float": float(np.max(np.abs(LB_f - LB_f.T))),
                   "symmetry_float": float(np.max(np.abs(S - S.T))),
                   "min_abs_imag_u": float(np.min(np.abs(fm.w)) * fm.nu - np.max(np.abs(s0.imag))),
                   "secs": time.time() - tf}
        LB_f_jac = DF_f
        if e1 > 1e-12 or e2 > 1e-12 or e3 > 1e-5:
            rep["fail"].append("model transcription fidelity")
        emit("[2/6] model fidelity (numerical): closed-form Phi vs FFT %.2e ; B vs apply_full %.2e ; "
             "Jacobian vs FD %.2e ; float lambda_max %.9f ; ||B-B^T||_inf %.2e"
             % (e1, e2, e3, lam_float, fid_rep["antisym_float"]))
    else:
        LB_f = closed_form(M, Nk, km, zf, fm.gRB, fm.CX, "lin")
        LB_f_jac = closed_form(M, Nk, km, zf, fm.gRB, fm.CX, "jac")
        S = 0.5 * (LB_f + LB_f.T)
        w, V = np.linalg.eigh(S)
        v = np.ascontiguousarray(V[:, -1])
        lam_float = float(w[-1])
        rq_float = lam_float
        fid_rep = {"lambda_max_float": lam_float, "rayleigh_float": rq_float}
    rep["fidelity"] = fid_rep

    if fid:
        # independent cross-check against the reference implementation of the model
        try:
            import importlib.util
            path = str(_repo_root() / "verify_1d_counterexample.py")
            spec = importlib.util.spec_from_file_location("vertex_ref", path)
            ref = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(ref)
            out = ref.solve_finite_m_lattice(lam=float(lam), m_omega=float(m_rat[0]) / m_rat[1],
                                             eta=float(eta), mu0=_r(mu0_rat), N_k=Nk, N_w=Nw)
            dl = abs(out["full"] - lam_float)
            fid_rep["reference_lambda"] = float(out["full"])
            fid_rep["reference_me"] = float(out["me"])
            fid_rep["reference_dlambda"] = dl
            emit("[2b/6] reference solve_finite_m_lattice: lambda_full %.9f (delta %.1e), lambda_me %.9f"
                 % (out["full"], dl, out["me"]))
            if dl > 1e-5:
                rep["fail"].append("reference lambda mismatch")
        except Exception as exc:  # pragma: no cover
            emit("[2b/6] reference cross-check unavailable: %s" % exc)
            fid_rep["reference_error"] = str(exc)

    if not interval:
        ok_f = (res_f <= 1e-12) and (lam_float >= 1.0 + theta_use)
        rep["float_only"] = True
        rep["checks"].append(("float_residual_small", res_f <= 1e-12))
        rep["checks"].append(("float_lambda_above_threshold", lam_float >= 1.0 + theta_use))
        if not ok_f:
            rep["fail"].append("float criterion (numerical stage)")
        rep["ok"] = not rep["fail"]
        rep["secs"] = time.time() - t0
        emit("[3/3] float-only: residual %.2e ; lambda_max %.9f vs 1 + %.1e -> %s ; wall %.1fs"
             % (res_f, lam_float, theta_use, ok_f, rep["secs"]))
        return rep

    # ------------------------------------------------------- 2) interval stage
    ti = time.time()
    am = ArbModel(lam, m_rat, eta, mu0_rat, Nk, Nw, cx_zero=cx_zero)
    if mutate == "rb_scale":
        am.gRB = am.gRB * arb(6) / arb(5)
    if mutate == "float_leak":
        am.gRB = arb(float(am.gRB))          # binary64 lift: an exact, radius-0 ball
        am.CX = arb(float(am.CX))
    sball = point_of(s0) if mutate == "point_not_ball" else ball_of(s0, rho_use)

    # ---------------- guards: the ball really is a ball, and the couplings are not floats
    rho_half = arb(rho_use) / arb(2)     # loose threshold: the guard separates rho from 0
    box_ok = True
    for n in range(M):
        for j in range(Nk):
            if not (bool(sball[n][j].real.rad() >= rho_half)
                    and bool(sball[n][j].imag.rad() >= rho_half)):
                box_ok = False
                break
    leak_ok = bool(am.gRB.rad() > 0) and bool(am.CX.rad() > 0)
    emit("      guards: ball radius >= rho in every component -> %s ; "
         "couplings computed in arb (radius > 0) -> %s" % (box_ok, leak_ok))
    if not box_ok:
        rep["fail"].append("box guard: the enclosure is not a ball of radius rho")
    if not leak_ok:
        rep["fail"].append("constant guard: a certified coupling was lifted from binary64")
    rep["checks"].append(("box_guard", box_ok))
    rep["checks"].append(("coupling_guard", leak_ok))
    zball, min_imu = am.z_ball(sball)
    rep["min_abs_im_denominator"] = lower_str(min_imu)
    if not bool(min_imu > 0):
        rep["fail"].append("propagator pole inside the box")

    # residual enclosure  A (Phi(Sigma_0) - Sigma_0)   -- evaluated AT the point Sigma_0,
    # so that the radius rho of the box does not enter the residual bound at all
    s0c = [[acb(arb(s0[n, j].real), arb(s0[n, j].imag)) for j in range(Nk)] for n in range(M)]
    z0c, _ = am.z_ball(s0c)
    PHI = am.phi(z0c, box_sgn)
    Dm = 1 << 20 if mutate == "res_scale" else 1
    res = [[(PHI[n][j] - s0c[n][j]) * Dm for j in range(Nk)] for n in range(M)]
    DFb = am.dphi(zball)
    DF0 = from_np(LB_f_jac)
    emit("[3/6] interval enclosure: Phi, DPhi over the box (%.1fs);  min |Im u| >= %s"
         % (time.time() - ti, rep["min_abs_im_denominator"]))

    # preconditioner A = inv(I - DPhi(Sigma_0)) from the float candidate (exact doubles)
    Amat = from_np(np.linalg.inv(np.eye(D) - LB_f_jac))
    rmat = acb_mat(D, 1)
    for i in range(D):
        rmat[i, 0] = res[i // Nk][i % Nk]
    rvec = Amat * rmat
    r_hi = arb(0)
    for i in range(D):
        r_hi = r_hi.max(rvec[i, 0].abs_upper())   # arb-overlap-safe maximum
    # E = I - A (I - DF(box))   in interval arithmetic
    E = eye(D) - Amat * (eye(D) - DFb)
    k_hi = rowsum_inf(E, D)
    # the same with the sum split (bound on the ball radius contribution) -- reported only
    k_split = rowsum_inf(Amat * (DFb - DF0), D)
    rep["contraction"] = {"r": upper_str(r_hi), "k": upper_str(k_hi), "k_split": upper_str(k_split)}
    rho_arb = arb(rho_use)
    cond_map = bool(r_hi + k_hi * rho_arb <= rho_arb)
    cond_ctr = bool(k_hi < 1)
    rep["krawczyk"] = {"maps_into_box": cond_map, "contraction": cond_ctr,
                       "fixed_point_radius": upper_str(r_hi / (1 - k_hi)) if cond_ctr else None}
    emit("[4/6] Krawczyk/Banach: ||A(Phi(S0)-S0)|| <= %s ; ||I - A(I - DPhi(box))|| <= %s "
         "(split bound %s)" % (upper_str(r_hi), upper_str(k_hi), upper_str(k_split)))
    emit("      criterion: r + k*rho <= rho  (%.3e <= %.3e) -> %s ;  k < 1 -> %s ;  unique fixed point "
         "in ||Sigma-Sigma_0||_inf <= %s"
         % (float(r_hi + k_hi * rho_arb), float(rho_arb), cond_map, cond_ctr,
            upper_str(r_hi / (1 - k_hi)) if cond_ctr else "inf"))
    rep["checks"].append(("krawczyk_maps_box", cond_map))
    rep["checks"].append(("krawczyk_contraction", cond_ctr))
    if not (cond_map and cond_ctr):
        rep["fail"].append("Krawczyk criterion")

    # ------------------------------------------------------- 3) Rayleigh stage
    Lmat = am.lb(zball, lb_sgn)
    if mutate == "diag_drop":
        for i in range(D):
            Lmat[i, i] = Lmat[i, i] - acb(arb("2e-2"))
    # B(Sigma) is identically real symmetric for every Sigma (CERT.md section 3, Prop. 1),
    # so the symmetric part is an exact enclosure of B over the box and the plain symmetric
    # Rayleigh quotient is a rigorous lower bound on lambda_max(B).  The antisymmetric part
    # of the *enclosure* is then only a diagnostic of the box width; it is bounded in a guard
    # below and is never subtracted (the shrunken bound of cert.v1.py is unsound, see REVIEW3).
    Ssym = (Lmat + Lmat.transpose()) * acb(arb(1) / arb(2))
    vc = acb_mat(D, 1)
    den = arb(0)
    for i in range(D):
        vv = arb(float(v[i]))
        vc[i, 0] = acb(vv)
        den = den + vv * vv
    Sv = Ssym * vc
    num = acb(0)
    for i in range(D):
        num = num + acb(arb(float(v[i]))) * Sv[i, 0]
    q_lo = num.real / den
    # diagnostic only: the antisymmetric part of the *enclosure* (not of B) must stay at the
    # width of the box, i.e. the enclosure itself confirms B = B^T away from the ball width
    delta_hi = arb(0)
    for i in range(D):
        s = arb(0)
        for j in range(D):
            s = s + (Lmat[i, j] - Lmat[j, i]).abs_upper()
        delta_hi = delta_hi.max(s)          # arb-overlap-safe maximum
    sym_ok = bool(delta_hi / arb(2) <= arb("1e-6"))
    bound = q_lo
    rep["rayleigh"] = {"q_lo": lower_str(q_lo), "delta_sh": upper_str(delta_hi / arb(2), 4),
                       "bound": lower_str(bound), "theta": theta_use, "symmetry_diagnostic": sym_ok}
    ok_ray = bool(bound >= arb(1) + arb(repr(theta_use)))
    emit("[5/6] Rayleigh: <v,B v>/<v,v> >= %s ;  antisym part of enclosure <= %s (diagnostic %s)"
         % (lower_str(q_lo), upper_str(delta_hi / 2, 4), sym_ok))
    emit("      criterion: bound >= 1 + %.1e  -> %s (margin >= %s)"
         % (theta_use, ok_ray, lower_str(bound - arb(1) - arb(repr(theta_use)))))
    rep["checks"].append(("rayleigh_bound", ok_ray))
    rep["checks"].append(("symmetry_diagnostic", sym_ok))
    if not ok_ray:
        rep["fail"].append("Rayleigh lower bound")
    if not sym_ok:
        rep["fail"].append("symmetry diagnostic")

    rep["secs"] = time.time() - t0
    rep["ok"] = not rep["fail"]
    emit("[6/6] %s  (%s)  wall %.1fs"
         % ("CERTIFIED" if rep["ok"] else "FAILED", ", ".join(rep["fail"]) or "all inequalities hold",
            rep["secs"]))
    return rep


def _r(rat):
    return float(rat[0]) / float(rat[1])


def _repo_root():
    """Directory holding verify_1d_counterexample.py (cert.py sits in its CERT/ subdir)."""
    import os
    from pathlib import Path
    env = os.environ.get("VERTEX_UNIV")
    if env:
        return Path(env)
    return Path(__file__).resolve().parent.parent


# ==================================================================================
#  driver
# ==================================================================================
def run(setname="main", models=None, mutate="none", rho=RHO, theta=THETA, prec=PREC,
        verbose=True, json_out=None, fid=True):
    if mutate not in ("none",) + MUTANTS:
        sys.exit("unknown mutation %r (choose from %s)" % (mutate, ", ".join(("none",) + MUTANTS)))
    names = models if models else SETS[setname]
    reps = []
    ok = True
    for nm in names:
        lam, m_rat, eta, mu0_rat, Nk, Nw = MODELS[nm]
        rep = certify_one(nm, lam, m_rat, eta, mu0_rat, Nk, Nw, mutate=mutate, rho_s=rho,
                          theta=theta, prec=prec, verbose=verbose, fid=fid,
                          interval=(nm in CERTIFIED_MODELS))
        reps.append(rep)
        if nm in CERTIFIED_MODELS:
            ok = ok and rep["ok"]
        else:
            if verbose:
                print("[%s] NUMERICAL diagnostic only (not part of the certified exit code): %s"
                      % (nm, "criteria met" if rep["ok"] else "criteria not met: " + ", ".join(rep["fail"])),
                      flush=True)
    if json_out:
        with open(json_out, "w", encoding="utf-8") as fh:
            json.dump({"mutate": mutate, "rho": rho, "theta": theta, "prec": prec,
                       "models": reps, "ok": ok}, fh, indent=2)
    if verbose and len(reps) > 1:
        print("\n=== summary (%s) ===" % mutate)
        for rep in reps:
            cert = rep["model"] in CERTIFIED_MODELS
            bound = rep.get("rayleigh", {}).get("bound", rep.get("fidelity", {}).get("lambda_max_float"))
            bound_s = bound if isinstance(bound, str) else "%s" % ("nan" if bound is None else "%.10f" % bound)
            print("%-6s Nk=%-3d Nw=%-2d D=%-4d lambda_max %s %-14s   %s"
                  % (rep["model"], rep["Nk"], rep["Nw"], rep["D"], ">=" if cert else "~",
                     bound_s,
                     ("CERTIFIED" if cert else "NUMERICAL") if rep["ok"]
                     else "FAILED: " + ", ".join(rep["fail"])))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", default="main", choices=sorted(SETS))
    ap.add_argument("--models", default=None, help="comma separated model names (overrides --set)")
    ap.add_argument("--mutate", default="none", choices=("none",) + MUTANTS)
    ap.add_argument("--rho", default=RHO)
    ap.add_argument("--theta", type=float, default=THETA)
    ap.add_argument("--prec", type=int, default=PREC)
    ap.add_argument("--json", default=None)
    ap.add_argument("--no-fidelity", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    models = a.models.split(",") if a.models else None
    return run(setname=a.set, models=models, mutate=a.mutate, rho=a.rho, theta=a.theta,
               prec=a.prec, verbose=not a.quiet, json_out=a.json, fid=not a.no_fidelity)


if __name__ == "__main__":
    sys.exit(main())
