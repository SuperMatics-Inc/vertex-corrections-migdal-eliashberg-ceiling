"""Exact-symmetric 1D/flat-sheet first-crossed-vertex counterexample verification.

Computes:
1. Continuum m -> 0 scaling-limit critical eigenvalue (Full vs ME) via the exact-symmetric
   lumped-tail L^2(|z|^2 dy) self-adjoint operator (sym_err < 1e-12 asserted).
2. lambda-continuation at eta = 1.000 and eta_crit root-finding at lambda = 4, 6, 8.
3. Continuum convergence table over (N_w, Y_max, N_y, dy) and Richardson-style extrapolation.
4. Exact-symmetric finite-m periodic-lattice sequence approaching the continuum from below,
   with bare RPA 2k_F charge susceptibility diagnostic and a CERTIFIED (padded binary64)
   50-digit Decimal inner product + conservative IEEE-754 forward-error Rayleigh-quotient
   lower bound on the fixed finite-(m, Nk, Nw) linearized operator at Sigma_num.
"""

from __future__ import annotations

import gc
import time
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator, eigsh

C_inf = 0.18272624777232952


def anderson_solve(F, x0, tol=1e-10, maxiter=60, m=8, beta=0.4):
    x = x0.copy()
    X_hist, R_hist = [], []
    for it in range(maxiter):
        fx = F(x)
        r = fx - x
        nr = float(np.max(np.abs(r)))
        if nr < tol:
            return x, nr, it
        X_hist.append(x.copy())
        R_hist.append(r.copy())
        if len(X_hist) > m:
            X_hist.pop(0)
            R_hist.pop(0)
        k = len(X_hist)
        if k == 1:
            x = x + beta * r
        else:
            dR = np.column_stack([R_hist[i + 1] - R_hist[i] for i in range(k - 1)])
            dX = np.column_stack([X_hist[i + 1] - X_hist[i] for i in range(k - 1)])
            gamma, *_ = np.linalg.lstsq(dR, R_hist[-1], rcond=1e-12)
            x = X_hist[-1] - dX @ gamma + beta * (R_hist[-1] - dR @ gamma)
    r = F(x) - x
    return x, float(np.max(np.abs(r))), maxiter


def make_1d_scaling_solver(N_w=12, Ymax=120.0, N_y=241):
    L = N_y // 2
    N_y = 2 * L + 1
    dy = 2.0 * Ymax / (N_y - 1)
    sqrt_dy = np.sqrt(dy)
    y = np.arange(-L, L + 1) * dy
    M = 2 * N_w
    n_idx = np.arange(-N_w, N_w)
    w_n = 2 * n_idx + 1
    N_fft = 1 << int(np.ceil(np.log2(6 * L + 4)))
    M_pad = 1 << int(np.ceil(np.log2(2 * M)))
    a_idx = np.arange(M)
    b_idx = np.arange(M)
    w_t = np.ones(N_y) * dy
    w_t[0] *= 0.5
    w_t[-1] *= 0.5
    idx_rev = 4 * L - np.arange(N_y)
    # Uniform-dy interior quadrature paired with exact-symmetric lumped tail degree of freedom.
    # Let S_p = sum_i sqrt(dy) * h_{p,i} + sqrt(tail_p) * d_p^tail, where h = |z| d_int.
    # The normal equation itself uses trapezoidal + tail integral for Sigma(n, y).

    def run(lam=8.0, eta=1.0, alpha_f=0.5, x_init=None, eig_tol=1e-8, check_sym=True):
        t0 = time.time()
        alpha_b = 1.0 - alpha_f
        w1 = 4.0 * alpha_f**2
        w23 = 4.0 * alpha_f * alpha_b
        t_val = eta * C_inf * np.sqrt(lam)
        ip_g, ia_g = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
        j_arr = np.arange(-(M - 1), M)
        K_mat = (1.0 / (1.0 + (2.0 * np.pi * t_val * j_arr) ** 2))[ip_g - ia_g + (M - 1)]

        def rb_op_normal(s):
            z = 1.0 / (1j * w_n[:, None] - y[None, :] - s)
            int_y = z @ w_t
            u_inf = w_n - 0.5 * (s[:, 0].imag + s[:, -1].imag)
            tail = -2j * np.sign(u_inf) * (0.5 * np.pi - np.arctan(Ymax / np.abs(u_inf)))
            return (lam / np.pi) * (K_mat @ (int_y + tail))[:, None] * np.ones((1, N_y))
        def cx_op_nonlinear(z):
            zf = np.fft.fft(z, n=N_fft, axis=1)
            z_rev_f = np.fft.fft(z[:, ::-1], n=N_fft, axis=1)
            pref = (lam / (2.0 * np.pi)) ** 2 * (dy**2)
            out1 = np.zeros((M, N_fft), dtype=complex)
            out23 = np.zeros((M, N_fft), dtype=complex)
            for p in range(M):
                U = K_mat[p, :, None] * zf
                Ur = K_mat[p, :, None] * z_rev_f
                Uh = np.fft.fft(U, n=M_pad, axis=0)
                Urh = np.fft.fft(Ur, n=M_pad, axis=0)
                idx_b = p + b_idx
                c1 = np.fft.ifft(Uh * Uh, axis=0)[idx_b]
                c23 = 2.0 * np.fft.ifft(Urh * Uh, axis=0)[idx_b]
                out1[p] = np.sum(c1 * z_rev_f, axis=0)
                out23[p] = np.sum(c23 * z_rev_f, axis=0)
            t1 = np.fft.ifft(out1, axis=1)[:, 2 * L : 2 * L + N_y]
            t23 = np.fft.ifft(out23, axis=1)[:, idx_rev]
            return pref * (w1 * t1 + w23 * t23)

        # Exact ME benchmark on same Matsubara cutoff
        u_me = w_n.astype(float).copy()
        for _ in range(100):
            u_me_n = w_n + lam * (K_mat @ np.sign(u_me))
            if np.max(np.abs(u_me_n - u_me)) < 1e-13:
                u_me = u_me_n
                break
            u_me = 0.5 * u_me_n + 0.5 * u_me
        inv_su = 1.0 / np.sqrt(np.abs(u_me))
        vals_me, vecs_me = np.linalg.eigh(lam * inv_su[:, None] * K_mat * inv_su[None, :])
        me_eig = float(vals_me[-1])
        v_me_freq = vecs_me[:, -1]

        if x_init is None:
            sig_me = np.zeros((M, N_y), dtype=complex)
            for _ in range(80):
                sn = rb_op_normal(sig_me)
                if np.max(np.abs(sn - sig_me)) < 1e-11:
                    sig_me = sn
                    break
                sig_me = 0.5 * sn + 0.5 * sig_me
            x_init = np.concatenate([sig_me.real.ravel(), sig_me.imag.ravel()])

        def eval_F_vec(xr):
            s = (xr[: M * N_y] + 1j * xr[M * N_y :]).reshape(M, N_y)
            z = 1.0 / (1j * w_n[:, None] - y[None, :] - s)
            fn = rb_op_normal(s) + cx_op_nonlinear(z)
            return np.concatenate([fn.real.ravel(), fn.imag.ravel()])

        x_sol, res, n_it = anderson_solve(eval_F_vec, x_init, tol=1e-10, maxiter=60, m=8, beta=0.4)
        sig = (x_sol[: M * N_y] + 1j * x_sol[M * N_y :]).reshape(M, N_y)
        min_u = float(np.min((w_n[:, None] - sig.imag) / w_n[:, None]))
        z = 1.0 / (1j * w_n[:, None] - y[None, :] - sig)
        az = np.abs(z)
        u_inf = np.abs(w_n - 0.5 * (sig[:, 0].imag + sig[:, -1].imag))
        tail_fac = 2.0 * (0.5 * np.pi - np.arctan(Ymax / u_inf)) / u_inf
        sqrt_tail = np.sqrt(tail_fac)

        # Precompute linear crossed-vertex response matrices for fast exact-symmetric matvec
        zf = np.fft.fft(z, n=N_fft, axis=1)
        z_rev_f = np.fft.fft(z[:, ::-1], n=N_fft, axis=1)
        zcf = np.fft.fft(np.conj(z), n=N_fft, axis=1)
        zcf_rev = np.fft.fft(np.conj(z)[:, ::-1], n=N_fft, axis=1)
        z_rev_flip_h = np.fft.fft(z_rev_f[::-1], n=M_pad, axis=0)
        A1 = np.zeros((M, M, N_fft), dtype=complex)
        A3 = np.zeros((M, M, N_fft), dtype=complex)
        B1 = np.zeros((M, M, N_fft), dtype=complex)
        B23 = np.zeros((M, M, N_fft), dtype=complex)
        for p in range(M):
            U = K_mat[p, :, None] * zf
            Ur = K_mat[p, :, None] * z_rev_f
            Uc = K_mat[p, :, None] * zcf
            Ucr = K_mat[p, :, None] * zcf_rev
            Uh = np.fft.fft(U, n=M_pad, axis=0)
            Urh = np.fft.fft(Ur, n=M_pad, axis=0)
            Uch = np.fft.fft(Uc, n=M_pad, axis=0)
            Ucrh = np.fft.fft(Ucr, n=M_pad, axis=0)
            idx_a = M - 1 + p - a_idx
            A1[p] = K_mat[p, :, None] * np.fft.ifft(Uh * z_rev_flip_h, axis=0)[idx_a]
            A3[p] = K_mat[p, :, None] * np.fft.ifft(Urh * z_rev_flip_h, axis=0)[idx_a]
            idx_b = p + b_idx
            B1[p] = np.fft.ifft(Uh * Uch, axis=0)[idx_b]
            B23[p] = np.fft.ifft(Urh * Uch + Uh * Ucrh, axis=0)[idx_b]

        pref_cx = (lam / (2.0 * np.pi)) ** 2 * (dy**2)

        def cx_lin(h):
            hf = np.fft.fft(h, n=N_fft, axis=1)
            h_rev_f = np.fft.fft(h[:, ::-1], n=N_fft, axis=1)
            s_ch1 = np.einsum("par,ar->pr", A1, hf)
            s_ch2 = np.einsum("par,ar->pr", A1, h_rev_f)
            s_ch3 = np.einsum("par,ar->pr", A3, hf)
            term_A = (
                w1 * np.fft.ifft(s_ch1, axis=1)[:, 2 * L : 2 * L + N_y]
                + w23 * np.fft.ifft(s_ch2 + s_ch3, axis=1)[:, idx_rev]
            )
            sb_ch1 = np.einsum("pbr,br->pr", B1, h_rev_f)
            sb_ch23 = np.einsum("pbr,br->pr", B23, h_rev_f)
            term_B = (
                w1 * np.fft.ifft(sb_ch1, axis=1)[:, 2 * L : 2 * L + N_y]
                + w23 * np.fft.ifft(sb_ch23, axis=1)[:, idx_rev]
            )
            return pref_cx * (2.0 * term_A.real + term_B.real)

        dim_int = M * N_y
        dim_tot = dim_int + M

        def P_even_sym(v):
            a = v[:dim_int].reshape(M, N_y)
            t = v[dim_int:]
            return np.concatenate([(0.5 * (a + a[::-1, :])).ravel(), 0.5 * (t + t[::-1])])

        def apply_full(v):
            vp = P_even_sym(v)
            d_int = vp[:dim_int].reshape(M, N_y)
            d_tail = vp[dim_int:]
            h_int = az * d_int
            rb_sum = np.sum(h_int, axis=1) * sqrt_dy + sqrt_tail * d_tail
            rb_val = (lam / np.pi) * (K_mat @ rb_sum)
            cx_val = cx_lin(h_int)
            out_int = (az * sqrt_dy) * rb_val[:, None] + az * cx_val
            out_tail = sqrt_tail * rb_val
            return P_even_sym(np.concatenate([out_int.ravel(), out_tail]))

        sym_err = 0.0
        if check_sym:
            rng = np.random.default_rng(42)
            u_r = P_even_sym(rng.standard_normal(dim_tot))
            v_r = P_even_sym(rng.standard_normal(dim_tot))
            Au = apply_full(u_r)
            Av = apply_full(v_r)
            sym_err = abs(np.dot(Au, v_r) - np.dot(u_r, Av)) / (np.linalg.norm(Au) * np.linalg.norm(v_r))
            assert sym_err < 1e-12, f"Operator asymmetry {sym_err}"

        phi_me = v_me_freq * np.sqrt(np.abs(u_me))
        v0 = P_even_sym(np.concatenate([(az * sqrt_dy * phi_me[:, None]).ravel(), sqrt_tail * phi_me]))
        full_eig, _ = eigsh(
            LinearOperator((dim_tot, dim_tot), matvec=apply_full, dtype=float),
            k=1,
            which="LA",
            v0=v0,
            tol=eig_tol,
            maxiter=50,
        )
        return {
            "x_sol": x_sol,
            "res": res,
            "nit": n_it,
            "me": me_eig,
            "full": float(full_eig[0]),
            "min_u": min_u,
            "sym_err": sym_err,
            "dt": time.time() - t0,
        }

    return run


def iv_rayleigh_lower_bound(v_rat: np.ndarray, Av_f64: np.ndarray, err_bound: float) -> str:
    """CERTIFIED (padded binary64) 50-digit Decimal lower bound on <v_rat, A_exact v_rat> / <v_rat, v_rat>."""
    with localcontext() as ctx:
        ctx.prec = 50
        ctx.rounding = ROUND_CEILING
        norm2_hi = Decimal(0)
        for val in v_rat:
            d = Decimal.from_float(float(val))
            norm2_hi += d * d
        ctx.rounding = ROUND_FLOOR
        num_lo = Decimal(0)
        for vi, avi in zip(v_rat, Av_f64):
            num_lo += Decimal.from_float(float(vi)) * Decimal.from_float(float(avi))
        ctx.rounding = ROUND_CEILING
        err_dec = Decimal.from_float(float(err_bound)) * norm2_hi
        ctx.rounding = ROUND_FLOOR
        rq_lo = (num_lo - err_dec) / norm2_hi
        return str(rq_lo)[:12]


def solve_finite_m_lattice(lam=8.0, m_omega=0.015, eta=1.0, mu0=-0.4, N_k=1026, N_w=12, certify_rq=False):
    t0 = time.time()
    dk = 2.0 * np.pi / N_k
    k = -np.pi + (np.arange(N_k) + 0.5) * dk
    xi = 2.0 * np.abs(k) / np.pi - 1.0 - mu0
    M = 2 * N_w
    n_idx = np.arange(-N_w, N_w)
    w_n = 2 * n_idx + 1
    tau = eta * C_inf * np.sqrt(lam) * m_omega
    nu = np.pi * tau

    ip_g, ia_g = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    j_arr = np.arange(-(M - 1), M)
    K_vec = m_omega**2 / (m_omega**2 + (2.0 * np.pi * tau * j_arr) ** 2)
    K_mat = K_vec[ip_g - ia_g + (M - 1)]

    M_pad = 1 << int(np.ceil(np.log2(2 * M)))
    b_idx = np.arange(M)
    idx_neg = (-np.arange(N_k)) % N_k

    def rb_op(F):
        s_k = np.mean(F, axis=1)
        return (2.0 * lam * tau) * (K_mat @ s_k)[:, None] * np.ones((1, N_k))

    def cx_op(F, G, H):
        Ff = np.fft.fft(F, axis=1)
        Hf = np.fft.fft(H, axis=1)
        G_rev_f = np.fft.fft(G[:, idx_neg], axis=1)
        pref = (2.0 * lam * tau / N_k) ** 2
        out = np.zeros((M, N_k), dtype=complex)
        for p in range(M):
            Uf = np.fft.fft(K_mat[p, :, None] * Ff, n=M_pad, axis=0)
            Vf = np.fft.fft(K_mat[p, :, None] * Hf, n=M_pad, axis=0)
            conv_ac = np.fft.ifft(Uf * Vf, axis=0)[p + b_idx]
            out[p] = np.sum(conv_ac * G_rev_f, axis=0)
        return pref * np.fft.ifft(out, axis=1)

    def P_even(v):
        a = v.reshape(M, N_k)
        a = 0.5 * (a + a[::-1, :])
        a = 0.5 * (a + a[:, ::-1])
        return a.ravel()

    sig_me = np.zeros((M, N_k), dtype=complex)
    for _ in range(100):
        z_me = 1.0 / (1j * w_n[:, None] * nu - xi[None, :] - sig_me)
        sn = rb_op(z_me)
        if np.max(np.abs(sn - sig_me)) < 1e-12:
            sig_me = sn
            break
        sig_me = 0.5 * sn + 0.5 * sig_me
    z_me = 1.0 / (1j * w_n[:, None] * nu - xi[None, :] - sig_me)
    az_me = np.abs(z_me)
    me_eig, vecs_me = eigsh(
        LinearOperator(
            (M * N_k, M * N_k),
            matvec=lambda v: P_even((az_me * rb_op(az_me * P_even(v).reshape(M, N_k))).ravel()),
            dtype=float,
        ),
        k=1,
        which="LA",
        tol=1e-7,
        maxiter=40,
    )

    def eval_F_vec(xr):
        s = (xr[: M * N_k] + 1j * xr[M * N_k :]).reshape(M, N_k)
        z = 1.0 / (1j * w_n[:, None] * nu - xi[None, :] - s)
        fn = rb_op(z) + cx_op(z, z, z)
        return np.concatenate([fn.real.ravel(), fn.imag.ravel()])

    x0 = np.concatenate([sig_me.real.ravel(), sig_me.imag.ravel()])
    x_sol, res, nit = anderson_solve(eval_F_vec, x0, tol=1e-10, maxiter=60, m=8, beta=0.4)
    sig = (x_sol[: M * N_k] + 1j * x_sol[M * N_k :]).reshape(M, N_k)
    min_u = float(np.min((w_n[:, None] * nu - sig.imag) / (w_n[:, None] * nu)))
    z = 1.0 / (1j * w_n[:, None] * nu - xi[None, :] - sig)
    az = np.abs(z)

    def apply_full(v):
        d = P_even(v).reshape(M, N_k)
        h = az * d
        return P_even((az * (rb_op(h) + 2.0 * cx_op(h, z, z).real + cx_op(z, h, np.conj(z)).real)).ravel())

    rng = np.random.default_rng(7)
    u_r = P_even(rng.standard_normal(M * N_k))
    v_r = P_even(rng.standard_normal(M * N_k))
    Au, Av = apply_full(u_r), apply_full(v_r)
    sym_err = abs(np.dot(Au, v_r) - np.dot(u_r, Av)) / (np.linalg.norm(Au) * np.linalg.norm(v_r))
    assert sym_err < 1e-12, f"Finite-m asymmetry {sym_err}"

    full_eig, vecs_full = eigsh(
        LinearOperator((M * N_k, M * N_k), matvec=apply_full, dtype=float),
        k=1,
        which="LA",
        v0=P_even(vecs_me[:, 0]),
        tol=1e-7,
        maxiter=50,
    )

    rq_lo_str = None
    if certify_rq:
        v_q = P_even(np.round(vecs_full[:, 0] * 16384.0) / 16384.0)
        Av_q = apply_full(v_q)
        u_mach = float(np.finfo(float).eps / 2.0)
        k_inf = float(np.max(np.sum(K_mat, axis=1)))
        az_max = float(np.max(az))
        op_abs = (2.0 * lam * tau) * az_max**2 * k_inf + 3.0 * (2.0 * lam * tau) ** 2 * az_max**4 * k_inf**2
        err_bound = 200.0 * M * np.sqrt(N_k) * u_mach * op_abs
        rq_lo_str = iv_rayleigh_lower_bound(v_q, Av_q, err_bound)

    q_idx = int(round(0.3 * N_k))
    chi_2kf = float(-2.0 * tau * np.mean(np.sum((z * np.roll(z, -q_idx, axis=1)).real, axis=0)))
    rpa_2kf = 2.0 * lam * chi_2kf
    return {
        "m": m_omega,
        "Nk": N_k,
        "Nw": N_w,
        "me": float(me_eig[0]),
        "full": float(full_eig[0]),
        "rq_lo": rq_lo_str,
        "min_u": min_u,
        "res": res,
        "nit": nit,
        "sym_err": sym_err,
        "rpa_2kf": rpa_2kf,
        "dt": time.time() - t0,
    }


def richardson_extrap(vals_1d: dict[tuple[int, float, int], float]) -> dict[str, float]:
    """Richardson extrapolation in 1/Ymax, dy^2, and 1/Nw^2 around (Nw=16, Ymax=240, dy=1)."""
    # 1/Ymax extrapolation: (Nw=16, Ymax=240) - (Nw=16, Ymax=120)
    # dy^2 extrapolation: (4/3)*((Nw=16, Ymax=120, dy=0.5) - (Nw=16, Ymax=120, dy=1.0))
    # 1/Nw^2 extrapolation: (9/5)*((Nw=24, Ymax=180) - (Nw=16, Ymax=180))
    base = vals_1d[(16, 240.0, 481)]
    d_ymax = vals_1d[(16, 240.0, 481)] - vals_1d[(16, 120.0, 241)]
    d_dy = (4.0 / 3.0) * (vals_1d[(16, 120.0, 481)] - vals_1d[(16, 120.0, 241)])
    f_nw24 = vals_1d[(24, 180.0, 361)]
    f_nw16 = vals_1d[(16, 180.0, 361)]
    # Extrapolate 1/Nw^2: (24^2 * f24 - 16^2 * f16)/(24^2 - 16^2) - f16 = (9/5)*(f24 - f16)
    d_nw = (9.0 / 5.0) * (f_nw24 - f_nw16)
    extrap = base + d_ymax + d_dy + d_nw
    return {"base_16_240": base, "d_ymax": d_ymax, "d_dy": d_dy, "d_nw": d_nw, "extrap": extrap}


if __name__ == "__main__":
    print("=== 1) Continuation in lam at eta=1.000 (Nw=12, Ymax=120, Ny=241, exact-sym) ===", flush=True)
    run12 = make_1d_scaling_solver(N_w=12, Ymax=120.0, N_y=241)
    for lam in [0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]:
        out = run12(lam=lam, eta=1.0)
        print(
            f"lam={lam:4.1f} | ME={out['me']:.6f} | Full={out['full']:.6f} | "
            f"min_u={out['min_u']:.4f} | sym={out['sym_err']:.1e} | res={out['res']:.1e} ({out['dt']:.1f}s)",
            flush=True,
        )
        gc.collect()

    print("\n=== 2) eta_crit root finding (exact-sym continuum) ===", flush=True)
    for lam, (e_lo, e_hi) in [(4.0, (1.005, 1.015)), (6.0, (1.040, 1.055)), (8.0, (1.058, 1.076))]:
        for nw, ym, ny in [(12, 120.0, 241), (16, 180.0, 361)]:
            solver = make_1d_scaling_solver(N_w=nw, Ymax=ym, N_y=ny)
            last_x = [None]
            def _obj(e):
                res = solver(lam=lam, eta=e, x_init=last_x[0], check_sym=False)
                last_x[0] = res["x_sol"]
                return res["full"] - 1.0
            eta_c = brentq(_obj, e_lo, e_hi, xtol=1e-5)
            rc = solver(lam=lam, eta=eta_c, x_init=last_x[0], check_sym=True)
            print(
                f"lam={lam:3.1f} Nw={nw:2d} Ymax={ym:5.1f} Ny={ny:3d} | eta_crit={eta_c:.6f} | "
                f"Full={rc['full']:.6f} | ME={rc['me']:.6f} | sym={rc['sym_err']:.1e} | res={rc['res']:.1e}",
                flush=True,
            )
            gc.collect()

    print("\n=== 3) Convergence table & Richardson extrapolation at eta=1.000 (lam=4, 6, 8) ===", flush=True)
    grid_cfgs = [
        (12, 120.0, 241),
        (12, 240.0, 481),
        (16, 120.0, 241),
        (16, 120.0, 481),
        (16, 180.0, 361),
        (16, 240.0, 481),
        (24, 180.0, 361),
    ]
    for lam in [4.0, 6.0, 8.0]:
        vals = {}
        print(f"-- lam={lam:.1f} --", flush=True)
        for nw, ym, ny in grid_cfgs:
            out = make_1d_scaling_solver(N_w=nw, Ymax=ym, N_y=ny)(lam=lam, eta=1.0)
            vals[(nw, ym, ny)] = out["full"]
            dy = 2.0 * ym / (ny - 1)
            print(
                f"  Nw={nw:2d} Ymax={ym:5.1f} Ny={ny:3d} dy={dy:.2f} | ME={out['me']:.6f} | "
                f"Full={out['full']:.6f} | min_u={out['min_u']:.4f} | sym={out['sym_err']:.1e} | res={out['res']:.1e} ({out['dt']:.1f}s)",
                flush=True,
            )
        rx = richardson_extrap(vals)
        print(
            f"  Richardson (lam={lam:.1f}): base(16,240)={rx['base_16_240']:.6f}, "
            f"d_ymax={rx['d_ymax']:+.6f}, d_dy={rx['d_dy']:+.6f}, d_nw={rx['d_nw']:+.6f} => extrap={rx['extrap']:.6f}",
            flush=True,
        )

    print("\n=== 4) Exact-symmetric finite-m periodic lattice & CERTIFIED (padded binary64) Rayleigh lower bound (mu0=-0.4) ===", flush=True)
    for lam, m_om, nk, nw, cert in [
        (6.0, 0.012, 1538, 12, False),
        (6.0, 0.008, 2050, 12, True),
        (8.0, 0.020, 1026, 12, False),
        (8.0, 0.012, 1538, 12, False),
        (8.0, 0.008, 2050, 12, False),
        (8.0, 0.005, 3074, 12, True),
    ]:
        out = solve_finite_m_lattice(lam=lam, m_omega=m_om, eta=1.0, mu0=-0.4, N_k=nk, N_w=nw, certify_rq=cert)
        rq_txt = f" | rq_lo={out['rq_lo']}" if out["rq_lo"] else ""
        print(
            f"lam={lam:.1f} m={out['m']:.4f} Nk={out['Nk']:4d} Nw={out['Nw']:2d} | ME={out['me']:.6f} | "
            f"Full={out['full']:.6f}{rq_txt} | min_u={out['min_u']:.4f} | rpa_2kf={out['rpa_2kf']:.4f} | "
            f"sym={out['sym_err']:.1e} | res={out['res']:.1e} ({out['dt']:.1f}s)",
            flush=True,
        )
        gc.collect()
