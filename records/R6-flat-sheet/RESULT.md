# R6: Flat-Sheet Counterexample and Universal Curved-FS First-Crossed-Vertex Bound

**Verdict: GEOMETRIC DICHOTOMY (counterexample `NUMERICAL` converged beyond discretization error + `CERTIFIED (padded binary64)` Rayleigh bound on fixed numerical normal branch; positive theorem `proof draft` + `CERTIFIED` constants).** Bounded density of states ($A < \infty$) and a momentum-independent bare phonon interaction $V_j(\mathbf{q}) = K_j / N_0$ at fixed chemical potential are **insufficient** on their own for the adiabatic Migdal-Eliashberg ceiling $k_B T_{\mathrm{lin}} \le (1 + o(1)) C_\infty \sqrt{A M}$ as $m = \Omega/W \to 0$ ($C_\infty = 0.18272624777232952$): on 1D periodic and 2D parallel-flat-sheet bands, the first-crossed skeleton retains an $O(1)$ scaling limit that pushes the linearized pairing onset $T_{\mathrm{lin}}$ above $T_{\mathrm{ceil}} = C_\infty\sqrt{AM}$ at $\lambda \ge 4$. Conversely, under a quantitative strip-overlap decay condition (which Lemma 1 derives as a `proof draft` on $d \ge 2$ compact $C^2$ dispersions with strictly convex energy-collar surfaces, nonvanishing Fermi velocity, and non-degenerate distinct-component tangencies), the crossed terms vanish as $O(x^\alpha \log^{3+\beta}(1/x))$ as $x = \pi k_B T / W \to 0$, yielding the small-$x$ curved-FS onset bound under the stated normal-branch and symmetry hypotheses.

All results below govern the linearized pairing instability temperature $T_{\mathrm{lin}}$ on a sign-preserving normal branch of the rainbow plus first-crossed skeleton (Schrodi-Oppeneer-Aperis, arXiv:1911.12872v2, Eq. (4)), with bare phonon, fixed bare chemical potential, and no Coulomb term. None of these results is a measured $T_c$, an exact many-body theorem, or a confirmed superconductor.

---

## 1. Headline: 1D / 2D-Flat-Sheet Counterexample Above $T_{\mathrm{ceil}}$

- **Model & Scaling Limit (`proof draft`):** On the 1D periodic dispersion $\xi(k) = W(2|k|/\pi - 1) - \mu_0$ on $[-\pi, \pi)$ (or the 2D parallel-flat-sheet dispersion $\xi(k_x, k_y) = W(2|k_x|/\pi - 1) - \mu_0$ on $[-\pi, \pi)^2$) with exact flat DOS $A = 1$, $N_0 = 1/(2W)$, $|\mu_0| < W$, $4 k_F \notin 2\pi\mathbb{Z}$ (e.g. $\mu_0 = -0.4 W$), and $V_j(\mathbf{q}) = K_j/N_0$, three of the four sheet combinations satisfy linear energy conservation and converge as $m = \Omega/W \to 0$ to an $m$-independent $O(1)$ continuum first-crossed skeleton (after absorbing the $n$-independent constant $\operatorname{Re}\tilde{\Sigma}_{\mathrm{RB}}$ by $y$-translation invariance on $\mathbb{R}$). On the bare propagator $\tilde{z}_0(n, y) = (iw_n - y)^{-1}$, exact Cauchy residue integration gives non-vanishing channels $I_1$ and $I_{23}$ (`NUMERICAL` cross-check against 2D FFT convolution on $[-2000, 2000]$ matches to $2.391 \times 10^{-4}$ in `replay.py`; `CERTIFIED` flat-sheet overlap ratio $\Phi_3(\delta)/(A N_0 \delta)^2 = 9/4$).
- **Exact-Symmetric Continuum Crossing (`NUMERICAL`, `sym_err <= 2.5e-17`):** Using the exact-symmetric lumped-tail $L^2(|z|^2 dy)$ inner product in `verify_1d_counterexample.py`, the top eigenvalue `Full(eta=1.000)` of the self-adjoint linearized pairing operator crosses unity at $\lambda \ge 4$ on a smooth normal branch ($\min_{n,y} u_n/w_n > 1.10$, Anderson residual $\le 9.6 \times 10^{-11}$):

| $\lambda$ | `Full(16, 240, 481)` | `Full(24, 180, 361)` | Richardson extrap | `ME(24, 180)` | $\eta_{\mathrm{crit}}$ (`16, 180, 361`) |
|---|---|---|---|---|---|
| `4.0` | `1.003971` | `1.003043` | `1.004245` | `0.968965` | `1.008871` ($+0.89\%$) |
| `6.0` | `1.014061` | `1.013527` | `1.013678` | `0.984313` | `1.045617` ($+4.56\%$) |
| `8.0` | `1.015859` (`1.016371` at `16,120,dy=0.5`) | `1.015679` | `1.014628` | `0.990572` | `1.065128` ($+6.51\%$) |

- **Exact-Symmetric Finite-$m$ Periodic Lattice Sequence (`NUMERICAL` convergence from below; `CERTIFIED (padded binary64)` Rayleigh bound on fixed numerical normal branch):** At $\mu_0 = -0.4 W$, $\eta = 1.000$, $N_w = 12$ (`sym_err <= 5.1e-18`):
  - $\lambda = 6.0$: `Full = 1.005499` ($m=0.012, N_k=1538$) $\to$ `1.009087` ($m=0.008, N_k=2050$), with **CERTIFIED (padded binary64) Rayleigh lower bound** $\lambda_{\max}(\mathcal{K}_{\mathrm{sym}}[\Sigma_{\mathrm{num}}]) \ge 1.0090706553$ (via 50-digit directed-rounding `Decimal` dot products and Higham IEEE-754 float64 FFT/matmul forward-error bound on an explicit dyadic-rational trial vector).
  - $\lambda = 8.0$: `Full = 0.996985` ($m=0.020, N_k=1026$) $\to$ `1.005204` ($m=0.012, N_k=1538$) $\to$ `1.009412` ($m=0.008, N_k=2050$) $\to$ `1.012327` ($m=0.005, N_k=3074$), with **CERTIFIED (padded binary64) Rayleigh lower bound** $\lambda_{\max}(\mathcal{K}_{\mathrm{sym}}[\Sigma_{\mathrm{num}}]) \ge 1.0123025473$.
- **Why End-to-End Certification of Theorem 1 Remains `NUMERICAL`:** Two error terms are not closed by rigorous interval arithmetic: (1) normal-branch existence and defect $\|\Sigma_{\mathrm{exact}}^{(N_w)} - \Sigma_{\mathrm{num}}\|_\infty$ (since the phase-blind Jacobian norm $\|D\mathcal{F}_X\|_{l^\infty \to l^\infty} \approx 2.96 > 1$ defeats naive Banach contraction without inverting $I - D\mathcal{F}$ in $\sim 7.4 \times 10^4$ dimensions), and (2) the $N_w \to \infty$ Matsubara tail (phase-blind rainbow tail bound $\sim 0.027$ at $N_w=12$ exceeds the $+0.012$ eigenvalue margin, even though the signed $N_w=12 \to 24$ shift is $< 10^{-3}$) plus $N_k \to \infty, m \to 0$ discretization error.
- **CDW / Peierls Caveat (1D and 2D Parallel Flat Sheets):** In both the 1D band and the 2D parallel-flat-sheet band $\xi(k_x, k_y) = W(2|k_x|/\pi - 1) - \mu_0$, the bare static RPA charge susceptibility at nesting $\mathbf{q} = (2k_F, 0)$ satisfies `rpa_2kf = 8.41–10.37 > 1` across $\lambda = 6, 8$. Thus Theorem 1 is a skeleton-level counterexample for the pairing ceiling hypotheses, not a claim of physical superconductivity winning over CDW/Peierls order.

---

## 2. Positive Theorem on Curved Non-Nested Fermi Surfaces ($d \ge 2$)

- **Hypotheses (`RIGOR.md` §3):**
  - **(H1)** Compact $d \ge 2$ momentum space $\mathcal{B}$, inversion-even dispersion $\xi(-\mathbf{k}) = \xi(\mathbf{k}) \in [-W, W]$, bounded DOS $N(\xi) \le A N_0$ ($0 < N_0 = N(0) < \infty$).
  - **(H2)** Nonnegative bare phonon spectrum $d\lambda(\Omega)$ on $(0, \Omega_{\max}]$, $0 \le V_j(\mathbf{q}) = V_j(-\mathbf{q}) = V_{-j}(\mathbf{q}) \le v_{\max} K_j / N_0$.
  - **(H3)** Uniform three-line strip overlap decay $\sup_{\mathbf{k}\in\mathcal{B}} I_{\mathbf{k}}(u_a, u_b, u_c) \le C_{\mathrm{geom}} (A N_0)^2 u_a u_b u_c W^{-\alpha} \min(u_a, u_b, u_c)^{\alpha-1} (1 + \log_+(W/\min))^\beta$.
- **Lemma 1 (Coarea-Curvature Sufficiency, `proof draft`):** If $\xi \in C^2(\mathcal{B})$ ($d \ge 2$) has a collar $\{|\xi| \le \varepsilon_0\}$ on which $|\nabla\xi| \ge v_{\min} > 0$, $\|D^2\xi\|_\infty \le M_2$, each level set $\Sigma_\xi$ is a finite disjoint union of strictly convex $C^2$ hypersurfaces with principal curvatures $\ge \kappa_{\min} > 0$, and distinct-component parallel/antiparallel tangencies satisfy $\det(\Pi_b \mp \Pi_c) \ne 0$, then (H3) holds with $\alpha = 1, \beta = 1$ ($d=2$) and $\alpha = 1, \beta = 0$ ($d \ge 3$). *(The unproved general degeneracy-order-$p$ claim from v1 was dropped.)*
- **Lemma 2 & Theorem 2 (`proof draft` + `CERTIFIED` constants in `replay.py`):**
  1. **All-$x$ Young-ball closure ($A R \le 1/9$ under $V_j(\mathbf{q}) = K_j/N_0$ and particle-hole-even DOS $N(-\xi) = N(\xi)$):** Under these two extra hypotheses (matching the `R3` record), $\operatorname{Re}\Sigma_{\mathrm{RB}}^{(0)} \equiv 0$, so triple-Hölder contraction in $\|\Sigma\|_\infty \le \nu/3$ (`u_floor_exact = 2/3`, `c_Holder_hi < 3.648`, `young_img_hi < 0.3187 < 1/3`, `young_lip_hi < 0.9340 < 1` all `CERTIFIED` in `replay.py`) closes unconditionally for all $x > 0$ without (H3), and $A\lambda \le 1/26^2$ gives stability above $\tau_Y = \frac{117}{667}\sqrt{AM} < C_\infty\sqrt{AM}$. *(Without $N(-\xi)=N(\xi)$ and constant $V_j$, $\operatorname{Re}\Sigma_{\mathrm{RB}}^{(0)} = \Theta(\nu A v_{\max} R \log(1/x))$ leaves the ball as $x \to 0$, so the all-$x$ $A v_{\max} R \le 1/9$ claim from v1 was shrunk accordingly.)*
  2. **Adiabatic small-$x$ curved-FS window:** Assuming a sign-preserving normal branch exists with $u_n(\mathbf{k}) \ge c_0|\omega_n|$ ($c_0 \in (0, 1]$) and $|\operatorname{Re}\Sigma_n(\mathbf{k})| \le \Delta_R$ ($\Delta_R \le c_1 \nu$ or $\Delta_R \le C_{\mathrm{Re}}\nu(1 + \log_+(1/x))$), Lemma 2 (`C_LC = 3 C_geom + 34`, `24*sqrt(2) < 34` `CERTIFIED` in `replay.py`) gives $|F_X[y]| \le 10 L_X Q_0[y]$ with $L_X = O((A v_{\max} R)^2 \gamma^3 c_0^{\alpha-1} x^\alpha (1 + \log_+(1/(c_0 x)))^{3+\beta}) \to 0$ as $x \to 0$ ($\gamma = 1 + \Delta_R/(c_0\nu)$). For $V_j = K_j/N_0$, the $|Q|$-weighted Fermi-surface mean/variance split $y = \bar{y} + y_\perp$ yields $F[y] \ge (1 - \sup\sigma S(b) - 10 L_X) Q_0[\bar{y}] + (1 - 10 L_X) Q_0[y_\perp]$:
     - On an even-DOS branch within $O(\nu L_X)$ of $\operatorname{Re}\Sigma_{\mathrm{RB}} \equiv 0$, comparison with Theorem (a) of the earlier curved-surface record gives $T_{\mathrm{lin}}^{\mathrm{vertex}} \le C_\infty\sqrt{AM}(1 + C_3 B^2 \rho^\alpha(1+\log_+(1/\rho))^{3+\beta})$ in the small-$x$ window $\rho = \Omega_{\max}/W \to 0$.
     - On a general asymmetric DOS ($N(-\xi) \ne N(\xi)$) or general $V_j(\mathbf{q}) \le v_{\max} K_j / N_0$, comparison holds against Theorem (c) of the earlier curved-surface record ($\frac{1}{\pi\sqrt{2}}\sqrt{AM}$) and against $T_{\mathrm{lin}}^{\mathrm{ME}}$ on that branch.

---

## 3. What Remains Open

1. **Non-perturbative intermediate-$x$ window for $\lambda \ge 1/9$** without spherical symmetry (`Open Boundary (1)` in `RIGOR.md`).
2. **Normal-branch existence and $C_\infty\sqrt{AM}$ ceiling without particle-hole symmetry $N(-\xi)=N(\xi)$ and constant $V_j$** for $\lambda > 1/9$ (or for all-$x$ $A v_{\max} R \le 1/9$ as $x \to 0$).
3. **Full end-to-end computer-assisted `CERTIFIED` proof** of the 1D/flat-sheet counterexample (requiring a Newton-Kantorovich inversion of $I - D\mathcal{F}$ and rigorous Matsubara/continuum tail bounds) and **exact supremum** of $T_{\mathrm{lin}}^{\mathrm{vertex}}/\sqrt{AM}$ in the 1D scaling limit over $\lambda > 0$.
4. **General degenerate non-convex Fermi surfaces** (where $\det(\Pi_b \mp \Pi_c) = 0$ or principal curvatures vanish).

---

## 4. Files and Exercised Checks

- `RESULT.md`, `RIGOR.md`, `replay.py`, `verify_1d_counterexample.py` and `verify_1d.log`, the five files of this record.
- `python -B replay.py` -> **exit 0** (`CERTIFIED` constants/flat-sheet overlap + `NUMERICAL` 1D FFT residue err `2.391e-04` and 2D circle decay)
- `python -B replay.py --mutate layer_cake_factor` -> **exit 1**
- `python -B replay.py --mutate u_floor` -> **exit 1**
- `python -B replay.py --mutate residue_sign` -> **exit 1**
- `python -B replay.py --mutate no_decay` -> **exit 1**
- `python -B verify_1d_counterexample.py > verify_1d.log 2>&1` -> **exit 0** (wall time `461.58s`)
