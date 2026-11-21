# RESULT.md -- final labels for the CERT/ upgrade of the VERTEX-UNIV counterexample

Scope: this directory upgrades **Theorem 1 of `../RESULT.md`** (the flat-sheet / strictly 1D vertex-corrected
linearized gap kernel crossing $\lambda_{\max}>1$ at the Migdal-Eliashberg ceiling $T_{\rm ME}$, $\eta=1$) from
`NUMERICAL` to a computer-assisted proof for a precisely stated finite model. Nothing outside `CERT/` was
modified.

---

## 1. Labels

| statement | label | evidence |
|---|---|---|
| **CERTIFIED THEOREM C.** For the finite truncated model $\mathcal{M}(N_k,N_w;\lambda,m,\eta,\mu_0)$ of `CERT.md` sec. 2 with $(N_k,N_w)=(24,12)$, $\lambda=8$, $\eta=1$, $\mu_0=-2/5$ and $m=1/100$ **and** $m=1/200$: the normal-state map $\Phi$ has a **unique** fixed point $\Sigma^\star$ with $\|\Sigma^\star-\Sigma_0\|_\infty\le10^{-12}$ (Krawczyk/Banach radii $r\le5.83329104451\cdot10^{-17}$, $k\le3.25923499198\cdot10^{-10}$; resp. $r\le3.63702630028\cdot10^{-17}$, $k\le4.50507433578\cdot10^{-10}$), and for **every** $\Sigma$ in that ball the vertex-corrected linearized gap operator satisfies $\lambda_{\max}(B(\Sigma))\ge1.00579751250$ ($m=1/100$) and $\ge1.01993371122$ ($m=1/200$). | **CERTIFIED** | `cert.py` exit 0 (`cert_main.log`); 10/10 guarded mutants exit 1 (`mutants/`); all deciding inequalities in `python-flint` ball arithmetic; second, independently written checker `independent.py` exit 0 and re-derives both inequalities for the `main` model itself at $N_k=24$ (unique fixed point, $\lambda_{\max}\ge1.0057975122$; `independent_nk24.log`) and for $N_k=8$ ($\lambda_{\max}\ge1.0388180068$; `independent_nk8.log`); and re-derives them for `m=1/200` at $N_k=24$ as well ($\lambda_{\max}\ge1.0199337108$, unique fixed point; `independent_nk24_m200.log` exit 0) |
| The same crossing on the coarser lattice $N_k=8$, $m=1/100$: $\lambda_{\max}\ge1.03881800695$ with a certified unique fixed point ($r\le9.074\cdot10^{-17}$, $k\le7.11\cdot10^{-10}$) | **CERTIFIED** | `independent.py --Nk 8` exit 0 |
| The `main` crossing at $N_k=24$ ($\lambda_{\max}\ge1.00579751243$, $r\le1.4347\cdot10^{-16}$, $k\le1.5337\cdot10^{-9}$) | **CERTIFIED** | `independent.py --Nk 24` exit 0 (`independent_nk24.log`) |
| The `m200` crossing at $N_k=24$: $\lambda_{\max}\ge1.01993371112$ ($m=1/200$) on the same certified fixed-point ball ($r\le1.00667718670\cdot10^{-16}$, $k\le2.06306136514\cdot10^{-9}$), reproduced in full by the second checker | **CERTIFIED** via both interval programs, `cert.py` ($\ge1.01993371122$, also at 400 bits) and `independent.py` ($\ge1.01993371112$) | `cert_m200.log`, `cert_prec400_m200.log`, `independent_nk24_m200_rerun.log` (`[J] PASS`, exit 0, wall 2020.8 s); `independent_nk24_m200.log` (ends at `[H]`) is the killed first attempt, kept as a record |
| Continuum $N_w\to\infty$ at fixed $(N_k,m)$: phase-blind rainbow tail $0.027$ against a $+0.012$ margin, no phase/decay bound converted into an enclosure | **NUMERICAL** (gap 2 of `../REVIEW2.md` not closed) | `../RIGOR.md`, `../verify_1d.log` |
| Continuum $N_k\to\infty$ and adiabatic $m\to0$ | **NUMERICAL** (gap 3 not closed) | `../RESULT.md` convergence tables |
| The flat-sheet/universal-curvature dichotomy as a statement about real materials | **MODEL STATEMENT** | model-level throughout; no measured $T_c$ |

## 2. Margin of the certified crossing

The certified statement is $\lambda_{\max}(B(\Sigma))\ge1+\theta$ for all $\Sigma$ in the certified
fixed-point ball, with $\theta=5\cdot10^{-3}$ *demanded*; the certified bounds exceed $1+\theta$ by

* $m=1/100$ (`main`): $\lambda_{\max}\ge1.00579751250$, i.e. $+5.80\cdot10^{-3}$ above $1$ and $+7.975\cdot10^{-4}$ above $1+\theta$;
* $m=1/200$ (`m200`): $\lambda_{\max}\ge1.01993371122$, i.e. $+1.99\cdot10^{-2}$ above $1$ and $+1.493\cdot10^{-2}$ above $1+\theta$.

`m200` is the recommended headline: larger margin, and the pure-rainbow column alone gives
$\lambda_{\rm ME}=0.994520242<1$ (**NUMERICAL**, no certified upper bound) there, so the crossing is a genuine
vertex effect (at $m=1/100$, $\lambda_{\rm ME}=0.980079149$, also **NUMERICAL**; at $m=1/400$ the rainbow alone
would already cross, which is why $m=1/400$ is kept **NUMERICAL**).

## 3. What changed relative to `../RESULT.md`

* `../RESULT.md` labeled the 1D crossing `NUMERICAL` and its Rayleigh number
  `CERTIFIED (padded binary64)`: a *point* estimate on a *numerically selected* normal branch, with a forward
  error estimate rather than an enclosure.
* Here the branch is proved to exist and to be unique (Krawczyk/Banach radii polynomial; no contraction of
  $\Phi$ is needed -- the phase-blind $\|D\Phi\|>1$ quoted in `../REVIEW2.md` gap 1 is not used anywhere), and
  the Rayleigh bound is taken over the whole certified ball, in ball arithmetic.
* A second, separately written checker (`independent.py`, different transcription, different interval
  arithmetic) reproduces both key inequalities on an independent lattice ($N_k=8$, exit 0), at $N_k=24$ for
  $m=1/100$ (exit 0) and at $N_k=24$ for $m=1/200$ (exit 0, $\lambda_{\max}\ge1.01993371112$; the first attempt at
  $m=1/200$ was killed by the laptop memory guard before its Rayleigh stage and the relaunched run completed).
  Every certified inequality of both models therefore has two independent checks, plus the 400-bit reproduction
  with `cert.py`.

* The certification is a single command: `python -B cert.py` (default `--set main`, both models) exits **0**
  (`cert_default.log`, 77 s); single-model runs, the 400-bit reproduction and the mutants are listed in
  `CERT.md` section 8.

## 4. Honest limits

The certified object is a **truncated finite model**: $M=2N_w=24$ Matsubara frequencies, $N_k=24$ momenta, with
the crossed-index restriction of the model definition. It is not the scaling limit and not the thermodynamic
limit. `../RESULT.md` keeps `NUMERICAL` for the continuum limits; what is CERTIFIED here is exactly: the
vertex-corrected linearized gap kernel of this fully specified finite model has $\lambda_{\max}>1$ at the
Migdal-Eliashberg ceiling on the unique normal branch, with validated margins $+5.80\cdot10^{-3}$ ($m=1/100$)
and $+1.99\cdot10^{-2}$ ($m=1/200$).
