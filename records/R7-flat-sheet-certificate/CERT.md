# CERT.md -- computer-assisted proof (CERTIFIED) for the VERTEX-UNIV 1D flat-sheet counterexample

**Status: CERTIFIED (rigorous interval / ball arithmetic, second independently written checker, guarded mutants).**
The certified object is the *finite* truncated model $\mathcal{M}(N_k,N_w;\lambda,m,\eta,\mu_0)$ defined below,
with all of its Matsubara frequencies, at $N_k=24$, $N_w=12$, $\lambda=8$, $\eta=1$, $\mu_0=-2/5$, and at the two
phonon-mass ratios $m=\Omega/W=1/100$ and $m=1/200$ (the certified set is `main` $+$ `m200`).
The continuum statements ($N_w\to\infty$, $N_k\to\infty$, $m\to0$) of `../RESULT.md` are **not** covered and stay
**NUMERICAL**; see section 7.

Revision note: this is the **v2** certificate. `cert.v1.py` / `CERT.v1.md` / `RESULT.v1.md` are kept unchanged for
audit. Their certified claims were **withdrawn**: `REVIEW3.md` found that v1's box enclosure was not a box
(`ball_of` shifted the centre by $\rho$ and kept a rounding-size radius), that a binary64 `1.0/Nk` and
`float(1.2)` leaked into the enclosures, that `arb` balls were compared with strict `<` in the max/min
reductions, and that the nonsymmetric Rayleigh perturbation inequality
$\lambda_{\max}(B)\ge\lambda_{\max}(\tfrac12(B+B^{\mathsf T}))-\tfrac12\|B-B^{\mathsf T}\|$ used there is
**false**. All five are fixed below (sections 3 and 9); the v1 numbers are not used anywhere.

Files: `cert.py` (main certificate), `independent.py` (second, separately written checker), `CERT.md`,
`RESULT.md`, logs `cert_main.log`, `cert_m200.log`, `cert_prec400.log`, `independent_nk8.log`,
`independent_nk24.log`, `independent_nk24_m200.log`, `independent_nk24_m200_rerun.log`, `ladder.log`,
`mutants/<mutant>.log`.

---

## 1. Theorem (exactly what is certified)

$\mathcal{M}$ is the finite-dimensional object of section 2: $M=2N_w=24$ Matsubara frequencies, $N_k=24$ momenta,
$D=MN_k=576$ unknowns, exact rational constants, temperature fixed at the Migdal-Eliashberg ceiling
$\tau=\eta C_\infty\sqrt\lambda\,m$. Let $\Phi:\mathbb{C}^{MN_k}\to\mathbb{C}^{MN_k}$ be its normal-state map
(rainbow $+$ first-crossed skeleton, bare phonon, fixed bare chemical potential, no Coulomb term) and
$B(\Sigma)\in\mathbb{C}^{D\times D}$ its linearised gap operator.

The input $\Sigma_0$ is the float fixed-point candidate of section 6 (deterministic: produced by the solver
recorded in `cert.py`; every entry a binary64 value that is then lifted *exactly* into the ball type).

> **Theorem C.** For each of the two parameter sets $(\lambda,m,\eta,\mu_0,N_k,N_w)=
> (8,\tfrac1{100},1,-\tfrac25,24,12)$ and $(8,\tfrac1{200},1,-\tfrac25,24,12)$:
>
> **(N) Fixed point.** $\Phi$ has **exactly one** fixed point $\Sigma^\star$ in the sup-norm ball
> $\{\Sigma:\|\Sigma-\Sigma_0\|_\infty\le\rho\}$, $\rho=10^{-12}$, with
> $\|\Sigma^\star-\Sigma_0\|_\infty\le r/(1-k)$ where
> $$m=\tfrac1{100}:\; r\le5.83329104451\cdot10^{-17},\; k\le3.25923499198\cdot10^{-10};\qquad
>   m=\tfrac1{200}:\; r\le3.63702630028\cdot10^{-17},\; k\le4.50507433578\cdot10^{-10}.$$
> The proof is a Krawczyk/Banach radii argument ($r+k\rho\le\rho$, $k<1$) with $A=(I-D\Phi(\Sigma_0))^{-1}$;
> **no** contraction of $\Phi$ itself is needed and none is asserted ($\|D\Phi\|$ is not below $1$ in the
> phase-blind norm quoted in `../REVIEW2.md` gap 1). $D\Phi$ is enclosed over the **whole ball** (each component
> an interval of radius $\ge\rho$; the enclosure is guarded, section 3).
>
> **(K) Crossing on the whole ball.** For **every** $\Sigma$ in that ball,
> $$\lambda_{\max}\big(B(\Sigma)\big)\ \ge\ 1+\theta,\quad\theta=5\cdot10^{-3},\qquad\text{with}$$
> $$\lambda_{\max}\ \ge\ 1.00579751250\ \ (m=1/100),\qquad \lambda_{\max}\ \ge\ 1.01993371122\ \ (m=1/200).$$
> Both exceed $1$ at the ceiling temperature with a validated margin
> ($\ge7.975\cdot10^{-4}$ and $\ge1.493\cdot10^{-2}$ above $1$, i.e. above $1$ by far more than the ball width).
> The bound is a rigorous Rayleigh quotient on $B$ itself: $B(\Sigma)$ is provably real symmetric for every
> $\Sigma$ (Proposition 1), so no antisymmetry penalty is needed or claimed.
>
> **(V) The vertex correction is what crosses here** -- label **NUMERICAL**, *not* certified. On the same
> lattices the pure-rainbow column of the audited implementation gives
> $\lambda_{\rm ME}=0.980079149$ ($m=1/100$) and $\lambda_{\rm ME}=0.994520242$ ($m=1/200$); these are
> binary64 `eigsh` values from `../verify_1d_counterexample.py::solve_finite_m_lattice`, i.e. a numerical
> cross-check that the rainbow kernel alone does not cross at the ceiling. No certified upper bound
> $\lambda_{\rm ME}<1$ has been computed, so vertex-decisiveness is **NUMERICAL**. (At $m=1/400$ the rainbow
> column alone would already cross, $\lambda_{\rm ME}=1.0254$; that model is a NUMERICAL diagnostic too.)

> **(S) Second checker -- both models.** The two inequalities of **(N)** are re-derived for
> **both** models by the separately written `independent.py` (section 5), and the Rayleigh inequality **(K)** is
> re-derived for $m=1/100$: $N_k=8$ `--Nk 8` exit **0** ($\lambda_{\max}\ge1.03881800695$,
> $r\le9.074\cdot10^{-17}$, $k\le7.115\cdot10^{-10}$) and $N_k=24$ `--Nk 24 --m 1/100` exit **0**
> ($\lambda_{\max}\ge1.00579751243$, $r\le1.435\cdot10^{-16}$, $k\le1.534\cdot10^{-9}$); for $m=1/200$
> (`--Nk 24 --m 1/200`) the run reproduced **(N)** ($r\le1.00667718670\cdot10^{-16}$, $k\le2.06306136514\cdot10^{-9}$,
> containment and contraction both true) but was **killed by `sm-mem-guard`** (laptop free commit fell to
> $15.6$ GB; `D:/dev/.cache/memguard/events.jsonl`, 2026-11-21T08:06:09Z, `python -B independent.py --Nk 24 --m 1/200`)
> during its Rayleigh stage, before any $q_{\rm lo}$ was printed (`independent_nk24_m200.log`, kept as the
> record of that attempt, ends at `[H]`). The relaunched command completed: `[I]` $\lambda_{\max}\ge1.01993371112$
> (enclosure width $\|B-B^{\mathsf T}\|_\infty\le6.853\cdot10^{-10}$, slack $0.0149337111253$ over $1+\theta$),
> `[J] PASS`, `exit=0`, wall $2020.8$ s (`independent_nk24_m200_rerun.log`). Both models therefore have their
> fixed point and their crossing re-derived by the second checker; for $m=1/200$ **(K)** rests on two interval
> programs, plus the 400-bit reproduction with `cert.py` (section 5).

The statement is *model-level*: no measured $T_c$, no real material, nothing about the $N_w\to\infty$ tail.

---

## 2. The model (closed form used by `cert.py`)

Strictly 1D periodic band, two flat sheets (spin-degenerate, $W=1$):

| object | definition |
|---|---|
| momenta | $k_j=-\pi+(j+\tfrac12)\tfrac{2\pi}{N_k}$, $j=0..N_k-1$ |
| dispersion | $\xi_j=2\lvert2j+1-N_k\rvert/N_k-1-\mu_0$, $\mu_0=-2/5$ (two equal-length branches, uniform DOS $A=1$, i.e. $N(\xi)=N_0=1/(2W)$ for the two branches combined) |
| Matsubara | $\omega_n=2n+1$, $n=-N_w..N_w-1$, $M=2N_w$ |
| ceiling $T$ | $\tau=\eta\,C_\infty\sqrt\lambda\,m$, $C_\infty=0.18272624777232952$, $\nu=\pi\tau$ |
| phonon kernel | $K(d)=m^2/(m^2+(2\pi\tau d)^2)$ (Einstein phonon $\omega_q=\Omega=Wm$, $\Lambda=2\pi\tau$) |
| couplings | $g_{RB}=2\lambda\tau$, $C_X=(2\lambda\tau/N_k)^2$ (exact rationals) |
| propagator | $z(n,j)=1/(\,i\omega_n\nu-\xi_j-\Sigma(n,j))$ |

$$\boxed{\;\Phi(\Sigma)(n,j)=g_{RB}\sum_a K(n-a)\,\frac1{N_k}\sum_m z(a,m)
\;+\;C_X\sum_{\substack{a,c\\0\le a+c-n<M}}K(n-a)K(n-c)\sum_{m,m'}z(a,m)\,z(c,m')\,
z(a+c-n,\;m+m'-j\;(\mathrm{mod}\,N_k))\;}$$

with $b=a+c-n$ and terms with $b\notin[0,M)$ omitted (that omission is part of the definition of the truncated
model: all $M=24$ frequencies are kept; only the "crossed" index $b$ is restricted, as in the audited
implementation). There is **one** factor $1/N_k$ (inside the rainbow mean) and **no** further $1/N_k^2$: the
$1/N_k^2$ sits inside $C_X$.

Linearised gap operator (az gauge, exactly as `verify_1d_counterexample.py::solve_finite_m_lattice`,
$\mathrm{az}(n,j)=\lvert z(n,j)\rvert$, $u=(j_0-j)\bmod N_k$, $b'=p+a_0-n$, $w=(j_0+j)\bmod N_k$):

$$B(\Sigma)_{(n,j),(a_0,j_0)}=\mathrm{az}(n,j)\,\mathrm{Re}\Big\{\mathrm{az}(a_0,j_0)\Big[
\frac{g_{RB}K(n-a_0)}{N_k}
+2C_XK(n-a_0)\!\!\sum_{\substack{p\\0\le b'<M}}\!\!K(n-p)\sum_m z(p,m)\,z(b',m+u)\Big]
+C_X\mathrm{az}(a_0,j_0)\!\!\sum_{\substack{a,c\\0\le a_0-a+n<M\\c=a_0-a+n}}\!\!K(n-a)K(n-c)\sum_m z(a,m)\overline{z(c,w-m)}\Big\}.$$

Every ingredient of section 2 is a verbatim transcription of `solve_finite_m_lattice` (see section 6 for the
numerical fidelity checks). The map is *not* projected onto inversion-even $\Sigma$ (`P_even` of the audited
code is dropped, which leaves the top eigenvalue unchanged).

### Proposition 1 ($B(\Sigma)$ is real symmetric, for every $\Sigma$)

For any $z\in\mathbb{C}^{M\times N_k}$ (no symmetry of $\Sigma$, hence of $z$, is assumed) the matrix of the
boxed formula satisfies $B(\Sigma)=B(\Sigma)^{\mathsf T}$. Write
$B_{(n,j),(a_0,j_0)}=\mathrm{az}(n,j)\,\mathrm{az}(a_0,j_0)\,[T_1+T_2+T_3]$ with

* $T_1=g_{RB}K(n-a_0)/N_k$: symmetric because $K(-d)=K(d)$ and $\mathrm{az}$ is real;
* $T_2=2C_XK(n-a_0)\sum_{p}K(n-p)\sum_m z(p,m)z(p+a_0-n,m+u)$, $u=(j_0-j)\bmod N_k$: in the transposed entry
  one has $u'=(j-j_0)=-u$ and, substituting $p'=p+a_0-n\Leftrightarrow p=p'+n-a_0$ and then $m\to m+u$,
  $$\sum_{p'}K(a_0-p')\sum_m z(p',m)z(p'+n-a_0,m-u)=\sum_{p}K(n-p)\sum_m z(p+a_0-n,m+u)z(p,m),$$
  which is $T_2$ again;
* $T_3=C_X\sum_{a+c=a_0+n}K(n-a)K(n-c)\sum_m z(a,m)\overline{z(c,w-m)}$ with $w=(j_0+j)\bmod N_k$ (symmetric in
  $(j,j_0)$): exchanging $(n,j)\leftrightarrow(a_0,j_0)$ maps $(a,c)\mapsto(c,a)$, keeps the pair set
  $a+c=n+a_0$ and the weight $K(n-a)K(n-c)$, and $m\mapsto w-m$ turns
  $\sum_m z(c,m)\overline{z(a,w-m)}$ into the complex conjugate of $T_3$'s sum, whose real part is equal.

Since $\mathrm{Re}\,\overline{s}=\mathrm{Re}\,s$, every entry satisfies $B_{(n,j),(a_0,j_0)}=B_{(a_0,j_0),(n,j)}$.
Checked numerically as a sanity check (not part of the proof): on random **non-inversion-even** perturbations of
$\Sigma_0$ ($N_k=8$ and $N_k=10$, $M=24$, $\|\delta\Sigma\|_\infty\approx10^{-3}$) the float transcription gives
$\max|B-B^{\mathsf T}|=5.55\cdot10^{-17}$ at operator scale $0.42$/$0.56$.

---

## 3. Method

**Existence/uniqueness of the normal branch (Krawczyk radii / Banach).** $\Sigma_0$ satisfies
$\|\Phi(\Sigma_0)-\Sigma_0\|_\infty\approx4\cdot10^{-17}$ in binary64. With $A=\mathrm{inv}(I-D\Phi(\Sigma_0))$
computed in binary64 and lifted exactly to $\mathbb{Q}$ (`arb` holds every double exactly, so $A$ is a fixed,
exactly known matrix), put
$$T(\Sigma)=\Sigma-A(\Phi(\Sigma)-\Sigma),\quad
r=\|A(\Phi(\Sigma_0)-\Sigma_0)\|_\infty,\quad
k=\|I-A(I-D\Phi(\mathrm{ball}))\|_\infty.$$
$r$ is enclosed from a *point* evaluation of $\Phi$ at $\Sigma_0$ (so $\rho$ does not enter $r$); $D\Phi$ is
enclosed **over the whole ball**. $r+k\rho\le\rho$ gives $T(\mathrm{ball})\subseteq\mathrm{ball}$, $k<1$ gives
contraction, hence a unique fixed point with $\|\Sigma^\star-\Sigma_0\|\le r/(1-k)$ (Krawczyk/Banach). The
matrix $A$ need not be accurate: its accuracy enters only through $r$ and $k$, both *enclosed* quantities
computed from the exactly lifted $A$.

**The box is a box (guard).** `ball_of` builds each complex component as the hull
$\big[\mathrm{Re}\,\Sigma_0-\rho,\mathrm{Re}\,\Sigma_0+\rho\big]$ of the two endpoints, and the run *checks* that
the resulting radii are $\ge\rho/2$; a mutant that certifies a point instead (`point_not_ball`) is rejected by
this guard. v1 added $\rho$ to the centre with a rounding-size radius, which is why it never enclosed $D\Phi$
over the ball (see the revision note).

**Second guard: no binary64 constants.** The certified couplings $g_{RB}$, $C_X$, $1/N_k$ and the kernel are
computed from exact integers/rationals inside the ball arithmetic, not lifted from binary64; the run checks
that their radii are $>0$ (a radius-0 "enclosure" would mean a rounded float was lifted as if it were exact).
Mutant `float_leak` trips this guard.

**The crossing (symmetric Rayleigh quotient on the ball).** By Proposition 1, $B(\Sigma)$ is real symmetric for
every admissible $\Sigma$; therefore for any real $v\neq0$
$$\lambda_{\max}(B(\Sigma))\ \ge\ \frac{v^{\mathsf T}B(\Sigma)v}{v^{\mathsf T}v}.$$
`cert.py` encloses $B$ entrywise over the ball, takes the fixed real vector $v$ (float eigenvector of the float
operator), and evaluates $q=\langle v,(\tfrac12(B+B^{\mathsf T}))v\rangle/\langle v,v\rangle$ in ball arithmetic
(the symmetrised matrix is an *exact* enclosure of $B$ by Proposition 1, so this is $q$ for $B$ itself); the
deciding inequality is $q_{\rm lo}\ge1+\theta$ on the lower endpoint. No antisymmetry penalty is subtracted --
v1's penalty term came from an inequality that is false (REVIEW3). The antisymmetric part of the *enclosure*
($\le1.6\cdot10^{-10}$ / $2.4\cdot10^{-10}$) is reported only; the run also checks it stays $\le10^{-6}$ so a
symmetric-enclosure failure would be caught.

**Precision.** Doubling the working precision to 400 bits reproduces the same certified bounds to all printed
digits ($1.00579751250$, margin $7.97512500371\cdot10^{-4}$, $r\le5.83329104451\cdot10^{-17}$,
$k\le3.25923499198\cdot10^{-10}$; $m=1/200$: $1.01993371122$, margin $1.49337112263\cdot10^{-2}$,
$r\le3.63702630028\cdot10^{-17}$, $k\le4.50507433578\cdot10^{-10}$; `cert_prec400.log` and `cert_prec400_m200.log`, both exit 0).

**Arithmetic and reporting.** All deciding inequalities are decided in `python-flint` ball arithmetic
(`arb`/`acb`, 200 bits) in `cert.py`, and with `gmpy2` directed-rounding MPFR intervals (own `RI`/`CI` classes,
validated by `--selftest` against `fractions.Fraction`) in `independent.py`. Binary64 is used **only** to
*select* $\Sigma_0$, $A$ and $v$. Every reported number is a *guaranteed* decimal: lower bounds are printed from
`bound.lower()` truncated toward $-\infty$ (`lower_str`), upper bounds from `upper_str`, never from a ball
midpoint (REVIEW3 finding 7). The decisions themselves are taken on the ball endpoints (`bool(q >= thr)`,
`bool(r + k*rho <= rho)`), not on the printed decimals.

---

## 4. Guarded mutants (`cert.py --mutate ...`, all ten exit **1**)

Each mutant is a deliberate fault injected so that the *certified* criteria must fail (all with `--models main`,
$N_k=24$). Two of them attack the machinery added in v2 (the box guard and the constant guard); the rest attack
the mathematics.

| mutant | fault | exit | where it is caught (verbatim from `mutants/<mutant>.log`) |
|---|---|---|---|
| `no_vertex` | the $C_X$ vertex term is removed from the model itself | **1** | certified Rayleigh bound $\ge0.9800791489<1+\theta$ (margin $\ge-2.49\cdot10^{-2}$), plus reference-lambda mismatch and the constant guard ($C_X$ is exactly 0, radius 0) |
| `m_threshold` | the certified model moved to $m=2/100$ (below the crossing) | **1** | Rayleigh bound $\ge0.9969667560<1+\theta$ (margin $\ge-8.03\cdot10^{-3}$) |
| `rho_zero` | verifies a point instead of a ball ($\rho=0$) | **1** | Krawczyk containment $r+k\rho\le\rho$ fails ($5.83\cdot10^{-17}\not\le0$) |
| `res_scale` | the residual enclosure is **inflated** by $2^{20}$ (so the test cannot pass by under-reporting) | **1** | Krawczyk containment: $r=6.1\cdot10^{-11}$ exceeds $\rho$ |
| `rb_scale` | rainbow kernel scaled by $6/5$ (and $\Phi$ again by $6/5$ through `box_sgn`), so $\Phi$ is inflated by $36/25$ against $D\Phi$: an inconsistent model | **1** | Krawczyk **containment** fails: $r\le8.946\cdot10^{-2}>10^{-12}$ (while $k\le0.3854<1$); the Rayleigh bound itself rises to $1.2327$ |
| `ray_thresh` | demands $\lambda_{\max}\ge1.01$ instead of $1+5\cdot10^{-3}$ | **1** | Rayleigh: $1.00579751250<1.01$ (margin $\ge-4.20\cdot10^{-3}$) |
| `diag_drop` | under-reports the diagonal of $B$ by $2\cdot10^{-2}$ | **1** | Rayleigh: $\ge0.9857975125<1+\theta$ (margin $\ge-1.92\cdot10^{-2}$) |
| `drop_lb_t3` | term 3 of the linearised operator dropped (vertex term kept in $\Phi$) | **1** | Rayleigh: $\ge0.6967262045<1+\theta$ |
| `point_not_ball` | certifies a point instead of a box: radius 0 while $\rho=10^{-12}$ | **1** | box guard: "the enclosure is not a ball of radius rho" |
| `float_leak` | lifts the certified couplings from binary64 (radius-0 balls) | **1** | constant guard: "a certified coupling was lifted from binary64" |

The per-mutant wall time (current logs, $N_k=24$, $M=24$) is $38$--$52$ s; the ten-mutant battery plus the
three certificate runs takes about $8$ min. Logs: `mutants/<mutant>.log` (verbatim output of each run). Note that `no_vertex` and
`drop_lb_t3` are caught by the Rayleigh criterion *even though* their float databases are also perturbed.

---

## 5. All bounds, with numbers

### `main`: $m=1/100$, $N_k=24$, $N_w=12$, $\lambda=8$, $\eta=1$, $D=576$

| quantity | certified value |
|---|---|
| $\|\Phi(\Sigma_0)-\Sigma_0\|_\infty$ (binary64 candidate, not certified) | $4.178\cdot10^{-17}$ |
| $r=\|A(\Phi(\Sigma_0)-\Sigma_0)\|_\infty$ | $\le5.83329104451\cdot10^{-17}$ |
| $k=\|I-A(I-D\Phi(\mathrm{ball}))\|_\infty$ | $\le3.25923499198\cdot10^{-10}$ (split part $\le3.25923153877\cdot10^{-10}$) |
| $k<1$ | holds |
| containment $r+k\rho\le\rho$ ($\rho=10^{-12}$) | $5.833\cdot10^{-17}+3.259\cdot10^{-22}\le10^{-12}$ holds |
| fixed-point radius $r/(1-k)$ | $\le5.83329104641\cdot10^{-17}$ |
| $\min\lvert\mathrm{Im}\,(i\omega_n\nu-\xi-\Sigma)\rvert$ over the ball | $\ge0.0565776818037$ (no pole) |
| $q_{\rm lo}\le\langle v,Bv\rangle/\langle v,v\rangle\le\lambda_{\max}$ | $\ge1.00579751250$ |
| antisymmetric part of the enclosure $\tfrac12\|\cdot\|_\infty$ (diagnostic) | $\le1.608\cdot10^{-10}$ |
| **certified** $\lambda_{\max}$ | $\ge1.00579751250$; margin above $1$: $\ge7.97512500371\cdot10^{-4}$ |

### `m200`: $m=1/200$, $N_k=24$, $N_w=12$, $\lambda=8$, $\eta=1$, $D=576$

| quantity | certified value |
|---|---|
| $\|\Phi(\Sigma_0)-\Sigma_0\|_\infty$ (binary64 candidate) | $2.861\cdot10^{-17}$ |
| $r$ | $\le3.63702630028\cdot10^{-17}$ |
| $k$ | $\le4.50507433578\cdot10^{-10}$ |
| containment | $3.637\cdot10^{-17}+4.505\cdot10^{-22}\le10^{-12}$ holds |
| fixed-point radius | $\le3.63702630192\cdot10^{-17}$ |
| $\min\lvert\mathrm{Im}\,(i\omega_n\nu-\xi-\Sigma)\rvert$ | $\ge0.0384152648904$ |
| $q_{\rm lo}$ | $\ge1.01993371122$ |
| antisymmetric part of the enclosure (diagnostic) | $\le2.443\cdot10^{-10}$ |
| **certified** $\lambda_{\max}$ | $\ge1.01993371122$; margin above $1$: $\ge1.49337112263\cdot10^{-2}$ |

### Second checker, `independent.py` (separately written, different interval arithmetic)

`independent.py` re-derives $\Phi$, $\partial\Phi$ and $B$ from the model definition above with its own
transcription (different contraction order, different table layout, no code shared with `cert.py`), its own
arithmetic (own `gmpy2` MPFR directed-rounding interval classes at 160 bits, validated by `--selftest` against
`fractions.Fraction`), its own box construction (endpoints, not centres) and its own Krawczyk construction
($E=I-D\Phi$ assembled explicitly, $A$ from its own float solve). Two things are *not* independent and are
declared: the matrix products in its Krawczyk step use `python-flint`'s `acb_mat` (the same library `cert.py`
uses; the inequalities they feed are independent), and both checkers share the model definition of section 2.

| quantity | $N_k=8$, $m=1/100$ | $N_k=24$, $m=1/100$ | $N_k=24$, $m=1/200$ | `cert.py` $m=1/100$ | `cert.py` $m=1/200$ |
|---|---|---|---|---|---|
| $\|\Phi(\Sigma_0)-\Sigma_0\|_\infty$ (float) | $5.594\cdot10^{-17}$ | $7.076\cdot10^{-17}$ | $2.964\cdot10^{-17}$ | $4.178\cdot10^{-17}$ | $2.861\cdot10^{-17}$ |
| $\partial\Phi$ vs central differences | $3.04\cdot10^{-10}$ | $4.19\cdot10^{-10}$ | $1.80\cdot10^{-10}$ | $5.53\cdot10^{-10}$ | $2.56\cdot10^{-10}$ |
| `solve_finite_m_lattice` reference $\lambda_{\rm full}$ | $1.038818007$ (delta $3.6\cdot10^{-11}$) | $1.005797512$ (delta $2.2\cdot10^{-10}$) | $1.019933712$ (delta $4.0\cdot10^{-10}$) | $1.005797512$ | $1.019933712$ |
| $\min|\mathrm{Im}\,u|$ over the ball | $\ge9.782252\cdot10^{-2}$ | $\ge5.657768\cdot10^{-2}$ | $\ge3.841526\cdot10^{-2}$ | $\ge5.6578\cdot10^{-2}$ | $\ge3.8415\cdot10^{-2}$ |
| $r$ | $\le9.07419094177\cdot10^{-17}$ | $\le1.43466873225\cdot10^{-16}$ | $\le1.00667718670\cdot10^{-16}$ | $\le5.833\cdot10^{-17}$ | $\le3.637\cdot10^{-17}$ |
| $k$ | $\le7.11496004021\cdot10^{-10}$ | $\le1.53363171691\cdot10^{-9}$ | $\le2.06306136514\cdot10^{-9}$ | $\le3.259\cdot10^{-10}$ | $\le4.505\cdot10^{-10}$ |
| $r+k\rho\le\rho$, $k<1$ | both true | both true | both true | both true | both true |
| $q_{\rm lo}$ | $\ge1.03881800695$ | $\ge1.00579751243$ | $\ge1.01993371112$ | $\ge1.00579751250$ | $\ge1.01993371122$ |
| antisymmetric part of enclosure | $\le2.161\cdot10^{-10}$ | $\le5.033\cdot10^{-10}$ | $\le6.853\cdot10^{-10}$ | $\le1.6\cdot10^{-10}$ | $\le2.4\cdot10^{-10}$ |
| **certified** $\lambda_{\max}\ge$ | $1.03881800695$ | $1.00579751243$ | $1.01993371112$ | $1.00579751250$ | $1.01993371122$ |
| slack above $1+\theta$ | $3.4\cdot10^{-2}$ | $7.98\cdot10^{-4}$ | $1.49\cdot10^{-2}$ | $7.98\cdot10^{-4}$ | $1.49\cdot10^{-2}$ |
| verdict / wall | PASS / 250 s | PASS / 2551 s | PASS / 2020.8 s (relaunch; the first attempt was killed at `[H]`, 2106 s) | CERTIFIED / $\approx40$ s | CERTIFIED / $\approx40$ s |

The two $N_k=24$ columns are the decisive ones: the second checker reproduces the certified crossing of the
$m=1/100$ model and of the $m=1/200$ model itself (the $m=1/200$ column is the relaunched run of section 8 row
17b; see **(S)**). Its enclosure is wider than `cert.py`'s (entrywise $\mathrm{acb}$ at 160 bits, own assembly),
which is why its $N_k=24$, $m=1/100$ bound reads $1.00579751243$ against `cert.py`'s $1.00579751250$: both are
rigorous lower bounds for the same $\lambda_{\max}$, both clear $1+\theta$ by $>7.9\cdot10^{-4}$, and the headline
bound in section 1 quotes the larger (`cert.py`) value. `independent.py` certifies the *inequalities* only; it does not certify the
vertex-decisiveness sub-claim **(V)**, which stays NUMERICAL.

## 6. Model fidelity (numerical stage; not part of the certified decision)

| model | $\Phi$ vs FFT map | $B$ vs `apply_full` | $\partial\Phi$ vs central differences | float $\|B-B^{\mathsf T}\|_\infty$ | reference `solve_finite_m_lattice` | pure-rainbow column $\lambda_{\rm ME}$ |
|---|---|---|---|---|---|---|
| `main` ($m=1/100$) | $8.89\cdot10^{-17}$ | $6.66\cdot10^{-16}$ | $5.53\cdot10^{-10}$ | $1.67\cdot10^{-16}$ | $1.005797512$ (delta $2.2\cdot10^{-10}$) | $0.980079149$ |
| `m200` ($m=1/200$) | $7.51\cdot10^{-17}$ | $4.44\cdot10^{-16}$ | $2.56\cdot10^{-10}$ | $1.11\cdot10^{-16}$ | $1.019933712$ (delta $4.0\cdot10^{-10}$) | $0.994520242$ |

This is the evidence that the certified object *is* the object of `../RESULT.md`; it is numerical evidence and
carries the NUMERICAL label. The last column is the NUMERICAL basis of **(V)** ($\lambda_{\rm ME}<1$ while the
certified full kernel is $>1$), not a certified statement.

---

## 7. What is NOT certified

| item | label | why |
|---|---|---|
| unique normal fixed point and $\lambda_{\max}>1$ for $\mathcal{M}(24,12;m=1/100)$ and $\mathcal{M}(24,12;m=1/200)$ | **CERTIFIED** | sections 1, 3, 5; two checkers; ten mutants |
| vertex-decisiveness: pure-rainbow $\lambda_{\rm ME}<1$ on the same lattice | **NUMERICAL** | binary64 `eigsh` only; no certified upper bound |
| Matsubara tail $N_w\to\infty$ at fixed $(N_k,m)$ | **NUMERICAL** | the phase-blind rainbow tail is $0.027>0.012$ margin and no phase/decay argument has been converted into an enclosure (`../RIGOR.md`); the certified statement uses exactly the $M=24$ frequencies of the model |
| $N_k\to\infty$ (continuum momentum), $m\to0$ (adiabatic limit) | **NUMERICAL** | `../RESULT.md` convergence tables; certification is at $N_k=24$ |
| the `--set ladder` models | **NUMERICAL** | printed as diagnostics only (`ladder.log`) |
| "the vertex correction lifts the ME ceiling for a real material" | **MODEL STATEMENT** | model-level; no measured $T_c$ |

`cert.py --set ladder` (float stage only, no interval certificate):

| model | $N_k$ | $N_w$ | $\lambda_{\max}$ (float) | reading |
|---|---|---|---|---|
| `m400` ($m=1/400$) | 24 | 12 | $1.0681283523$ | crosses; but $\lambda_{\rm ME}=1.0254>1$, so it does not isolate the vertex effect |
| `nk48` | 48 | 12 | $1.0074033857$ | the crossing survives $N_k=24\to48$ |
| `nk64` | 64 | 12 | $1.0073349348$ | and $N_k=64$ (converging, still $>1$) |
| `nw16` | 24 | 16 | $1.0054800598$ | survives adding four Matsubara frequencies (no tail bound) |
| `m20` ($m=2/100$) | 24 | 12 | $0.9969667561$ | below the crossing, as expected (same point as the `m_threshold` mutant) |
| `lam6` ($\lambda=6$) | 24 | 12 | $1.0061798062$ | crosses at $\lambda=6$ as well |

---

## 8. Runtime, commands and exit codes

Python 3.13 via `uv run --no-project --with mpmath --with numpy --with scipy --with python-flint --with gmpy2`;
single process, $<2$ GB RSS, laptop CPU only, temp files under `C:/t/VertexCert` (deleted afterwards).
`UV=uv run --no-project --with mpmath --with numpy --with scipy --with python-flint --with gmpy2`,
cwd = this directory:

| # | command | exit | wall |
|---|---|---|---|
| 0 | `$UV python -B cert.py` (default `--set main` = `main,m200`, one process) | **0** | 77 s |
| 1 | `$UV python -B cert.py --models main` | **0** | 47 s |
| 2 | `$UV python -B cert.py --models m200` | **0** | 40 s |
| 3 | `$UV python -B cert.py --models main --prec 400` | **0** | 37 s |
| 3b | `$UV python -B cert.py --models m200 --prec 400` | **0** | 44 s |
| 4 | `$UV python -B cert.py --models main --mutate no_vertex` | **1** | 41 s |
| 5 | `$UV python -B cert.py --models main --mutate m_threshold` | **1** | 39 s |
| 6 | `$UV python -B cert.py --models main --mutate rho_zero` | **1** | 40 s |
| 7 | `$UV python -B cert.py --models main --mutate res_scale` | **1** | 43 s |
| 8 | `$UV python -B cert.py --models main --mutate rb_scale` | **1** | 42 s |
| 9 | `$UV python -B cert.py --models main --mutate ray_thresh` | **1** | 38 s |
| 10 | `$UV python -B cert.py --models main --mutate diag_drop` | **1** | 38 s |
| 11 | `$UV python -B cert.py --models main --mutate drop_lb_t3` | **1** | 46 s |
| 12 | `$UV python -B cert.py --models main --mutate point_not_ball` | **1** | 47 s |
| 13 | `$UV python -B cert.py --models main --mutate float_leak` | **1** | 52 s |
| 14 | `$UV python -B independent.py --selftest` | **0** | 2 s (4000 random interval trials vs `fractions.Fraction`: 0 containment failures, 0 width failures) |
| 14b | `$UV python -B independent.py --selftest-gap` | **0** | 6 s (non-even $\Sigma$: `gap_tab` vs the raw 4-index definition, max err $9.5\cdot10^{-18}$, rel $9.3\cdot10^{-16}$) |
| 15 | `$UV python -B independent.py --Nk 8` | **0** | 250 s (PASS) |
| 16 | `$UV python -B independent.py --Nk 24 --m 1/100` | **0** | 2551 s (PASS) |
| 17 | `$UV python -B independent.py --Nk 24 --m 1/200` | **1** (killed by `sm-mem-guard`, not by a failed inequality: `independent_nk24_m200.log` ends at `[H]` with containment and contraction both true) | 2106 s |
| 17b | `$UV python -B independent.py --Nk 24 --m 1/200` (relaunch) | **0** (`[J] PASS`, `lambda_max >= 1.01993371112`, all certified inequalities hold) | 2020.8 s |
| 18 | `$UV python -B cert.py --set ladder` (diagnostics, NUMERICAL) | **0** | 680 s |
| 19 | `$UV python -B independent.py --selftest` (rerun after the REVIEW4 edits) | **0** | 0.6 s (4000 trials, 0 failures) |
| 19b | `$UV python -B independent.py --selftest-gap` (rerun after the REVIEW4 edits) | **0** | 1.6 s (`PASS`) |
| 20 | `$UV python -B C:/t/VertexCert/micro4.py` (REVIEW4 micro check: `r/(1-k)` upper bound and `q_lo` lower bound, incl. a negative numerator and $[-1,1]/[2,3]$, against exact `mpq`) | **0** | 1 s |

Logs: `cert_main.log` (1), `cert_m200.log` (2), `cert_prec400.log` (3), `cert_prec400_m200.log` (3b), `mutants/<mutant>.log` (4--13),
`cert_default.log` (0), `independent_nk8.log` (15), `independent_nk24.log` (16), `independent_nk24_m200.log` (17),
`independent_nk24_m200_rerun.log` (17b), `ladder.log` (18).
Row 0 was run after the REVIEW4 fixes; it prints the multi-model summary
`main ... lambda_max >= 1.00579751250 CERTIFIED` and `m200 ... lambda_max >= 1.01993371122 CERTIFIED` and exits **0**
(before the fix the same command exited 1 in the summary printer, `REVIEW4.md` P1-1).

---

## 9. Changes from v1 (resolution of `REVIEW3.md`)

| REVIEW3 finding | fix |
|---|---|
| **F1 (P0)** `ball_of` returned a shifted point with a rounding-size radius, so $D\Phi$/$B$ were never enclosed over the ball | `ball_of` now builds each component as the hull of the two endpoints $s_0\pm\rho$; new **box guard** checks that every component radius is $\ge\rho/2$; new `point_not_ball` mutant is rejected by the guard; $\|B-B^{\mathsf T}\|$ over the box rose from $2.3\cdot10^{-58}$ (an artefact of the point evaluation) to $1.6\cdot10^{-10}$/$2.4\cdot10^{-10}$, as expected for a radius-$10^{-12}$ box |
| **F2 (P0)** binary64 `1.0/Nk` and `float(1.2)` leaked into the `acb` enclosures | exact `arb(1)/arb(Nk)` and rational `fac=(6,5)` in the certified stage; new **constant guard** (radii must be $>0$) and new `float_leak` mutant |
| **F3 (P0)** strict `<`/`>` between overlapping `arb` balls in five max/min reductions | all replaced by `arb.max()`/`arb.min()` (`rowsum_inf`, `r_hi`, `delta_hi`, `min_imu`; `acc.b > delta` in `independent.py` was already endpoint-exact) |
| **F4 (P0)** false nonsymmetric Rayleigh perturbation inequality | removed everywhere; replaced by **Proposition 1** ($B(\Sigma)=B(\Sigma)^{\mathsf T}$ identically, section 2) so the plain symmetric Rayleigh quotient is rigorous; the `antisym_inject` mutant (which tested the bogus penalty) is dropped and replaced by `diag_drop`, `point_not_ball`, `float_leak` |
| **F5 (P1)** momentum-shift sign bug in `independent.py::gap_tab` | fixed to `(j - j0) % Nk`; `--selftest-gap` now checks $B$ against a raw 4-index assembly on a **non-even** $\Sigma$ |
| **F6 (P1)** CERT.md model formulas did not match the code (double $1/N_k^2$, missing momentum sums, $A=2$) | section 2 rewritten from the code, term by term; $A=1$ |
| **F7 (P1)** midpoints printed as certified bounds | `lower_str`/`upper_str` print guaranteed decimals from the ball endpoints ($\pm8$ ulp margin); decisions are taken on the endpoints, not on the decimals |
| **F8 (P1)** mutant descriptions and backend descriptions wrong | corrected (`res_scale` inflates; `rb_scale` is caught by the Krawczyk criterion, not by $k<1$; `drop_lb_t3` drops term 3 only; `independent.py` uses its own `gmpy2` intervals plus a shared `acb_mat` product) |
| **F9 (P2)** `RI.exact` rounded integers in $(2^{53},2^{60})$ to zero width; `delta/2` computed under round-down | threshold lowered to $2^{53}$; the halving is gone with the penalty term |
| **Check 5/6** **(V)** was presented without a NUMERICAL label in `RESULT.md` | **(V)** is now explicitly `NUMERICAL` in section 1 and in `RESULT.md` |


## 10. Changes from v2 (resolution of `REVIEW4.md`; backups kept as `*.v2.*`)

| REVIEW4 finding | fix |
|---|---|
| **P1-1** default `python -B cert.py` crashed (`TypeError`) in the multi-model summary because the certified bound is stored as a decimal string | the summary now prints the string bound verbatim and formats a float bound only when it is not a string (`cert.py`, summary loop); the default invocation exits **0** (section 8 row 0) |
| **P2-1** stale v1 docstrings | `cert.py` header: the exact boxed $B(\Sigma)$ formula of section 2 (with $c=a_0-a+n$) replaces the v1 three-term form, and **(K)** now states the Proposition-1 symmetric Rayleigh quotient with no antisymmetry penalty; the `rb_scale` mutant comment states the actual $6/5$ and $36/25$ scalings; `independent.py` header: the stale "matrix-free $\| |A|W \|_\infty$" bullet is replaced by the actual signed `acb_mat` product |
| **P2-2** documentation inaccuracies | section 2 boxed formula now carries $c=a_0-a+n$ in the third sum; the `rb_scale` row reports the actual failure ($r\le8.946\cdot10^{-2}>10^{-12}$, containment; $k\le0.3854<1$); `RESULT.md` §2 marks $\lambda_{\rm ME}$ values **NUMERICAL** inline |
| **P2-3** `independent.py` printing/sign nits | the $r/(1-k)$ print uses the directed interval quotient $(r_\mathrm{iv}/(1-k_\mathrm{iv}))_\mathrm{up}$ instead of a binary64 division; `q_lo = q.re.a` (the directed interval quotient, correct for every sign) instead of `num.re.a/den.b`; both verified by the micro check of row 20 |
| **P1-2** `m200` second-checker Rayleigh stage incomplete | **closed after this revision**: the relaunched `independent.py --Nk 24 --m 1/200` finished with `[J] PASS`, `exit=0`, `lambda_max >= 1.01993371112`, wall `2020.8` s (`independent_nk24_m200_rerun.log`, section 8 row 17b); **(S)** in section 1 now records both models (`--Nk 24 --m 1/100` and `--Nk 24 --m 1/200`), and the killed `independent_nk24_m200.log` is kept as the record of the earlier attempt |
