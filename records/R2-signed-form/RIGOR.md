# First-crossed momentum form (infinite volume and Matsubara line)

## 1. Exact signed identity

Use energy units `τ=k_B T`, `ν=πτ`, `ω_n=(2n+1)ν`. A four-momentum `p=(k,n)` has fermionic index; an increment `s=(q,j)` has bosonic index. Integrals include counting over Matsubara indices and translation-invariant momentum measure. The positive interaction `V_s=V_{-s}` has units of energy for normalized momentum measure. Assume separately even frequency and momentum, as in the models below. Write `a=p-s`, `b=p-s-t`, `c=p-t`. Use the bare phonon, fixed chemical potential, no Coulomb term, dressed normal line, and precisely the rainbow and the first crossed skeleton of [Schrodi–Oppeneer–Aperis, Eq. (4)](https://arxiv.org/html/1911.12872v2#S2.E4). We do not sum vertex ladders or dress the phonon.

Write `C∞=√g(2)/(2π)` for the exact scalar spectral coefficient, not its rounded decimal; `g(2)≥1` by the first-coordinate soft trial, hence `C∞>1/8` since `π<4`. The sharper certified enclosure is in `../GAMMA/RESULT.md`. In this note `τ_lin` is the upper instability boundary on the specified regular normal branch: the supremum of temperatures with `sup σ(B)>1` (zero if none). A stability exclusion does not require the existence or uniqueness of a crossing.

Put `z_p=(iω_p-ξ_p-Σ_p)^(-1)=r_p-iQ_p` and `u_p=ω_p-Im Σ_p`, so `Q_p=u_p|z_p|²`. On a sign-preserving normal branch `ω_p Q_p>0`. With real even anomalous self-energy `φ`, let `y_p=φ_p/u_p` and `h_p=φ_p|z_p|²=Q_p y_p`. To first order, the Nambu line is `G_p=[[z_p,-h_p],[-h_p,-conj(z_p)]]`. The normal electron equation and real anomalous equation are

```
Σ_p = τ ∫ V_s z_a + τ² ∫ V_s V_t z_a z_b z_c,
φ_p = τ ∫ V_s h_a + τ² ∫ V_s V_t [Re(z_b z_c)h_a
                                  + Re(z_a conj(z_c))h_b + Re(z_a z_b)h_c].
```

For example, `σ3 G_a σ3 G_b σ3 G_c σ3` has linear 12 entry `h_a conj(z_b z_c)+z_a h_b conj(z_c)+z_a z_b h_c`. Interchanging `s,t` makes the integrated coefficient real for real `h`. Taking minus the imaginary part of the normal equation, multiplying the resulting `u_p` equation by `y_p`, subtracting the pairing equation, and testing with `Q_p y_p` gives `F=Q0+F_ME+F_X`, where

```
Q0   = ∫ ω_p Q_p |y_p|²,
F_ME = (τ/2) ∫ V_s Q_p Q_a |y_p-y_a|²,
F_X  = τ² ∫ V_s V_t [ -(Q_p Q_a Q_b Q_c/4) |y_p-y_a+y_b-y_c|²
                      +(1/2) Σ_{v=a,b,c} Q_p Q_v r_{v'} r_{v''} |y_p-y_v|² ].   (1)
```

`v',v''` are the other two among `a,b,c`. In the normal equation `-Im(z_a z_b z_c)=Q_a r_b r_c+r_a Q_b r_c+r_a r_b Q_c-Q_a Q_b Q_c`. The product of four `Q` multiplies `-y_p(y_p-y_a+y_b-y_c)` before reanchoring. Average over the four anchors using `(p,s,t)`, `(a,-s,t)`, `(b,-s,-t)`, `(c,s,-t)`. Reanchoring to `v` for each two-real-line term gives the factor 1/2. The same proof polarizes to the Hermitian form on complex trial functions. Start with bounded trials with finite frequency support; the applications below provide summable form bounds, and hence extension by density. In particular, retaining only the alternating square is wrong, and neither the square coefficient nor the pair coefficients have a universal sign. The alternating square vanishes for zero **four**-transfer, not merely for a static phonon or for zero momentum transfer at nonzero bosonic frequency.

The same elimination gives `F[y]=<d,(I−B)d>` for the self-adjoint raw pairing kernel `B` in `d_p=|z_p|u_p y_p`. The underlying squared norm is `∫ u_p Q_p |y_p|²`. All normal branches used here have `|u_p/ω_p|` bounded above and away from zero (by the sup-norm branch bounds relative to `ν`), so a positive gap relative to `Q0` implies a gap in the pairing spectrum.

## 2. Forward model: no onset above the ceiling in the ultranarrow limit

On a 1D zone `[-π,π)` with measure `dk/(2π)`, use `ξ=W(2|k|/π-1)`. It is affine near the two Fermi points `±π/2`, with speed `v_F=2W/π`; the corners at band edges are not stationary Fermi points. Pushforward to energy is exactly `N(ξ)=1/(2W)` on `(-W,W)` (per spin), so `A=1` and the bare band-bottom `E_F=W`. Let `V_j(q)=4λW 1_{|q|≤q_c} Ω²/[Ω²+(2νj)²]`, with `0<q_c<π`. The conventional double Fermi-surface average `N(0) E_FS[V_0(k-k')]` is `λ`: of the four equal-weight ordered point pairs, the two with `q=0` are inside the window. Thus the Einstein moment in the ceiling is `M=λΩ²`, not a mass-enhancement substitution. The integrated zero-frequency weight is `B0=∫V_0=2λ v_F q_c`.

Let `B=Σ_j ∫ V_j`, `r=B/(π²τ)=τ B/ν²`. Positivity, monotone integral comparison, and the Einstein kernel give `B/B0≤1+Ω/(2τ)`; the separate tail estimate `B/B0-1≤(Ω/ν)²` follows by `Σ_{j≥1} j^(-2)≤2`. Both are independent of Matsubara truncation. In the `L∞` ball `||Σ||≤ν/2`, `|z|≤2/ν` and `|z(Σ)-z(Σ')|≤4 ||Σ-Σ'||/ν²`. The normal map has image bound `||Σ_new||/ν≤2r+8r²` and Lipschitz bound `4r+48r²` by telescoping the product. At `r≤1/16` these are at most `5/32` and `7/16`; hence the ball has a unique fixed point by contraction. The solution continues from zero interaction and preserves inversion and frequency-conjugation symmetries by uniqueness, with `sgn u=sgn ω`. In the symmetric pairing variable `d_p=|z_p|φ_p`, each pairing matrix entry has two external absolute lines and each crossed term two more. Its absolute row sums are bounded by `4r+48r²`; by reversal of the rainbow or reanchoring each of the three crossed terms, the kernel is symmetric (each real product is invariant), so the same bound applies to columns and thus to `L²` by Schur. The norm is strictly below 1. This is a bound on this normal branch, not on disconnected solutions of the truncated nonlinear equations.

At fixed finite `λ>0`, if `v_F q_c/Ω→0`, the bound on `r` tends uniformly to zero for `τ≥c sqrt(λ) Ω` at any fixed `c>0`. Thus the ceiling, and in fact any such positive fraction of it, cannot be attained on this branch as `m=Ω/W→0`. The perfect-forward model in [Rademaker et al.](https://arxiv.org/html/1507.03967v2) holds fixed the integrated `g₀²/Ω²` (the leading-order perfect-forward mass-renormalization parameter), giving `B0/Ω=2g₀²/Ω²`. In the 1D affine window, put `s=v_F q_c/Ω`. With `V_0(q)=2|g_q|²/Ω`, the integrated strength `∫|g_q|² dq/(2πΩ²)=λs`, whereas the bare zero-temperature one-loop frequency derivative at a Fermi point is `∫|g_q|²/(Ω+|ξ(k_F−q)|)² dq/(2π)=λs/(1+s)` by direct integration across the window. Keeping the integrated strength nonzero as `s→0` requires `λ` to diverge. Neither parameter here is a fully dressed mass or velocity enhancement. An exact transfer delta is not an ordinary function on the Fermi-surface product measure.

## 3. A signed-form obstruction that does not vanish as m → 0

The same bounded-DOS model with width proportional to `Ω` has a negative crossed form of order *independent of `m`* relative to `Q0`. Precisely, fix `λ≥1` and `t=τ/Ω≥2^12`, let `δ=2^-16`, `ε=δ/λ`, and choose `v_F q_c=εν`. As `m→0`, `q_c<π/16` and `λΩ/E_F→0`. In this family `r0=τ B0/ν²=2δ/π`, `s*=B/B0-1≤2^-24`, `r=r0(1+s*)≤2δ`; thus §2 gives the normal branch and `||Σ||/ν≤2r+8r²≤5δ`. For the witness exactly at the ceiling, take `λ=2^32` and `t=C∞·2^16>2^12` with the exact spectral `C∞` from §1. This is a large but fixed formal coupling; small `λΩ/E_F` alone does not control omitted diagrams or phonon renormalization.

Take `y(p)=sgn(ω_n)` for `n=0,-1` and momenta within `q_c/16` of either Fermi point, and zero elsewhere. The support has finite positive momentum-times-counting measure `H`; `φ=u y` is real, inversion-even and frequency-even. For double-static transfers (`j_s=j_t=0`), if a corner belongs to this support then all four are within `(2+1/16)q_c` of its Fermi point and have `n=0` or `-1`. The affine dispersion and the branch bound imply `|e|/ν≤8δ` and `1-5δ≤|u|/ν≤1+5δ`, where `e=ξ+Re Σ`. Hence `1/(2ν)≤|Q|≤2/ν`, `|r_p|≤16δ |Q_p|` there. Globally `|Q|,|r_p|≤2/ν`. Also `Q0≤2H`.

The alternating double-static term has nonpositive sign everywhere. Restrict it to `p` in the support and to both transferred momenta in `[q_c/2,3q_c/4]` (positive, same sign); the other three corners are then outside the support. Each transfer has exactly `1/8` of the static weight, the squared alternating difference is 1, and the four lines have product at least `1/(16ν⁴)`. Thus the double-static alternating term is at most `-τ² B0² H/(4096 ν⁴)`. Double-static real-line terms are bounded above in absolute value by `24576 δ² τ² B0² H/ν⁴`: use `|y_p-y_v|²≤2(1_support(p)+1_support(v))`, reanchor the support, sum over three `v`, and bound each line product by `16(16δ)²/ν⁴` on the relevant sets.

For all transfer pairs with at least one nonstatic transfer, use `|y_p−y_a+y_b−y_c|²≤4Σ_{w=p,a,b,c}1_support(w)` (counting all four corner occurrences even when some coincide), the preceding pair bound, four global line bounds, and the same reanchorings. The result is `|F_X,nonstatic|≤160 τ²(B²−B0²) H/ν⁴≤480·2^-24 τ² B0² H/ν⁴`. The exact localized-witness margin is `1/4096−24576·2^-32−480·2^-24=55/262144>1/8192`. Therefore

```
F_X[y]/Q0[y] ≤ -r0²/16384 < 0,    independently of Ω/W.              (2)
```

This no-go concerns `F_X≥−o(1)Q0`, not the onset ceiling. In both witness families the pairing operator norm is at most `7/16` and `u_p/ω_p≥1/2`, so `F≥(9/16)∫u_p Q_p |y_p|²≥(9/32)Q0` for all admissible trials. Conversely, the absolute reanchoring estimate in (4) with all internal lines bounded by `2/ν` gives `|F_X|≤80 r² Q0`, so the witnessed order in `m` is constant, without a claim of an optimal coefficient. The reference `A` and the global adiabatic ratio alone are not a joint phase-space condition.

### 3b. The no-go persists on the momentum-independent gap sector

The angular variance cannot by itself pay for all of `F_X`: a separate test has zero angular variance. Keep the model and fixed `λ≥1`, `t≥2^12` from §3, but permit any sufficiently small fixed `δ>0` with `ε=δ/λ`. Take `y_n=sgn(ω_n)` only for `n=0,-1`, zero otherwise, independent of momentum. (The gap `Δ=ω y` is momentum-independent, and the anomalous self-energy `φ=u y` has the required parities.) Its `Q0` is positive and bounded by `16π N0 ν`: on the two active frequencies use `|u|≤2ν`, `|z|≤2/(ν²+ξ²)^(1/2)` and integrate the flat DOS. The normal branch still has `||Σ||/ν≤5δ`.

Write `k_j=V_j(q)/V_0(q)=1/[1+(2πtj)²]` on the interaction support. Divide `F_X` by `r0² N0 ν`, where `r0=τ B0/ν²`. To compute the limiting normalized form, write `x=ξ_p/ν` and integrate each transfer with its uniform probability measure on `[-q_c,q_c]`. Replace the dimensionless dressed corner lines by bare external-energy lines `z^0_n(x)=1/[i(2n+1)−x]`. The bare energy change at any corner is `|ξ_{k'}−ξ_k|/ν≤2ε` by the globally Lipschitz dispersion; adding the dressed shift gives `η=2ε+5δ`. With `D(x)=(1+x²)^(-1/2)` each replacement changes a line by at most `2η D(x)²`, and all lines are bounded by `2D(x)` for small `δ`. Telescoping four factors costs at most `64η D(x)^5`. Bound each square by a fixed multiple of the number of its active corner occurrences (`n=0,-1`), counting repeats. Summing over external frequencies gives at most eight active occurrences for each transfer pair; `Σ_j k_j` is finite. The line replacement therefore costs `O(δ)` uniformly in small `m`. Extending the bare energy integration from `[-W/ν,W/ν]` to the real line costs `O((ν/W)^3)` by the bare majorant `D(x)^4`.

For two static transfers the form is zero. With precisely one static transfer its alternating square vanishes, but the two real-pair terms combine, for ordered frequencies `n,m`, into

```
k_(n−m) sgn(ω_n ω_m) (y_n−y_m)^2 I(w_n,w_m),
I(a,b) = ∫_R a b x² / [(x²+a²)² (x²+b²)²] dx = π/[2(a+b)^3],
```

Here `w_n=|2n+1|`. Differentiate `∫ x²/[(x²+a²)(x²+b²)] dx=π/(a+b)` once in each parameter and divide by 4; local domination holds for `a,b>0`, including at equality. For each choice of static transfer, `(0,−1)` and its reverse contribute `−8 k_1 I(1,1)`. Pairs with exactly one active frequency give `4 Σ_{l≥1}(k_l−k_(l+1)) I(1,2l+1)≤4 k_1 I(1,3)`. Both static choices together give at most `−15π k_1/16`. With two nonstatic transfers the absolute sum is at most `10π (Σ_{j≠0} k_j)^2`: apply the bounds in (4) over corner occurrences (4 alternating and 6 pair contributions per active frequency), then integrate `D(x)^4`, whose integral is `π/2`. Since `Σ_{j≠0} k_j≤5 k_1` and `k_1<1/1000`, the limiting form is below `−π k_1/4`. At the chosen fixed `λ,t`, first choose `δ>0` sufficiently small and **fix it**, then take `m` sufficiently small; the bound is below `−π k_1/8` and `F_X/Q0≤−r0² k_1/128`. The magnitude need not be uniform in `t` or `δ`. The momentum variance is exactly zero under the actual dressed weights, but the ME term need not vanish. The full form still has the stability gap in §3.

## 4. Joint-overlap model: no onset on `t₀(λ) Ω≤τ≤x_λ W/π`

Take a UV-projected single continuum band, active momenta in the unit ball `D⊂R³`, dispersion `ξ(k)=W(2|k|³−1)` there, and absent lines outside (`z=0`). Use measure `d³k/vol(D)` for each integral, a constant-in-momentum Einstein interaction `V_j=K_j/N0`, and `K_j=λ Ω²/[Ω²+(2νj)²]`. The momenta conserve exactly; this is not a local DOS approximation. The sphere `|k|=2^(-1/3)` has nonzero velocity. The radial law `3R² dR` gives `N0=1/(2W)`, `A=1`, and `E_F=W`. The single-particle projection is part of the model, not a claim about a smooth periodic band edge.

**Joint overlap.** For fixed external `k`, draw `a,c` independently in `D` and set `b=a+c-k`, discarding `b` outside `D`. In the three energies, the pushforward measure is bounded above by

```
  dξ_a dξ_b dξ_c / (48 W³ R_a R_b R_c),   R_i=((1+ξ_i/W)/2)^(1/3).      (3)
```

At fixed radii of `a,c` and direction of `a`, `U=cos(angle(c,a-k))` is uniform on `[-1,1]`. Since `R_b²=R_c²+|a-k|²+2R_c|a-k| U`, the energy derivative in `U` is `6W R_b R_c |a-k|`. The angular average of `1/|a-k|` at radius `R_a` is `1/max(R_a,|k|)≤1/R_a`: integrate `(R_a²+|k|²−2R_a|k| U)^(-1/2)` directly against `dU/2`, with the limit at `|k|=0`. Multiplying by the two flat radial densities gives (3). Endpoint singularities are integrable; no atoms occur, since the exceptional zero radii and coincident vector sets have measure zero. This proof works with any dispersion/domain whose joint law has the same separable upper envelope (3), with its constant changed.

**Bounds from the overlap.** Write `x=ν/W`, `R=Σ_j K_j≤λ(1+Ω/(2τ))`, `D_n(ξ)=(ω_n²+ξ²)^(-1/2)`, and

```
J(x)=∫_{−1}^1 ((1+v)/2)^(-1/3) (x²+v²)^(-1/2) dv,
L=(2/(3π²)) R² x J(x)³.
```

On any set of branches or iterates with `|z_{n,k}|≤2D_n(ξ_k)`, (3) bounds `∫ |z_a z_b z_c| d³s d³t` at each frequency pair by `J(x)³/(6W³)` (with the normalized measures); change variables from `s,t` to the independent `a,c`. Consequently `||Σ_X||∞/ν≤L`, and the resolvent difference identity and product telescoping give a crossed Lipschitz bound `6L` in the sup norm of the self-energy. For `x≤1`, split at `v=-1/2`: the negative edge interval contributes `<3`; the remainder contributes `≤4 asinh(1/x)≤4 log(3/x)`. Thus `J≤8(1+log(1/x))`. For `x≥1`, `J≤3/x`. In particular `x J(x)^3→0` as `x→0` and `x J(x)^3≤4096` globally (maximize `e^(-s)(1+s)^3` on `s≥0`). For a square on `l` corner occurrences use `|Σ_i (±y_i)|²≤lΣ_i|y_i|²`, counting coincident occurrences separately. Reanchor each term separately in (1): the alternating square costs 4 and the three pairs cost 6 in total. This yields

```
 |F_X| ≤ 10 τ² ∫ |Q_p| |y_p|² ∫ V_s V_t |z_a z_b z_c| ≤ 10 L Q0.      (4)
```

Each corner carrying a squared trial value also carries a `Q` in (1); keeping that factor, rather than replacing it by `|z|`, is essential off shell. The joint-overlap assumption is on all three internal energies, not merely bounded marginal DOS.

**Existence at any fixed coupling in the adiabatic temperature window.** Fix finite `λ` and a lower temperature ratio `t0>0`; then `R≤R0=λ(1+1/(2t0))` for `τ≥t0 Ω`. First solve the symmetric finite-band rainbow normal equation for `f*` (a frequency-only self-energy). For `n≥0`, use the positive widths `u*_n≥νw_n`, `w_n=2n+1`, and the folded differences `K_{n−m}−K_{n+m+1}≥0`. The width map sends `νw_n≤u*_n≤ν(w_n+R0)` to itself, with weight `b_n^*=(2/π) arctan(W/u*_n)`. It is a sup-norm contraction on width offsets for small `x`, since `|db_n^*/du*_n|≤2/(π W)` and the kernel row sum is at most `R0`; reflect to negative `n` and set `f*_n=i(ω_n−u_n^*)`. The real rainbow average is zero by the symmetric DOS.

For the crossed problem write `Σ=f+g`, with `f=F(f,g)` the momentum-independent rainbow and `g=X(f,g)` the crossed map. Take `||f−f*||≤ν/4`, `||g||≤ν L`, with `L→0` uniformly in `τ≥t0 Ω` as `x→0`. In this box, `|z|≤2D_n` by the resolvent identity relative to the rainbow solution, and `sgn u=sgn ω`. To see that the split map contracts:

* For `g=0` the derivative of the flat-band momentum average with respect to `f_n` is an endpoint difference, bounded by `4N0/W` when `|Re f_n|≤W/2`. Thus the `f`-Lipschitz bound of `F(f,0)` is `4Rτ/W`.
* For fixed `f`, the `g`-Lipschitz bound of `F` is `4R`, using `∫ D_n²≤π N0/ν`. For fixed `g`, the difference between `F(f,g)−F(f',g)` and the `g=0` difference is bounded by `(32/π) R L ||f−f'||`: the per-line product difference of the two resolvents is `≤16 ||g|| D_n³`, and `∫ D_n³≤2N0/ν²`.
* The `X` image bound is `ν L` and its Lipschitz bound is `6L(||δf||+||δg||)`. The `F` image stays in its ball by `||F(f,g)−f*||≤(4Rτ/W)(ν/4)+4Rν L`. Use the norm `max(||δf||,C||δg||)`, `C=max(1,16R0)`: the Lipschitz row bounds are `4Rτ/W+(32/π)RL+4R/C` and `6L(C+1)`, both less than 1 for small `x`.

Contraction gives the normal branch in this box, inherits the reflection symmetries, and continues from zero coupling (all small-`x` guards can be taken uniformly for interactions up to `λ`). It also gives `||f−f*||≤8Rν L` by the rainbow Lipschitz bound with `g=0`. This is a small joint-overlap assumption verified by geometry, not an ansatz for the dressed widths.

**Strict scalar comparison in the window.** Put `I_n=∫_D |Q_{n,k}| dk`, `b_n=I_n/(πN0)` and use angular averages weighted by `|Q|`. Here “angular” means all momentum dependence at fixed frequency, including radius. The baseline `Q0+F_ME` decomposes into the mean and variance in that probability measure. The variance term is at least its share of `Q0`, since it has additional coefficient `d_n=(Σ_{m≥0}(K_{n−m}−K_{n+m+1}) b_m)/w_n≥0` relative to `Q0`. On the mean, after normalization to its share of `Q0`, the form is `I−S(b)`, with

```
 S(b)_{nm} = sqrt(b_n b_m) (K_{n−m}+K_{n+m+1}) / sqrt(w_n w_m) − δ_nm d_n .
```

Indeed the same-sign and opposite-sign terms in `Q_p Q_a` give the difference in `d_n`; in the bilinear pairing term `Q_p y_p` is even and their sum enters. The branch bound, `|Q−Q*|≤|z−z*|`, and `∫ D_n²≤πN0/ν` give `|b_n−b_n^*|≤2(1+8R)L`; the `b_n^*` tend to 1 at every fixed `n` as `x→0`, and never exceed 1. Thus `S(b)−S(1)→0` in norm uniformly for `τ≥t0 Ω` as `x→0`: the diagonal normal tail is `O(R0/w_n)` and each fixed row converges by a summable transfer majorant. The pairing part is dominated in Hilbert–Schmidt norm by the single Toeplitz/Hankel envelope `K_j(t0)` divided by `sqrt(w_n w_m)` (use `1/(w_n w_m)≤(w_n^(-2)+w_m^(-2))/2` and squared summability of both one-index envelopes). This argument avoids assuming that the high-frequency `b_n` are close to 1.

`S(1)` is the scalar flat Einstein operator of [Rens, Eqs. (7)–(8), Theorems 3–5](https://arxiv.org/html/2611.20496v1#S3). Its top is strictly below 1 at `t_c=C∞√λ` for fixed finite `λ`. Continuity gives a common positive gap on `[t_c−η,t_c+η]` for some `0<η<t_c`; beyond it the soft bound `sup σ S(1;t)≤(t_c/t)²` gives a common positive gap as well. Take `t₀(λ)=t_c−η`; by (4) and the norm convergence the full form is strictly positive on the regular branch for `t₀(λ)Ω≤τ≤x_λ W/π` for some `x_λ>0`. This excludes violations at temperatures comparable to `Ω` as `m→0` with fixed `λ`, without any inflation of the ceiling. It is not uniform as `λ→∞` and does not cover `τ/W` of order one.

**An all-temperature weak-coupling instance.** For `0<λ≤2^-20` and any `τ≥√λ Ω/8`, `R≤λ+4√λ<1/128`. Use the ball `||Σ||≤ν/2` centered at zero; the bound `|z|≤2D_n` continues to hold. The bare rainbow has `||F(0)||≤ν R`, and its Lipschitz bound in this ball is `4R` by the same two-line estimate. Since `π>3` and `x J^3≤4096`, `L≤304R²`, so the full image bound relative to `ν` is `3R+L<1/2` and the Lipschitz bound is `4R+6L<1`. For the pairing form in `d=|z|φ`, apply weighted Schur with weight `|z|` on `D`; its row criterion is precisely the absolute raw pairing row. The rainbow row is `≤4R`, the crossed rows sum to `≤6L` by bounding the additional line `|z_v|≤2/ν` in each summand. The kernel is symmetric as in §2. Therefore `τ_lin≤√λ Ω/8<C∞√λ Ω` on the regular branch, with no restriction on `m`. This argument, unlike the finite-coupling adiabatic window, covers arbitrarily high temperatures. It proves only the stated small coupling range.

## 5. Trust base, literature gate, and replay

The skeleton and scalar flat-DOS spectral gap are taken from the cited sources. All momentum bounds, reanchorings, branch constructions, and asymptotic tests here are analytic. `replay_form.py` independently differentiates the literal Nambu matrix by exact odd interpolation, and compares with (1) and a separately implemented unfolded/folded ME identity. It rejects changing the square factor, dropping the real-line terms, changing the middle anomalous sign, omitting the normal subtraction, treating the whole double-static form as zero, and dropping the momentum variance penalty. `replay_bounds.py` checks fixed rational margins and geometric box extrema, clears the margin denominator independently, recomputes the angular-mean integral by partial fractions rather than parameter differentiation, and checks derivative factors and contraction guards. Its mutants alter the nonstatic tail, static witness subset, relative Matsubara signs, and both contraction thresholds. These are exact replay checks, not finite-grid extrapolation, machine verification of the full analytic argument, or evidence of an onset counterexample.

Day-one decision: **KILL** the ultranarrow fixed-Fermi-surface-coupling mechanism; **GO** on the finite-width crossover and on joint-overlap control. [Grimaldi–Pietronero–Scattoni](https://arxiv.org/html/cond-mat/9806012v1#S2.E5) distinguish kinematic limits of a single triangle; [Wang et al., §§2.3–2.4](https://arxiv.org/html/1602.00656v1#S2.SS3) explicitly distinguish the physical mass enhancement from the conventional Fermi-surface average, and treat the elastic perfect-forward vertex by a ladder sum, which is not this truncation. [Gor’kov](https://arxiv.org/abs/1510.03327) treats screened interface LO phonons, not a fixed-`M` crossed-truncation ceiling. The 60-day scan (2026-09-19 through 2026-11-18) located the flat scalar theorem above, [Mishra et al.](https://arxiv.org/abs/2609.30564) on material-specific anharmonic/vertex corrections, and [Chubukov et al.](https://arxiv.org/html/2604.14294v3) on strong-coupling polaron/bipolaron breakdown, but no preemption for the signed identity or these fixed-coupling model statements. This is a targeted literature check, not an exhaustive novelty claim. Paper-reader delegation was unavailable; only legal arXiv HTML/abstracts were used.

What remains: a universal onset comparison even for bounded `A` is open. At finite but not weak coupling, the overlap proof does not cover the full interval of order `W`. The crude high-temperature contraction of §2 does apply once `τ B_eff/ν²≤1/16`, where `B_eff=R/N₀=2WR`, so `τ B_eff/ν²=2R/(πx)→0`. This is an effective weight, not the integral of the constant interaction over unrestricted momentum transfers: the rainbow integrates over `a∈D`, and the crossed term over independently normalized `a,c∈D` after changing variables, with `b∈D` dropped in the crude upper estimate. For general anisotropy and fixed `v_F q_c/Ω`, even the phonon-scale onset is not settled. No comparison here is with an experimental transition, a fluctuation-corrected transition, a different momentum-dependent definition of `M`, or the all-orders electron–phonon problem.
