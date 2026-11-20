# Momentum-resolved first-crossed onset: signed projection and weak-coupling closure

Work in energy units \(\tau=k_BT\), \(\nu=\pi\tau\), \(\omega_n=(2n+1)\nu\). Use the bare phonon, no Coulomb term, fixed bare chemical potential, and precisely the rainbow and first-crossed normal and anomalous skeleton of [G1-SIGNED, §1](../G1-SIGNED/RIGOR.md). In particular all three crossed anomalous insertions are included. The regular branch below means the indicated contraction solution continued from zero interaction at the same temperature, not every possible positive-width solution. An onset bound excludes a pairing spectral value at or above one in the stated temperature region. It does not prove that there is a crossing, that the state is stable to other channels, or anything about the exact transition of the Hamiltonian.
We work in the inversion-even, even-frequency singlet sector, allowing the trial to vary over momentum. A stability bound on this sector in particular covers a momentum-independent gap.

For a nonzero positive spectrum \(d\lambda(\Omega)\) supported in \((0,\Omega_{\max}]\), set
\[
K_j=\int \frac{\Omega^2}{\Omega^2+(2\nu j)^2}\,d\lambda,\quad R=\sum_{j\in\mathbb Z} K_j,\quad M_i=\int\Omega^i\,d\lambda,\quad M=M_2,\quad \rho=\frac{\Omega_{\max}\sqrt\lambda}{W},\quad x=\frac{\nu}{W}.
\]
Here \(\lambda=K_0\), \(M_1\) has units of energy and \(M\) of energy squared. The momentum-independent transfer is \(V_j=K_j/N_0\); the lines still conserve momentum in the crossed skeleton. Two separate bounds are
\[
R\le \lambda+M_1/(2\tau),\qquad R\le\lambda+M/(12\tau^2).
\tag{1}
\]
The first follows from monotone integral comparison on positive transfers, the second by summing \(M/(2\nu j)^2\). We use \(M_1\le\sqrt{\lambda M}\), \(M\le\lambda\Omega_{\max}^2\), and the true spectral \(9/50<C_\infty<1/5\), not its rounded decimal.

## A. PROVED: a linear signed error in the spherical overlap model

Take the same UV-projected unit ball \(D\subset\mathbb R^3\) as G1-SIGNED §4: \(\xi(k)=W(2|k|^3-1)\) on \(D\), absent lines off \(D\), momentum measure Lebesgue divided by \({\rm vol}(D)\). Then \(N_0=1/(2W)\), \(A=1\), and the Fermi-surface coupling is \(K_0\). Let \(\tau_c=C_\infty\sqrt M\), take \(B\ge\max(1,\lambda+3\sqrt\lambda)\), and set
\[
x_*=2^{-80}B^{-8},\qquad C_0=13000B^4.
\]
For \(\tau\ge\tau_c\) and \(0<x\le x_*\) the regular branch exists. On it, decompose the trial \(y=\phi/u\) as \(\bar y+y_\perp\), its momentum mean and remainder in the \(|Q_n|\)-weighted probability measure at each frequency. Here \(z=r-iQ=(i\omega-\xi-\Sigma)^{-1}\), \(u=\omega-\operatorname{Im}\Sigma\), and \(Q_0[y]=\int \omega Q |y|^2\), as in G1-SIGNED. With \(S(1)\) the flat scalar canceled ME matrix at the *same* spectrum and temperature, the full form satisfies
\[
F[y]\ \ge\ [1-\sup\sigma S(1)-C_0 x]Q_0[\bar y]+\tfrac38 Q_0[y_\perp]. \tag{2}
\]
Consequently there is no onset in this window when \(\tau\ge\tau_c(1+52000B^4\rho)\). The constants are sufficient, not sharp. Neither the temperature window nor the constant is uniform as \(B\to\infty\). A4 extends the same estimates to a class of nonconstant momentum transfers.

### A1. Geometry and the signed three-line bound

Write \(v_i=\xi_i/W\), \(R_i=((1+v_i)/2)^{1/3}\); here the subscript distinguishes the radii from the transfer sum \(R\). For fixed external \(k\), let \(a,c\) be independently drawn momenta in \(D\) and \(b=a+c-k\), dropping \(b\notin D\). At \(K=|k|>0\) the density in \(d\xi_a d\xi_b d\xi_c\) is exactly
\[
\frac{T}{48W^3 R_aR_bR_c},\qquad T=\frac{|[|R_a-K|,R_a+K]\cap[|R_b-R_c|,R_b+R_c]|}{2K}. \tag{3}
\]
In this formula the bars around the intersection are its length. To derive it, at fixed \(R_a,K\) the distance \(d=|a-k|\) has law \(d\,dd/(2R_aK)\) on the first interval. Given that distance, the cosine of the angle between \(c\) and \(a-k\) is uniform on \([-1,1]\), and \(d\xi_b/d\cos=6W R_bR_c d\). The admissible distances for a given \(R_b,R_c\) are the second interval. Multiply the two flat radial densities \(1/(2W)\) and integrate in \(d\). Exceptional zero denominators occur on null sets; all estimates use the integrable envelope \(T\le1\), with the bound of G1-SIGNED §4 at \(k=0\). No independent-line approximation is used.

Call \(|v|\le1/2\) the core. At an external core momentum, \(K>1/2\). The intersection length is 2-Lipschitz in each of the three radii by the min/max formula, while \(|dR_i/dv_i|\le 2/3\) in the core. For \(H=T\prod_i R_i^{-1}\), its Lipschitz constant in a single core coordinate \(v_i\), with all others fixed, is at most \(7\prod_{j\ne i}R_j^{-1}\): there \(R_i^{-1}\le2\), its derivative is bounded by 3, and the slope bound of \(T\) is at most 2. This is a one-coordinate estimate and holds also at triangle boundaries.

Use any centered lines \(z_i^*=(iu_i^*-\xi_i)^{-1}\) with fixed nonzero frequency widths, arbitrary in size. Each real factor \(W r_i^*\) is odd and bounded by \(1/|v_i|\). After taking absolute values, its exterior integral with weight \(R_i^{-1}\) is at most 6, since \(\int_{-1}^1 R_i^{-1}dv_i=3\). The weighted integral of \(W|Q_i^*|\) is at most \(2\pi+3<11\): on the core use the full Lorentzian integral and \(R_i^{-1}\le2\); outside, \(W|Q_i^*|\le1/(2|v_i|)\le1\).

In the product of three line components, split real coordinates into core and exterior. For the \(h\ge1\) core real coordinates, fold over all \(2^h\) signs, keeping the odd sign factors; use the difference in the coordinate with the smallest positive \(v_i\). The alternating coefficient is then bounded by \(7\,2^{2h-1}\min_i v_i\prod_{j\notin \mathrm{fold}}R_j^{-1}\). This follows by pairing the \(2^{h-1}\) sign choices, using a difference over distance \(2v_i\), and bounding the other folded radius factors by 2. Layer cake gives
\[
\int_{(0,1/2)^h}\frac{\min_i v_i}{\prod_i v_i}\,d^h v=\int_0^{1/2}\log^h((1/2)/t)\,dt=h!/2.
\]
Thus the fold costs for \(h=0,1,2,3\) are \(1,7,56,672\). Summing the choices of which real coordinates are in the core, the total (for 0, 1, 2, 3 real factors) is at most \(1331,1573,1936,2652\), respectively. For every external core momentum and every assignment of frequencies and components \(r_i^*,Q_i^*\), the absolute value of the *signed spatial integral* of their product is therefore bounded by \(63/W^3\), using \(3000/48<63\). For the complex normal product of three lines use its eight monomials. Taking absolute values of real factors before folding would not prove this bound.

### A2. Branch and normal-width control

For \(\tau\ge\tau_c\), (1) gives \(R\le\lambda+3\sqrt\lambda\le B\). Use \(D_n(\xi)=(\omega_n^2+\xi^2)^{-1/2}\), the symmetric rainbow self-energy \(f^*\) of G1-SIGNED §4, its centered widths \(\nu w_n\le |u_n^*|\le\nu(w_n+R)\), \(w_n=|2n+1|\), and the split crossed branch \(\Sigma=f+g\). With the \(J(x)\) of that proof let
\[
L=\frac{2}{3\pi^2}R^2 x J(x)^3 \le 40 B^2 x(1+\log(1/x))^3
   \le 2^{17} B^2 x^{3/4},\qquad L^2\le x/64\quad(x\le x_*).
\tag{4}
\]
For the middle estimate maximize \(e^{-s/4}(1+s)^3\) for \(s\ge0\); it is at most \(12^3\). The G1-SIGNED contraction estimates apply: the box \(\|f-f^*\|\le\nu/4,\ \|g\|\le\nu L\) has \(|z|\le2D_n\), and the split Lipschitz rows in the norm weighted by \(16B\) are bounded by \(4R\tau/W+(32/\pi)RL+4R/(16B)\) and \(6L(1+16B)\). Both are less than one by (4), and the \(f\)-row image bound \(\nu[(4R\tau/W)/4+4RL]\le\nu/4\) also holds. The rainbow branch itself exists by the positive folded map with width Lipschitz constant \(2R\tau/W\). All bounds hold along the coupling homotopy up to the target coupling. In particular \(\|f-f^*\|\le8R\nu L\) and \(d=\|\Sigma-f^*\|/\nu\le(1+8B)L<1/8\). The ratios \(u/\omega\) are bounded above and away from zero, so coercivity in \(Q_0\) excludes a pairing spectral value \(\ge1\).

The resolvent identity gives \(|z-z^*|\le2\nu d D_n^2\). Telescoping any internal triple product (including the product of complex lines in the normal equation) costs at most \(3d\) times the absolute majorant \(8D_aD_bD_c\), as \(\nu D_n\le1\). After interactions and sums, this is at most \(3dL\) in units of \(\nu\). On the core the signed estimate gives \(|g|/\nu\le256 B^2 x+3dL\le Gx\), \(G=257B^2\). For the conversion, \(\tau^2(R/N_0)^2/(\nu W^3)=4R^2 x/\pi^2\), and \(8\cdot63\cdot4/\pi^2<256\); the improvement \(3dL/x\le3(1+8B)/64\le 27B/64\) fits the final unit of \(B^2\). The global bound is still \(|g|/\nu\le L\).

Uniformly in frequency, \(\int_D D_n^2\le\pi N_0/\nu\) and \(\int_{\rm ext}D_n^2\le 1/W^2\), using the flat DOS and the integral of \(\xi^{-2}\) on \(W/2<|\xi|\le W\). The latter integral is at most \(2x/\pi\) times the former *upper bound*. The rainbow map difference from setting \(g=0\) is bounded by \(4R\nu(G+1)x\) (use \(4D_n^2|g|\) and \(L\le1\)); the rainbow map at \(g=0\) has Lipschitz constant \(4R\tau/W\le1/2\) on the branch box. Hence \(\|f-f^*\|/\nu\le8B(G+1)x\). Applying \(|Q-Q^*|\le|z-z^*|\le2D_n^2(|f-f^*|+|g|)\) with the same split gives
\[
|b_n-b_n^*|\le 2(8B+1)(G+1)x\le5000B^3x=:\varepsilon,\qquad b_n=\frac{\int_D|Q_n|}{\pi N_0} . \tag{5}
\]
These estimates use the dressed weights, not a replacement of the crossed branch by the rainbow branch in the gap equation.

### A3. Projection and coercivity

The normalized ME part is the direct sum of \(I-S(b)\) on momentum means and a variance term bounded below by its share of \(Q_0\), by G1-SIGNED §4. For \(n,m\ge0\), set \(P_{nm}=K_{n-m}-K_{n+m+1}\ge0\). The canceled matrix has off-diagonal weights \(\sqrt{b_nb_m}(K_{n-m}+K_{n+m+1})/\sqrt{w_nw_m}\), and diagonal \([2K_{2n+1}b_n-\sum_{m\ne n}P_{nm}b_m]/w_n\). Since \(b^*\le1\), (5) gives a positive off-diagonal comparison with excess at most \(\varepsilon\) times the flat weights. Including the positive diagonal excess, its Schur bound is \(R\varepsilon\) (\(2K_{2n+1}\le K_0+K_{2n+1}\)); the negative normal diagonal contributes a further \(R\varepsilon\) plus the rainbow deficit. Precisely, for the latter use \(1-b_m^*\le 2\tau(w_m+R)/W\) and the [V1 row identities](../G1-VERTEXBOUND/RIGOR.md), \(\sum_{m\ne n}P_{nm}/w_n\le K_1\) and \(\sum_{m\ne n}P_{nm}w_m/w_n=2\sum_{j\ge1}K_j+K_{2n+1}\le R\). Then
\[
\sup\sigma S(b)\le\sup\sigma S(1)+2R\varepsilon+2(\tau/W)(R+R^2)\le\sup\sigma S(1)+11000 B^4 x . \tag{6}
\]
Here the comparison of off-diagonals is after replacing the trial coordinates by their absolute values, using positivity of their reference weights. The sign of the normal mass term is not dropped.

On \(\bar y\), the trial factors in each crossed square are fixed by the frequencies alone. For each transfer pair one can take the absolute value *after* the spatial integrals and use the four-square and two-square estimates of G1-SIGNED §4. Reanchor each squared occurrence within the integrated coefficient before bounding it: the new external corner still carries a \(Q\) factor. This gives the same multiplicity 10 as in G1-SIGNED (four alternating and six pair contributions). On the external core the summed signed internal coefficient divided by \(\nu\) is at most \(28 B^2x+3dL\), by A1 for a monomial and A2 for its replacement error. Outside the core the absolute bound is \(L\), multiplied in the form normalized by \(Q_0\) by the exterior fraction of \(\int |Q_n|\) divided by \(w_n\). That fraction divided by \(w_n\) is at most \(14(B+2)x\). Indeed, if \(|u_n^*|\le W\), then \(\int|Q_n|\ge \tfrac12\int |Q_n^*|\ge\pi N_0/4\) (use \(d<1/8\) in the resolvent ratio), and the exterior integral is at most \(4\nu(w_n+B+1)/W^2\). If \(|u_n^*|>W\), use the trivial fraction one and \(1/w_n<(1+B)x\). Thus
\[
|F_X[\bar y]|\le 10[28 B^2 x+3dL+14(B+2)x L]\,Q_0[\bar y]\le1000B^2x Q_0[\bar y]. \tag{7}
\]
This operation does not claim a pointwise bound on a crossed triangle or a small absolute kernel in the momentum-independent trial sector.

For arbitrary trials the absolute form bound is \(|F_X[y]|\le\delta Q_0[y]\) with \(\delta=10L\), by the reanchoring bound of G1-SIGNED §4. Polarization and the norm bound for the self-adjoint form give a mixed term bounded by \(2\delta\sqrt{Q_0[\bar y]Q_0[y_\perp]}\). Absorb it by spending \(\tfrac12 Q_0[y_\perp]\) plus \(2\delta^2 Q_0[\bar y]\). The remaining variance coefficient is at least \(3/8\) since \(\delta\le1/8\), and \(2\delta^2\le4x\) by (4). Combining (6) and (7) gives (2). The estimates start on finitely supported frequency trials and extend by density in \(Q_0\), with no Matsubara cutoff.

For the inflation statement, use \(\sup\sigma S(1)\le(\tau_c/\tau)^2\) by the scalar soft-kernel comparison used in V1. At \(\tau\ge2\tau_c\) the gap is at least \(3/4\) before the error \(C_0x\le C_0x_*<1/4\). For \(\tau\le2\tau_c\), \(x\le2\pi C_\infty\rho<2\rho\), and if \(\tau\ge(1+4C_0\rho)\tau_c\) then \(1-(\tau_c/\tau)^2>(4C_0\rho)/2\), using \(4C_0\rho\le1\). This proves the claimed window bound. The earlier \(9/32\) margin in the forward example of G1-SIGNED is model-specific and cannot replace the variance estimates here.

### A4. A nonconstant transfer with balanced angular averaging

On \(D\times\mathbb T\) use normalized Haar measure on the spectator circle \(\mathbb T\), dispersion \(\xi(k,\theta)=\xi(k)\), and \(V_j(q,q_\theta)=K_j v(q_\theta)/N_0\) with \(v\ge0\) even, \(\int_{\mathbb T}v=1\), and \(\eta=\|v-1\|_1\le1/(16B)\). For example \(v(\theta)=1+\beta\cos\theta\) is permitted for \(0<\beta\le1/(16B)\). The extra momentum still conserves at each vertex. Its Fermi-surface double average of \(v\) is one, so it does not change \(\lambda\) or \(M\). This separable geometry is an extra assumption, not a model for general anisotropic dispersions.

For the split normal map the rainbow part \(f\) now depends on \(\theta\) but not on \(k\). The contraction estimates of A2 are unchanged: take the spatial estimates at fixed angular variables, uniformly in those variables, then integrate one or two \(v\)'s against the normalized angular measures, with unit mass. By translation symmetry and uniqueness in the box the resulting normal lines are independent of \(\theta\), and coincide with the lines of A2. In the crossed form the full-space mean uses both momentum variables; (7) is unchanged on this mean since its trial is independent of both variables, and the absolute form bound is unchanged for all trials by reanchoring and the unit angular masses.

For the ME part the angular convolution \(T_v\) is the identity on constants and has norm at most \(\eta\) on zero-mean functions, by applying the \(L^1\)-\(L^2\) convolution inequality to \(v-1\). The pairing bilinear term factors into this convolution and the positive folded frequency matrix with entries \(\sqrt{b_nb_m}(K_{n-m}+K_{n+m+1})/\sqrt{w_nw_m}\), applied to the spatial \(|Q_n|\)-weighted means at each angle. Taking the spatial mean is a contraction. The frequency norm is at most \((1+\varepsilon)R\) by Schur; \(\varepsilon\le1\) in (5). The normal diagonal is independent of angle and adds a nonnegative term. As a result the ME form still separates the full mean from the variance, with the same \(I-S(b)\) on the mean, and with variance coefficient at least \(1-\eta(1+\varepsilon)R\ge7/8\). In the mixed crossed bound of A3 spend \(3/8\) of the variance and \((8/3)\delta^2\le5x\) on the mean. After the direct variance loss \(\delta\le1/8\), the remainder is again \(3/8\). This proves (2) and the same inflation bound for these nonconstant transfers.

The weak estimate of B also applies to this separable model: the triple Hölder estimate is applied to the three spatial lines with each angle fixed, followed by integration of the unit-mass angular weights. The two spherical temperature regions in B2 can therefore be patched also in this class for the uniform weak-coupling result. The original \(t_0(\lambda)\) is used only in the strictly spherical Einstein subcase.

## B. PROVED: weak-coupling high-temperature closure by a different estimate

Take an inversion-invariant active support on a continuum projected momentum domain (or an additive momentum group with Haar measure), with \(\xi(-k)=\xi(k)\) almost everywhere. Assume the even marginal DOS has total measure one and is bounded by \(A N_0\), using the finite positive reference \(N_0=N(0)\). Use the same momentum-independent \(V_j=K_j/N_0\), and extend lines by zero when projected out. No joint density assumption is required here. At any \(\tau>0\) with \(r=A R\le1/9\), the normal skeleton is a contraction in \(\|\Sigma\|_\infty\le\nu/3\), continued from zero, and the absolute Schur bound for the symmetrized pairing operator is strictly below one.

To prove it use \(h=1/3\) and \(|z|\le D_n/(1-h)\) in the ball. The normal rainbow image divided by \(\nu\) is at most \(r/(1-h)\): the bare rainbow image is at most \(r\) by DOS symmetry and the Lorentzian integral, and its difference from the dressed image uses \(|z-z_0|\le h\nu D_n^2/(1-h)\). Its Lipschitz constant is at most \(r/(1-h)^2\). For the crossed term, for nonnegative \(u(a),v(c),w(a+c-k)\) the spatial triple integral is at most \(\|u\|_{3/2}\|v\|_{3/2}\|w\|_{3/2}\). This follows directly by applying Hölder with exponents 3 to the three factors \((uv)^{1/2},(uw)^{1/2},(vw)^{1/2}\) and changing variables in the latter two pair integrals. For any of the bare factors, \(\int D_n^{3/2}\le A N_0 J_{3/2}/\sqrt\nu\), where \(J_{3/2}=\int_{\mathbb R}(1+t^2)^{-3/4}dt<6\) by splitting at \(|t|=1\) and integrating \(|t|^{-3/2}\) outside. Thus, with \(c=(J_{3/2}/\pi)^2<4\), the crossed image divided by \(\nu\) is bounded by \(c r^2/(1-h)^3\), and telescoping the difference of three resolvents yields the Lipschitz bound \(3c r^2/(1-h)^4\). These scalings include the two powers of \(\tau\) and two powers of \(N_0^{-1}\). At \(r=1/9\), replacing \(c\) by 4 gives the limiting total image bound \(1/6+1/6=1/3\), and limiting total Lipschitz bound \(1/4+3/4=1\). Since \(c<4\), both true bounds are strict. Uniqueness in the ball preserves Matsubara conjugation and any momentum symmetry of the bare problem. For the pairing kernel in the variable \(|z|\phi\), Schur with weight \(|z|\) gives the same bound as this normal Lipschitz estimate: the rainbow row has one internal \(|z|^2\); each of the three crossed rows has the triple product times one extra \(|z_v|\le1/[(1-h)\nu]\). Normal and anomalous BCS factors are not interchanged. This proves the assertion.

### B1. Explicit all-temperature weak-coupling bound

If \(A\lambda<1/9\), the criterion holds at every \(\tau\ge \tau_Y=9A M_1/[2(1-9A\lambda)]\) by (1). In particular \(\tau_{\rm lin}\le\tau_Y\) on the regular branch, with no small-bandwidth-ratio condition. For \(A\lambda\le 1/26^2\), \(\tau_Y\le (9/2)\sqrt{A\lambda}\sqrt{AM}/(1-9A\lambda) < (9/50)\sqrt{AM}<C_\infty\sqrt{AM}\) by Cauchy and the elementary rational guard at \(1/26\). This also works for multispectral phonons. It does not exclude pairing at temperatures below \(\tau_Y\).

### B2. Completion of the spherical window above the onset bound

For the Einstein model in G1-SIGNED §4, at fixed \(0<\lambda<1/9\) choose a positive \(x_h\) below the upper end of its small-\(x\) stability window so small that \((1+8B_0)L<1/9\) there, with \(B_0\) the bound on \(R\) for \(\tau\ge t_0(\lambda)\Omega\). For sufficiently small \(\Omega/W\), \(t_0(\lambda)\Omega<x_h W/\pi\) and (1) ensures \(R<1/9\) for \(x\ge x_h\). At the overlap the low-\(x\) branch obeys \(\|\Sigma\|/\nu\le R+(1+8B_0)L<1/3\); it therefore coincides with the Young-ball solution. The same strict estimates persist in a neighborhood of the joining temperature. This closes the upper window for all \(\tau\ge t_0(\lambda)\Omega\), not just temperatures of order \(\Omega\); the bound is strictly below the ceiling. It says nothing about a crossing below \(t_0\Omega\).

There is also an uninflated, uniform-in-\(\lambda\) statement on any fixed bounded weak-coupling range \(0<\lambda\le\Lambda<1/9\) in the spherical model, with general spectra. Use the V1 scalar strict finite-\(\lambda\) gap at \(\tau_c\), uniform over spectra at fixed \(\lambda,M\) by the scalar Einstein optimality and monotonicity in temperature. The resulting \(\delta_\lambda=1-\lambda h((2\pi C_\infty\sqrt\lambda)^{-1})>0\) is continuous for each positive \(\lambda\) (the scalar operator is norm-continuous on compact positive temperature ratios). It tends to one as \(\lambda\to0\) since the positive folded pairing norm at the Einstein ceiling is at most \(R\le\lambda+3\sqrt\lambda\), before normal subtraction. Thus \(\delta_*:=\inf_{(0,\Lambda]}\delta_\lambda>0\). In (2) take \(B=1+\Lambda+3\sqrt\Lambda\) and \(x_h=\min(x_*,\delta_*/(2C_0))\), making the low-\(x\) form positive on \(x_c\le x\le x_h\), where \(x_c=\pi\tau_c/W\). For the upper criterion impose
\[
\rho < \frac{x_h}{\pi}\sqrt{12(1/9-\Lambda)}, \tag{8}
\]
Since \(M/W^2\le\rho^2\), (1) gives \(R<1/9\) for \(x>x_Y=\pi\rho/\sqrt{12(1/9-\Lambda)}\). Under (8), \(x_Y<x_h\), and \(x_c\le\pi C_\infty\rho<x_h\) since \(C_\infty\sqrt{12(1/9-\Lambda)}<1\). Thus the strict overlap
\[
\max(x_c,x_Y)<x<x_h
\]
is nonempty. Throughout it the small-\(x\) normal branch has sup norm below \(\nu/3\) by A2 and \(R<1/9\), including along the coupling homotopy. Young-ball uniqueness therefore identifies the same branch. Hence \(\tau_{\rm lin}\le\tau_c\) on that branch for the range (8). This is not a claim uniform in unbounded coupling.

## Trust boundary and adversarial replay

- **PROVED, analytic.** The imports are the momentum skeleton, signed identity and absolute form/contraction estimates from G1-SIGNED; and the centered scalar comparison, strict Einstein gap and V1 row identities from G1-VERTEXBOUND. Sections A1 and B prove new phase-space estimates in different ways: odd signed integration on a regular spherical joint law versus absolute \(L^{3/2}\) control with only the marginal DOS. All factor, bandwidth and temperature scalings concern the specified projected continuum or Haar-normalized models.
- **CERTIFIED finite checks only.** `replay.py` compares the signed density normalization with a separately implemented distance-slicing calculation, checks the core-fold arithmetic and row identity by exact finite sums, verifies the rational guards and the mixed-term bound, and rejects targeted mutants. It is not a continuum or normal-branch certificate. The constants and conclusions above do not rely on a numerical eigenvalue or a finite grid.
- **OPEN.** For arbitrary momentum-dependent \(V_j(q)\) and bounded DOS alone, the signed joint law (3), the scalar/variance decomposition, and the Young argument with constant transfer do not apply. The forward example in G1-SIGNED rules out simply replacing A3 by \(F_X\ge-o(1)Q_0\), not the onset inequality. Even on the sphere, at \(\lambda\ge1/9\) the weak-coupling guard does not cover the intermediate \(x>x_*\) region, and the small-\(x\) estimates are not uniform as \(\lambda\to\infty\). These are gaps in the proof, not evidence of a violation.
