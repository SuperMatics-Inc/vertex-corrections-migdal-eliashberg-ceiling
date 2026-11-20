# A ceiling for the local crossed truncation, and a pointwise lattice obstruction

## 1. PROVED — self-consistent flat-DOS local skeleton

Use energies, including \(\tau=k_BT\). Normalize the per-spin per-cell DOS to one state and take \(N(\xi)=N_0\mathbf1_{[-W,W]}(\xi)\), \(N_0=(2W)^{-1}\), so \(A=1\). Set the chemical potential at half filling with any constant Hartree offset absorbed. Use no Coulomb term, a momentum-independent interaction and a prescribed nonnegative spectral measure \(d\lambda(\Omega)\) on \((0,\Omega_{\max}]\), of positive finite mass \(\lambda\). Define
\[
K_j=\int\frac{\Omega^2}{\Omega^2+(2\pi\tau j)^2}\,d\lambda(\Omega),\qquad
M=\int\Omega^2\,d\lambda(\Omega),\qquad V_j=K_j/N_0,\qquad L=\sum_{j\ge1}K_j.
\]
For an Einstein oscillator \(\lambda=2N_0g^2/\Omega\) and \(M=\lambda\Omega^2\). Thus \(V_j\) has units of energy, and \(M\) units of energy squared. In the **local skeleton** the crossed self-energy of [Schrodi–Oppeneer–Aperis, Eq. (4)](https://arxiv.org/html/1911.12872v2#S2.E4) is evaluated with three separately DOS-integrated Green functions (no remaining momentum correlations). Earlier work includes the local conserving treatment of [Miller–Freericks–Nicol, §II](https://arxiv.org/html/cond-mat/9805254v1#S2) and the numerical \(T_c\) study with a half-filled constant DOS and local vertex of [Paci et al.](https://doi.org/10.1103/PhysRevB.65.012512). Explicitly, with \(a=n-j,\ b=n-j-l,\ c=n-l\), the skeleton used here is
\[
 \Sigma_n=\tau\sum_j V_j\sigma_3\bar G_a\sigma_3
  +\tau^2\sum_{j,l} V_jV_l\sigma_3\bar G_a\sigma_3\bar G_b\sigma_3\bar G_c\sigma_3.
\tag{1}
\]
Here \(\bar G\) is the integral against \(N\), not a Green function divided by \(N_0\). This is the model definition; it is neither an equality for a generic finite-dimensional lattice nor a claim that the Fermi-surface factor in Schrodi et al. Eq. (20) equals \(N_0\). Normal and pairing equations both use (1), with the same fixed phonons.

**Theorem.** Restrict (1) to centered, particle-hole-symmetry-preserving normal branches with zero real shift \(\chi_n=0\) and positive reflected widths. In the even-frequency singlet channel, define onset as the supremum of temperatures of linear pairing instability (including a pairing eigenvalue equal to one does not change the upper bound). Set \(m=\Omega_{\max}/W\) and \(\rho=m\sqrt\lambda\). Then
\[
 k_BT_{\rm lin}\le C_\infty\sqrt{AM}\,[1+30\{(1+\lambda)\rho+\rho^2\}],\qquad A=1 .
\tag{2}
\]
A centered normal solution with positive reflected widths exists at each \(C_\infty\sqrt M\le \tau\le2C_\infty\sqrt M\) if \((1+\lambda)\rho\le1/100\). Neither uniqueness, global branch continuation, nor real-axis causality is asserted. At fixed \(\lambda\), the correction in (2) vanishes as \(m\to0\); the quantitative joint regime for a small bound is \((1+\lambda)\sqrt\lambda\,m\ll1\), not merely \(\lambda m\ll1\) if \(\lambda\) itself diverges.

### Exact elimination

Let \(w_n=2n+1\), \(u_n\) be the renormalized signed frequency, \(s_n=|u_n|/(\pi\tau)>0\), and \(\sigma_n=\mathrm{sgn}(w_n)=\mathrm{sgn}(u_n)\), with \(u_{-n-1}=-u_n\). The hypothesis \(\chi_n=0\) makes the normal \(\sigma_3\) integrand odd in \(\xi\). Its DOS integral vanishes, and (1) preserves the centered ansatz. Oddness of \(u_n\) under reflection alone does not force \(\chi_n=0\) or exclude other normal solutions. Put
\[
 b_n=\frac2\pi\arctan(W/(\pi\tau s_n)),\qquad q_n=\sigma_n b_n,\qquad y_n=\phi_n/u_n,\qquad \eta=\pi^2\tau N_0=\frac{\pi^2\tau}{2W}.
\]
The real even-frequency gap gives odd \(y\); the linearized operator complexifies. The integrated Nambu function is \(\bar G_n=-\pi N_0 q_n(iI+y_n\sigma_1)+O(\phi_n^2)\), and Dyson gives \(\Sigma_n=i(\pi\tau w_n-u_n)I+\phi_n\sigma_1\). Moving the \(\sigma_3\)'s to conjugate the first and third factors gives the following dimensionless normal and linear gap equations:
\[
 \frac{u_n}{\pi\tau}=w_n+\sum_j K_jq_a-\eta\sum_{j,l}K_jK_l q_aq_bq_c,
 \quad
 \frac{\phi_n}{\pi\tau}=\sum_j K_jq_a y_a-\eta\sum_{j,l}K_jK_l q_aq_bq_c(y_a-y_b+y_c).
\tag{3}
\]
In particular the middle insertion has the opposite sign from both end insertions. The linearized normal equations have no linear response in \(\phi\). Subtract the gap equation from \(y_n\) times the normal equation and test against \(q_n\bar y_n\). Summing gives the ME form plus a **signed** crossed term. Reanchoring the parallelogram at each of its four vertices (negating the corresponding step directions) preserves \(K_jK_l q_nq_aq_bq_c\); the alternating difference changes sign at each adjacent vertex. Average the four expressions to get exactly
\[
 \mathcal Q_X=-\frac\eta4\sum_{n,j,l}K_jK_l q_nq_aq_bq_c\,|y_n-y_a+y_b-y_c|^2 . \tag{4}
\]
All terms with either \(j=0\) **or** \(l=0\) vanish. This does not say the corresponding normal terms vanish, or that a zero-frequency but nonzero-momentum transfer in the lattice problem gives zero. For finite-support \(y\), all sums here are absolutely convergent. Use
\[
 |y_n-y_a+y_b-y_c|^2\le4\bigl(|y_n|^2+|y_a|^2+|y_b|^2+|y_c|^2\bigr).
\]
The four corners are counted with multiplicity, including when two indices coincide. Reanchor each of the four summands separately and bound each of the other three width factors by one. With \(\sum_{j,l\ne0}K_jK_l=4L^2\), this gives
\[
 |\mathcal Q_X|\le 16\eta L^2\sum_n b_n|y_n|^2\le4X\|x\|^2,\qquad X=4\eta L^2,\quad x_n=\sqrt{|w_n|b_n}\,y_n.
\tag{5}
\]
These forms extend by density; \(\|x\|\) in (5) is on the full signed line. On the odd subspace, the two halves have identical norm contributions and the same normalized bound holds on the positive half-line. The possible large terms \(K_0^2\) and \(K_0L\) were not bounded by absolute value in the eliminated pairing equation.

### Comparing the ME part on the *crossed normal branch*

On the positive half-line let \(D_{nm}=K_{|n-m|}\), \(H_{nm}=K_{n+m+1}\), \(P=D-H\ge0\) entrywise. With the \(x\)-coordinates above the ME form is \(I-S(b)\), with off-diagonals \(\sqrt{b_nb_m}(D_{nm}+H_{nm})/\sqrt{w_nw_m}\) and diagonal \((2H_{nn}b_n-\sum_{m\ne n}P_{nm}b_m)/w_n\) (cf. [G3, Eqs. (1)–(4)](../G3/proofs.md)). This uses width weights from (3), **not** weights from a separate ME normal solve. Off-diagonals are no larger than those of \(S(1)\). Its excess diagonal is at most \(w_n^{-1}\sum_{m\ne n}P_{nm}(1-b_m)\).

For \(n\ge0\), the terms in the crossed normal row with \(j=0\) or \(l=0\) contribute
\(-\eta\lambda b_n(\lambda b_n^2+2\sum_{j\ne0}K_j b_{n-j}^2)\le0\); the remaining normal row is bounded in absolute value by \(X\). Thus \(s_n\le w_n+\lambda+2L+X\). The exact flat-band integral gives \(1-b_n=(2/\pi)\arctan(\pi\tau s_n/W)\le2\tau s_n/W\). For the off-diagonal normal rows,
\[
 \sum_{m\ne n}P_{nm}\le w_n K_1,\qquad
 \sum_{m\ne n}P_{nm}w_m=w_n(2L+K_{2n+1}).
\]
The latter identity groups the upward and downward weights at each positive transfer (net \(2w_n\)), plus the omitted Hankel diagonal (an extra \(w_n\) at \(2n+1\)); the endpoint discrepancy tends to zero since \(jK_j=O(1/j)\). Hence by the diagonal bound and replacing a vector by its absolute values for the off-diagonals,
\[
 \sup\sigma S(b)\le\sup\sigma S(1)+D_*,\qquad
 D_*=\frac{2\tau}{W}\big[2L+K_1+(\lambda+2L+X)K_1\big]. \tag{6}
\]
The sharp scalar soft-kernel comparison of [Rens, arXiv:2611.20496](https://arxiv.org/abs/2611.20496) gives \(\sup\sigma S(1)\le C_\infty^2M/\tau^2\): apply the scalar comparison to \(K_j\) and its soft saturation \(M/(2\pi\tau j)^2\) for \(j\ge1\); \(K_0\) cancels in \(S(1)\). The soft operator scales linearly with \(M/\tau^2\). This is the only imported spectral bound. In particular, no comparison principle for signed crossed interactions is assumed.

### Uniform temperature control and normal existence

For \(\tau\ge \tau_0=C_\infty\sqrt M\), write \(r=\tau_0/\tau\in(0,1]\). The moment envelope and the Basel sum give \(L\le M/(24\tau^2)\le(4/3)r^2\), and \(\tau L/W\le \rho r/(24 C_\infty)\le\rho r/4\) by \(M\le\lambda\Omega_{\max}^2\) and \(C_\infty>9/50\). Using \(\pi^2<10\), \(X\le(20/3)\rho r^3\). From (5)–(6) and \(K_1\le L\),
\[
 D_*+4X\le r\left[\rho(\lambda/2+59/2)+(10/3)\rho^2\right]
 \le r\epsilon,\qquad \epsilon=30[(1+\lambda)\rho+\rho^2].
\]
Thus the normalized pairing form is bounded below by \((1-r^2-r\epsilon)I\). If \(\tau\ge\tau_0(1+\epsilon)\), then \(r(1+\epsilon)\le1\) and \(r<1\), so the lower bound is positive. The raw symmetrized gap form \(I-B\) is obtained by changing the coordinate from \(\sqrt{|w_n|b_n}y_n\) to \(\sqrt{s_nb_n}y_n\). This is a bounded invertible diagonal change at each fixed positive branch: (3) uniformly bounds \(s_n-|w_n|\) and each \(s_n>0\), so the width ratios tend to one. The estimates therefore exclude linear onset above (2) without assuming monotonicity in temperature.

At each fixed finite \(\lambda\), the same argument also recovers the **uninflated** ceiling for all sufficiently small \(m\). Let \(h(v)\) be the scalar Einstein top eigenvalue per unit coupling as in [Rens, Eqs. (6)–(7) and Thms. 4–6](https://arxiv.org/html/2611.20496v1#Thmtheorem5). Then \(\delta_\lambda=1-\lambda h((2\pi C_\infty\sqrt{\lambda})^{-1})>0\) by the strict scalar bound. Einstein optimality at matched \(M,\lambda\) and the decreasing scalar eigenvalue in temperature give \(\sup\sigma S(1)\le1-\delta_\lambda\) for \(\tau\ge\tau_0\), for every spectrum of those moments. Consequently \(\epsilon<\delta_\lambda\) implies stability for all \(\tau\ge\tau_0\). This additional criterion does not assume monotonicity of the crossed eigenvalue.

For the stated existence region, \(\eta\le 2\rho\) follows from \(\tau\le2\tau_0\), \(M\le\lambda\Omega_{\max}^2\), \(C_\infty<1/5\) and \(\pi^2<10\). Thus \(\eta(\lambda+4L)<1\) and \(X\le1/15\). In the normal map combine the term \(+\lambda b_n\) with the nonpositive crossed terms involving a zero transfer. The combined term is nonnegative. The remaining ME part is at least \(-K_1\), because its off-diagonal compressed entries are nonnegative. Since \(K_1+X<4/5+1/15<7/8\), the map takes the positive box \([w_n-7/8,w_n+\lambda+2L+X]\) into itself. To apply Brouwer, vary only the first finitely many positive coordinates in this box, fix the exterior at \(w_n\), and reflect. Diagonal subsequence extraction and dominated convergence with the summable one- and two-transfer envelopes pass to (3) on the full line. This establishes a centered normal solution at each temperature in the claimed window, not a choice of branch or real-axis causality.

## 2. PROVED — a bounded DOS does not imply a uniform pointwise Migdal factor

The nonvanishing dynamic forward limit is established in the nonadiabatic literature; see [Grimaldi–Pietronero–Strässler I](https://doi.org/10.1103/PhysRevB.52.10516), [II](https://doi.org/10.1103/PhysRevB.52.10530), and the box-DOS discussion of [Grimaldi–Pietronero–Scattoni, §II](https://arxiv.org/html/cond-mat/9806012v1#S2). See also [Chubukov et al.](https://arxiv.org/abs/2004.01281) for the distinction between special forward kinematics and typical transfers. Here we use the triangle of Grimaldi–Pietronero–Scattoni Eq. (5) to exhibit a fixed-\(\tau/\Omega\) obstruction and an explicit bare Cooper-block witness. Take an Einstein mode, fixed \(\lambda>0\), \(\tau=\Omega\), external fermionic indices \(-1,0\), and zero spatial transfer. For the flat example take \(\xi(k)=W(2|k|/\pi-1)\) on the torus \(-\pi\le k\le\pi\) with Haar measure. It is even and has \(A=1\), \(\epsilon_F=W\). The dimensionless triangle (positive \(V\) convention of (1)) is
\[
 P(0;-1,0)=\tau\sum_{p\in\mathbb Z} V_{-1-p}\int
 \frac{N(\xi)\,d\xi}{(iu_p-\xi)(iu_{p+1}-\xi)}. \tag{7}
\]
One may use bare \(u_p=\pi\tau w_p\) or the scalar flat-band ME positive normal branch.
The integral for \(p=-1\) is positive and contributes \(\lambda b_0/s_0\).
For all other \(p\), the frequencies have the same sign. Write the integral as a
divided difference of \(g(u)=-2iN_0\arctan(W/u)\), namely
\([g(u_p)-g(u_{p+1})]/[i(u_{p+1}-u_p)]\). Including the equal-width limit,
each such integral is nonpositive and has magnitude at most \(2N_0/W\), by the
derivative of \(g\). Because \(L\le\lambda\Omega/(4\tau)\) by integral comparison,
the negative sum has magnitude at most \(2\tau(\lambda+2L)/W\le3\lambda m\).
For bare lines \(b_0\ge1-2m\), giving \(P(0;-1,0)\ge\lambda(1-5m)\).
On the ME branch \(w_n\le s_n\le w_n+\lambda+2L\) since \(P_{nm}\ge0\);
for every fixed \(n\), \(b_n\to1\).
The first ME row gives \(s_0=1+\sum_{h\ge0}(K_h-K_{h+1})b_h\to1+\lambda\)
by dominated convergence and telescoping, so the dressed \(P\to\lambda/(1+\lambda)\).
A positive scalar ME branch exists by the same finite-box argument as above with
the crossed terms set to zero and the box \([w_n,w_n+\lambda+2L]\).
In either case \(\lambda\Omega/\epsilon_F=\lambda m\to0\).
At each \(m>0\), the triangle is continuous in spatial transfer by dominated convergence,
so the obstruction also applies to an essential-supremum estimate.

Only for the **bare flat-band example** we also show that a termwise absolute Schur bound on this individual end insertion, weighted by the Cooper lines, cannot vanish with \(m\). Let \(v=2W/\pi\), the Lipschitz constant of \(\xi(k)\). By the resolvent identity and Cauchy–Schwarz, for small \(q\),
\[
 |P(q;-1,0)-P(0;-1,0)|\le(\lambda+2L)\,\frac{v|q|}{\pi\tau}.
\]
Here one resolvent is estimated in \(L^\infty\) by \(1/(\pi\tau)\), and each of the remaining two in \(L^2\) by \((N_0/\tau)^{1/2}\). Take \(m\le1/20\) and the patch \(I\) around \(k_F=\pi/2\) on which \(|\xi|\le\pi\tau/16\). For \(k,k'\in I\), \(v|k-k'|\le\pi\tau/8\), and the previous lower bound and \(2L\le\lambda/2\) give \(|P(k-k';-1,0)|\ge\lambda/2\). On \(L^2(dk/(2\pi))\) the absolute symmetrized frequency block has kernel \(\tau V_1 |G_{-1}(k)|\,|P(k-k';-1,0)|\,|G_0(k')|\). Its norm is at least its Rayleigh quotient on \(\mathbf1_I |G_0|\), namely at least
\[
 \frac{\lambda V_1}{2}\,\tau\int_I |G_0|^2\,\frac{dk}{2\pi}
 =\frac{\lambda K_1}{2\pi}\arctan(1/16)>0,\qquad K_1=\frac{\lambda}{1+4\pi^2}.
\]
The bound is independent of \(m\): although the patch shrinks, the symmetrized pairing measure is not uniform. The dispersion can be extended by independent spectator momenta without changing this example. It does not give a dressed block bound or a lower bound after adding the other signed insertions and normal contribution.

The pointwise obstruction also holds in 3D. For bare lines in the simple cubic band \(\xi=(W/3)(\cos k_x+\cos k_y+\cos k_z)\), the DOS is a convolution of three arcsine densities in \(L^{3/2}\). Young's inequality followed by translation continuity in \(L^{3/2}\) makes it bounded and continuous; it is positive at zero, so \(A\) is fixed and finite and \(\epsilon_F=W\). At each fixed \(p\), \(g(u_p)\to-i\pi N_0\operatorname{sgn}(w_p)\) as \(m\to0\) at fixed \(W\), by evenness and the Poisson approximate identity. The divided difference in (7) has denominator \(2\pi i\tau\), so the same-sign terms tend to zero after multiplication by \(\tau\), whereas \(p=-1\) tends to \(\lambda\). The series limit is justified by Cauchy–Schwarz, \(\|G_p\|_{L^2}^2\le \pi \|N\|_\infty/|u_p|\), and the summable \(K_{-1-p}\) at \(\tau/\Omega=1\). Continuity in spatial transfer again rules out an essential-supremum estimate.

Thus bounded one-particle DOS alone does not imply the proposed pointwise factor for these bare or ME propagators. The cubic and ME-dressed flat examples are pointwise only. These examples rule out neither an estimate restricted to a fully crossed self-consistent lattice branch nor a bound on the signed summed pairing form. They do not disprove a lattice ceiling or prove an unbounded onset ratio.

## 3. Trust / adversarial map

- **PROVED, analytic:** the local bound imports the scalar comparison and soft normalization from the cited preprint, not an unproved bound on a vertex. Equations (3)–(6), the normal box, and the triangle obstruction are infinite-sum arguments with no Matsubara cutoff or interval certificate.
- **Vertex-specific contribution:** on the matched centered crossed normal widths, the three gap insertions plus the normal mass yield exact signed cancellation of either static transfer and a uniform form estimate. The finite-temperature absolute-block witness is separate. Neither the local model nor the nonvanishing dynamic forward limit is claimed as new, and the literature check does not establish priority.
- **Formal scope:** `Core.lean` isolates the Nambu polynomial, four-corner average, static identities, four-square estimate and final real-algebra guard. The Hilbert-space estimates, branch existence, dispersion and integral arguments remain analytic.
- **Independent replay:** `check_exact.py` verifies the corner identity with signed weights on a finite cyclic index set and the ME row arithmetic with a finite-support monotone kernel; it is not used as a substitute for the infinite argument. It also checks the rational constant guards and reproduces the G1 ratio from the rounded inputs.
- **Adversarial checks:** treating crossed pairing without crossed normal leaves a nonzero static diagonal; changing the middle sign spoils the cancellation even for constant \(y\). The alternating vector \((1,-1,1,-1)\) saturates the four-square estimate, so neither the factor of four in the estimate nor the factor of one quarter from reanchoring may be discarded by convention. The flat and cubic examples test the pointwise hypothesis independently of the local closure. No finite-model proof is used to infer an exact Hamiltonian result or the size of omitted higher diagrams.
