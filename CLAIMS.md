# CLAIMS.md, version 3.0.0

Statement-to-record map for the released manuscript. The manuscript is `paper/main.tex` (sha256
`9a4944297579fbdfdaa5590e59a0f7d06b504f7c729de60057fd7405c4401e41`, 118647 bytes) and its PDF is
`paper/main.pdf` (sha256 `ee9f15012993dd2b445c85cc8e0f0cc3868d7ee240fbc013a752509af3ea60d6`, 299526 bytes,
31 pages). The released version 2 manuscript, which is the baseline of this text, is `paper/main.tex` at tag
`v2.0.1` of this repository. The certificate record cited below is `records/R7-flat-sheet-certificate/`, and the
hash of every file of it is in `BUILD.md` section 3.

Labels have the same meaning as in the released versions:
- **PROVED**: a complete analytic proof that has passed the independent review recorded with the supplementary
  material, or an exact rational argument.
- **CERTIFIED**: a finite arithmetic enclosure, exact rational or interval, decided on endpoints.
- **CERTIFIED (padded binary64)**: a floating-point bound with an explicit forward-error bound, evaluated in
  IEEE 754 binary64 with round-to-nearest under stated accuracy assumptions.
- **NUMERICAL**: a computed value without such an enclosure.
- **OPEN**: a statement not established.

Numbering in the built PDF: `sec:certified` is Section 4.4, the fourth subsection of the results section, "A
certified crossing at two discretisation points" (page 14), `def:cert` is Definition 4.6, `prop:sym` is
Proposition 4.7, `thm:C` is Theorem 4.8, and the reproduction subsection is "Theorem 4.8" in Section 7. The record
calls the theorem "Theorem C"; the manuscript numbers it and refers to it only by `\ref`.

## 0. What this version adds to the map

The map of versions 1 and 2 is unchanged for every statement that this version does not touch: rows C1 to C8 of the
version 2 map and its open-problem rows stay valid, because the corresponding paragraphs, labels and record
references are byte-identical in `paper/main.tex`. The only exception is that the `main.tex` line numbers recorded
there are shifted by the insertions listed in `BUILD.md` section 4; identify statements by label, not by line.
Nothing is withdrawn: every statement that version 2 asserted is asserted here with the same label.

## 1. New statements of version 3

| claim | label | record and evidence | where in `main.tex` |
|---|---|---|---|
| **D1** The finite two-flat-sheet model: moments, dispersion, $A=1$, $N_0=1/2$, $\omega_n=2n+1$ for $n=-N_w,\dots,N_w-1$, $M=2N_w$, $\tau=\eta\Cinf\sqrt{\lambda}m$, $K(d)=m^2/(m^2+(2\pi\tau d)^2)$, $g_{\rm RB}=2\lambda\tau$, $C_X=(2\lambda\tau/N_k)^2$ as exact rationals, the map $\Phi$, and the linearised operator $B$ in the az gauge. | definition (no label of proof) | R7. `CERT.md:84-118` (section 2), a verbatim transcription of `solve_finite_m_lattice`; fidelity checks `CERT.md:283-292` ($\Phi$ vs the FFT map $8.89\cdot10^{-17}$/$7.51\cdot10^{-17}$, $B$ vs `apply_full` $6.66\cdot10^{-16}$/$4.44\cdot10^{-16}$, both \NUM) | `def:cert`, lines 404 to 429 |
| **P1** $B(\Sigma)=B(\Sigma)^{\mathsf T}$ for every $\Sigma\in\mathbb C^{M\times N_k}$, with no symmetry assumed on $\Sigma$; consequently a Rayleigh quotient of a fixed vector bounds $\lambda_{\max}(B)$, and the earlier nonsymmetric perturbation inequality is not used. | PROVED | R7. `CERT.md:119-140` (Proposition 1, proof); sanity check (\NUM) `CERT.md:138-140`, $5.55\cdot10^{-17}$ on non-inversion-even perturbations | `prop:sym` with its proof, lines 431 to 445 |
| **C9a (N)** Fixed point. $\Phi$ has exactly one fixed point in the sup-norm ball of radius $\rho=10^{-12}$ around $\Sigma_0$, with $\|\Sigma^\star-\Sigma_0\|_\infty\le5.83329104641\cdot10^{-17}$ ($m=1/100$) and $\le3.63702630192\cdot10^{-17}$ ($m=1/200$), from $r\le5.83329104451\cdot10^{-17}$, $k\le3.25923499198\cdot10^{-10}$ and $r\le3.63702630028\cdot10^{-17}$, $k\le4.50507433578\cdot10^{-10}$ with $r+k\rho\le\rho$, $k<1$. No contraction of $\Phi$ is used or asserted. | CERTIFIED | R7. `CERT.md:39-47` (statement), `:223-234` and `:238-248` (bounds tables), `cert_main.log`, `cert_m200.log`, `cert_default.log` | `thm:C` item (N), `eq:certN`, lines 452 to 466 |
| **C9b (K)** Crossing on the whole box: $\lambda_{\max}(B(\Sigma))\ge1.00579751250$ at $m=1/100$ and $\ge1.01993371122$ at $m=1/200$ for every $\Sigma\in\mathrm{Box}(\Sigma_0,\rho)$, i.e. $\ge1+\theta$ with $\theta=5\cdot10^{-3}$, clearing $1+\theta$ by at least $7.97512500371\cdot10^{-4}$ and $1.49337112263\cdot10^{-2}$, hence exceeding unity by at least $5.7975125\cdot10^{-3}$ and $1.993371122\cdot10^{-2}$. | CERTIFIED, by the symmetric Rayleigh quotient on the box with P1 | R7. `CERT.md:49-55` (statement), `:168-178` (method), `:179-182` (precision replay); margins are `bound - 1 - theta` as printed by `cert.py:784-785` | `thm:C` item (K), `eq:certK`, lines 457 to 466 |
| **C9c (S)** Second checker: an independently written program with its own transcription and its own 160-bit MPFR directed-rounding intervals reproduces (N) for both parameter sets ($r\le1.43466873225\cdot10^{-16}$, $k\le1.53363171691\cdot10^{-9}$; $r\le1.00667718670\cdot10^{-16}$, $k\le2.06306136514\cdot10^{-9}$) and (K) at $N_k=8$, $m=1/100$ ($\ge1.03881800695$) and at $N_k=24$, $m=1/100$ ($\ge1.00579751243$). At $N_k=24$, $m=1/200$ it also reproduces (K), with $\lambda_{\max}\ge1.01993371112$. | CERTIFIED for every reproduced item of both models, each on two interval programs | R7. `CERT.md:65-78`, `independent_nk8.log`, `independent_nk24.log`, `independent_nk24_m200_rerun.log` (`[J] PASS`, `exit=0`, wall 2020.8 s), `independent_nk24_m200.log` (the killed first attempt, stage `[H]`, kept as the record of it: `CERT.md:71-74` and the P1-2 row at `CERT.md:388`) | `thm:C` item (S), line 463 |
| **C9d (V)** Vertex-decisiveness: the pure-rainbow column of the same lattices gives $\lambda_{\rm ME}=0.980079149$ ($m=1/100$) and $0.994520242$ ($m=1/200$), both below unity, so the crossing of C9b is produced by the crossed skeleton and not by the rainbow equation. No certified upper bound on $\lambda_{\rm ME}$ is computed. | NUMERICAL, explicitly not certified | R7. `CERT.md:57-63`, `RESULT.md` section 2 (marked NUMERICAL), values from `records/R6-flat-sheet/verify_1d_counterexample.py::solve_finite_m_lattice`; `CERT.md:283-292` last column | `thm:C` item (V), line 464 |
| **C9e** Guard and refusal record: box guard (radii $\ge\rho/2$), constant guard (no radius-0 certified coupling), endpoint decisions everywhere, one guard for the pole-free strip ($\min|\mathrm{Im}\,(i\omega_n\nu-\xi-\Sigma)|\ge0.0565776818037$ and $\ge0.0384152648904$), ten broken variants refused with exit 1, and both bounds reproduced at 400 bits. | CERTIFIED (guards); the pole distances are certified, the mutant and precision replays are exit-code evidence | R7. `CERT.md:142-193` (method), `:194-218` (mutants), `:219-282` (all bounds), `:320-362` (commands and exit codes), `mutants/*.log`, `cert_prec400.log`, `cert_prec400_m200.log` | lines 470 to 494, 703, 710 |
| **C9f** Adjacent discretisations ($N_k=48$, $N_k=64$, $N_w=16$, $\lambda=6$, $m=2/100$, $m=1/400$): binary64 diagnostics only, not enclosures; the crossing survives $N_k\to48,64$ and $N_w\to16$; $m=1/400$ does not isolate the vertex effect because the rainbow column alone crosses there. | NUMERICAL | R7. `CERT.md:296-318` (both tables), `ladder.log` | line 492 |
| **C9g** Scope: Theorem C is a statement about the finite discretised model; it is not the continuum operator, carries no Matsubara tail, is not a measured $T_c$ and not a material, and does not address the charge and Peierls instabilities of the nested band. | statement of scope | R7. `CERT.md:7-8`, `:296-318`; band caveat unchanged from version 2 (`rem:cdw`) | lines 494 to 498 |

## 2. Changed status of existing statements

| statement | version 2 | version 3 | evidence |
|---|---|---|---|
| `prob:flatcert`, Gap 1 (fixed-point defect) | OPEN at every truncation, with the phase-blind supremum norm of the Jacobian at about $2.96$ blocking a contracted remainder bound | closed at $N_k=24$, $N_w=12$, $m=1/100$ and $m=1/200$, where the inversion of $I-DF$ is in $576$ dimensions (C9a); still OPEN at the truncations of `prop:flatcert` (dimension of order $10^5$) | `CERT.md:39-47`, `:296-318`; `main.tex` lines 398, 765 |
| `prob:flatcert`, Gap 2 (Matsubara tail) | OPEN | OPEN, unchanged | `CERT.md:302` |
| `prob:flatcert`, Gap 3 (continuum and adiabatic limits) | OPEN | OPEN, unchanged | `CERT.md:303` |
| `prop:flatcert` (padded binary64 lower bounds) | CERTIFIED (padded binary64) at $(\lambda,m)=(6,0.008)$, $N_k=2050$ and $(8,0.005)$, $N_k=3074$ | unchanged label and numbers; keeps its scope | `main.tex` lines 386 to 394 |
| `prop:flatcross` (crossing of the scaling limit) | NUMERICAL | NUMERICAL, unchanged | `main.tex` lines 364 to 384 |

## 3. Carried forward

Every other statement of versions 1 and 2 keeps its label, its record and its wording, including rows C1 to C8 of
the version 2 map: the flat-sheet scaling limit (PROVED, R6), the numerical crossing of that limit (NUMERICAL, R6),
the padded binary64 bounds (CERTIFIED (padded binary64), R6), the CDW caveat (NUMERICAL scope restriction, R6),
coarea sufficiency and the layer-cake reduction (PROVED with CERTIFIED constants, R6), the all-$x$ Young-ball
closure and the small-$x$ curved-surface window (PROVED, R6), and the supersession note. Record labels: R0 to R6 as
in version 2, R7 the finite-model certificate of this version.

## 4. Claims the manuscript does not make

- **(V) is not a certified separation.** No certified upper bound on $\lambda_{\rm ME}$ exists; the comparison of
  the full kernel with the rainbow column is a binary64 `eigsh` comparison on the same lattices.
- **No continuum claim.** The certificate is at fixed $N_k=24$, $N_w=12$; the $N_k=48,64$ and $N_w=16$ values are
  binary64 diagnostics and the Matsubara tail has no enclosure.
- **No material claim.** No measured $T_c$, no material prediction, and no statement that the vertex correction
  lifts a ceiling for any physical system. The band is nested, and the static charge response indicates that the
  chain and Peierls channels compete.
- **No theorem beyond the record.** The number of checkers, the refusal of the ten variants, the 400-bit replay
  and the two checker programs of the $m=1/200$ crossing stage are reported exactly as they stand in `CERT.md`.

## 5. Policy checks of this version

- Byline: Davis Rens, SuperMatics Inc., ORCID 0000-0001-5839-5197. Unchanged from version 2.
- "This preprint has not been peer reviewed." appears exactly once in `paper/main.tex`, in the section "Supporting
  material and data availability", which also carries the unchanged version 2 sentence about AI-generated
  proofs, code and manuscript drafts.
- No em dash anywhere in `main.tex` or in the rendered text; no internal project or process name.
- The released manuscript still carries the version 2 record identifiers in its DOI macros (`\MSPAPERDOI`
  10.5281/zenodo.24410274, `\MSSUPPDOI` 10.5281/zenodo.24410273), because the deposited PDF is the reviewed
  text and was not rebuilt for the release. The records of this version are 10.5281/zenodo.24415577 (paper) and
  10.5281/zenodo.24415572 (supplement); the difference is metadata only and is recorded in the release ledger.
  `README.md`, `CITATION.cff` and this file carry the version 3 identifiers.
- No placeholder marker remains in `main.tex`: the $m=1/200$ second-checker run passed (`[J] PASS`, `exit=0`) and
  its result is written into item (S) and into the `Version 3 note` paragraph.
- Every hash in this file is of a released file: `BUILD.md` sections 1 and 3 list the bytes, and the manifests
  verify the tree in CI.
