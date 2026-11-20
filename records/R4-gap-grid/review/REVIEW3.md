# G1-GAP phase 3 adversarial review

**Overall: ACCEPT-WITH-FIXES.** The analytic bounds and the pointwise certification criterion hold under the stated binary64 trust base. Qualify the verifier independence and the norm being bounded, and bind the receipt check to the published claims. There is no universal ceiling proof here.

| Claim | Verdict | Evidence and scope |
|---|---|---|
| Static-transfer equations and vertex | **ACCEPT** | The transfer profile converges in ℓ¹ to its zero-transfer part. The two ends give `2(r*r−Q*Q)` and the middle gives `r*r+Q*Q` after inversion, so the negative term survives. The Ward criterion is conditional. |
| `x≤2⁻⁴⁰`, gap `>0.00807313` | **ACCEPT** | The reviewed G1-SIGNED absolute estimate, uniform transfer bound and signed mass cancellation give a positive margin on the whole stated coupling and η ranges, throughout the homotopy. Independent rational gap ≥0.008080278229. No G1-ONSET signed estimate is imported. |
| `x≥19/2`, spectral top `<0.70468692` | **ACCEPT** | Independent unfolded bound ≤0.704684220003; exact triple-overlap volume `15/32`. The same ball works for every coupling scale. |
| Peak rectangles and phase 2 pointwise grid | **ACCEPT** | All ten rectangle anchors pass the peak guard; each positive `A k_j` and its products decrease in `x,η` and increase in `λ`. This extends the majorant, not the actual eigenvalue. Peak counts 15/45 and 16/45, plus 7 and 8 pointwise phase 2 closures. The earlier inverse-moment pre-quantization fix is present. |
| Phase 3 pointwise nodes and bounds | **ACCEPT-WITH-FIXES** | Exact receipt union matches the published 28 nodes (12 and 16), giving 34/45 and 40/45. First qualifying receipt maxima: 0.72827964 for pairing and 0.70929592 for the **preconditioned** normal map. The two paths are not end-to-end independent. |
| Static 45/45 and near-static transfer | **ACCEPT** | Both full static replays previously closed 45/45 at ≤0.70431954. The omitted mass is at most `δ=1/(12 C∞² λη²)` and the crossed omitted mass is ≤`2δ+δ²`. The uniform product-tube budget covers every Matsubara block and works for all `η≥2¹⁸` at the same 45 points; no interpolation. |
| Weak adiabatic edge | **ACCEPT** | Hölder gives the full-space `L^(3/2)` triple bound; the DOS secant gives `(J/π)²<3`. At `R≤3/25` the image is <0.29737610 and the common derivative/Schur bound is <0.78467306. At `λ≤7/60`, independently bounded `R<0.119439824` for `η≥30`, and `R<0.118167163` for `x≥2⁻⁴⁰, Ω/W≤2⁻⁴³`. Below that `x` the accepted small-`x` theorem closes the rest of `η≥1`. In the static limit `R=λ`, so `λ≤3/25` works at every `x>0`. |

## Phase 3 checks

- **Regularity.** At each slab the head correction is finite rank. Its independently bounded defect makes `A(I−JV)` invertible; `A` is then invertible because it is surjective and identity plus finite rank. The derivative bound for `G=S+A(F−S)`, plus the endpoint residual/curvature bound, gives a unique root within each strict tube and an invertible true Jacobian. Junctions share the exact center and hence nested balls. The first center is bare; positive widths and rotation symmetry are preserved.
- **Momentum.** The odd-height convolution integrates the continuous radial separation. Cell range and norm, signed insertion grouping, and the positive convolution against a squared line bound the unprojected absolute row. The final Schur weight is `|Z|`; every allowed angular sector is covered, not just ℓ=0. The pairing middle insertion reflects the coefficient line; the normal one does not.
- **Frequency.** Low resolved rows pay for both disjoint wings of the double-transfer square. The larger moment panel also pays its outer complement. High outputs split on electron width and bound all remaining rows using the width cap and forced transfer separation. The finite-rank rainbow subtraction separately pays the omitted head columns and all tail rows. The first bare-tail row bounds the unrepresented norm maxima. Negative frequencies are included by reflection.
- **Arithmetic.** The represented dyadic slab endpoints/centers are the inputs. Coupling endpoints are widened, the two paths check inverse defects rather than trusting a solve, and the signed products and residuals use pre-cancellation magnitude pads. The integer path charges quantization and conversion, and checks accumulator capacity. This is padded binary64 under the stated rounding assumptions, not a formal interval-library or proof-assistant certificate.
- **Independence.** The second path genuinely recomputes the folded convolution, DOS moments, signed rows, inverse defect, and transfer sums on another mesh/panel. It reuses `Replay.finish`, `Replay.pairing`, and the moment/averaged guard implementations; agreement does not detect a shared error in them. The three mutant fixtures exercise the second convolution, omitted wing, and inverse-defect code.
- **Receipt provenance.** `refined-weak512.log:2` is a valid two-verifier closure printed before cancellation; `PRIMARY` and `UNRESOLVED` records were not counted. A fresh easy node gave slightly tighter bounds than its historical row. The reported maxima are from the historical receipt union, not from a fresh 28-node rerun.

## Exact fixes

1. **`RESULT.md:25`, `RIGOR.md:128`.** Replace the unqualified independent-replay description with “second spatial/frequency implementation with shared slab, pairing-growth and optional moment bounds.” Cite the shared paths at `radial.py:318,374`, `radial_check.py:254–274`, `radial_moments.py:192–205`, and `radial_averaged.py:9–43`. Do not present this as two independent proofs.
2. **`RESULT.md:25`, `radial.py:467`, `grid_receipts.py:79`.** Label `q` as the contraction bound for the preconditioned normal map `G`; it is not an upper bound on `||DF||`. The invertibility conclusion for the true normal Jacobian is unchanged.
3. **`grid_receipts.py:55–79`.** For published/default receipts assert the exact 28-node set and both published upper bounds before printing a passing summary; `len(added)>=28` plus per-row `<1` accepts a substituted node or a worsened published maximum. Keep optional new receipts separate from the pinned baseline. `phase3_checks.py` independently checks the exact current sets and bounds.

## Numerical trend and scope

At `(λ,η)=(2,1)` the sampled pairing margins shrink from about 0.117 to 0.095 as `x` falls from 0.04 to 0.01. This is consistent with a conservative lower bound of 0.00807313 inside the much smaller proved window. The samples do not establish closure, and vanishing margin is ruled out inside that window. Neither the trend nor the unresolved tube at `(2,1,1/16)` decides the intervening region. The 28 added nodes are pointwise in physical parameters; the coupling homotopy is not an interpolation in temperature.

## Replay receipts

- Fresh `radial_check.py`: **0**, 3/3 mutants rejected. Fresh `radial_moments.py --lam 1/9 --eta 1 --xs 1/8 --cells 128 --head 24 --steps 30 --panel 8`: **0**, `CLOSED`, pairing ≤0.27391491, preconditioned `q≤0.09038406`.
- Fresh `grid_receipts.py D:/t/G1Gap`, `replay.py`, `independent_checks.py`, and `phase3_checks.py`: each **0**. Prior review ran both full static modes with exit **0**; read those receipts, did not rerun them in phase 3.
- Not rerun: the remaining radial nodes, full phase 2 finite grid, spectral sweep, Lean, builds, lint, or formatters. No source changes outside `G1-GAP-REVIEW/`.
