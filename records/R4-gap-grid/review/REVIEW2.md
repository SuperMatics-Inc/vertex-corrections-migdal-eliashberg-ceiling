# G1-GAP phase 2 adversarial review

**Verdict: ACCEPT-WITH-FIXES.** The 15 added finite nodes and the static/near-static grid reproduce. The independent finite checker has one inward-rounding path before integer quantization. An exact audit shows it does not affect the published witnesses, and an in-memory repair preserves all 15 and the bound.

## Earlier fixes

1. **Applied, `replay.py:30–34,74–81`.** Rational bounds are printed by exact decimal ceiling/floor, in the right directions, and the replay asserts the published small-x gap.
2. **Applied, `replay.py:240–269`.** Each of the ten rectangle anchors must pass `choose` before the DOS fallback; the log records guard and h. Per-slice counts and the large-x target are asserted. Monotonicity concerns the positive `A k_j` factors and the guard, not the true eigenvalue.
3. **Applied, `static.py:249–250,273–287`.** Duplicate inputs are rejected. Every requested point must reach `FULL`, and the reported grid must have 45 points and meet the claimed target. The specified `--steps 1` run emitted `UNRESOLVED`, exited 1, and emitted neither `FULL` nor a pass summary.
4. **Applied, `static.py:126–134,145–156`, `RIGOR.md:61–63`.** The width is measured from the represented coupling midpoint, and separate signed arithmetic pads use the pre-cancellation scales of the residual and vertex. The analytic convolution and resolvent errors are added separately under the explicit binary64 trust base.

## Added finite nodes

- **Tube and branch: ACCEPT.** The exact center rainbow uses the folded difference; the crossed defect is bounded in absolute value. Both paths check the curvature of the affine-center defect between coupling endpoints, not merely the endpoint residuals. The tube includes arbitrary allowed momentum dependence. Direct and preconditioned contractions have positive widths. The preconditioner distinguishes the full rainbow derivative norm from its norm on momentum-constant functions, bounds its inverse by the squared-map argument, and checks the error on the full complex tube. Shared endpoint centers give nested balls and identify the regular branch from zero.
- **Tail and spectrum: ACCEPT.** Head rows pay for omitted transfer mass; high rows split on the electron index and use a bare-width tail. Both low-row panels index the three internal frequencies separately and cover both wings outside the transfer square. Cauchy-Schwarz/peak envelopes cover all allowed momentum sectors; the two ends and the middle insertion are included. The final pairing row is an absolute Schur bound, not an eigenvalue of a truncated matrix.
- **Second verifier: ACCEPT-WITH-FIXES.** It shares the analytic inequalities, scalar interval normalization and proposed tube witness. It independently recomputes slabs with folded integer dots, different cutoffs/panels and a different arctangent implementation. It is an independent implementation check, not a second analytic proof. Exact ceiling comparisons on 5,031 inverse-moment inputs found no inward quantization on the published witnesses.
- **Mutation evidence.** The supplied reflection, cutoff-tail and middle-insertion fixtures are meaningful local checks: the unfolded row distinguishes the sign and omitted mass, and exact cubic differentiation distinguishes dropping one insertion. They compare explicit wrong formulas, not full runs of mutated verifiers. `phase2_checks.py` executed two in-memory source mutants. Removing the interpolation allowance gave a bound of 0.0000328521 for a midpoint defect above 0.06925756. Removing the junction check let a discontinuous center witness pass that the original rejected.

## Static and near-static grid: ACCEPT

For each fixed grid point, the exact missed transfer fraction is at most `δ = 1/(12 C∞² λ η²)`. The crossed mass omitted is bounded by `2δ+δ²`, including transfers in either slot. With `|z|≤2/x` throughout the product tube, these yield the implemented uniform residual and derivative bounds (`static.py:230–237`). The single budget applies to every slab, all frequencies and both normal and pairing derivative estimates. The static tail bounds decrease with odd Matsubara width and include the reflected half-line. Componentwise endpoint distances can use the static-row derivative plus the *global* off-transfer image bound, even though the finite model couples the blocks. The exact independent maxima are below `5.726332e−8` for residual/x and `3.406204e−7` for each derivative/Schur remainder. Both full author receipts report 45/45 and `0.70431954`; the fresh worst-point replay agrees. The region is the specified 45 `(λ,x)` points and all `η≥2^18`, not interpolation or continuation toward `η=1`.

## Exact fix list and scope

1. **`finite_check.py:202–203`, `RIGOR.md:87–89`:** Use a downward-enclosed center width and upward-enclose `weight/(v²+x⁻²)` before passing it to the integer dot. Ceil after a nearest-rounded quotient is not itself outward: for `w=1`, `x=1/8`, `d=1.9998346544058268` the scaled exact quotient exceeds 57457 by about `5.69e−13`, but the current ceiling is 57457. The outer sum pad may absorb the loss, but the per-input outward claim is false for valid represented widths. The proposed padding passed all 15 existing witnesses in memory at the same published spectral target.
2. **`RESULT.md:34`, `RIGOR.md:91`:** Say the three supplied formula mutants are rejected by local fixtures, rather than implying three complete mutated verifier runs. The new review script exercises whole checker methods in memory.

The `PROVED` labels are appropriate for the static map formulas, continuum windows and sufficient tube criterion. The `CERTIFIED` labels attach to their respective regions only, with the finite rounding fix above. The transform results remain `NUMERICAL`. No sentence implies a universal ceiling or interpolation of the added nodes; the monotone peak rectangles and uniform near-static η range are distinct proved extensions.

## Replay receipts

- Fresh local: `replay.py` exit **0**; `finite.py --published` exit **0**, 15/59, worst outward 0.50417617; `static.py --lam 2 --xs 1/16 --steps 1` exit **1** without a pass; `static.py --lam 2 --xs 1/16 --steps 80 --independent` exit **0**, 1/1, outward 0.70431954.
- Review scripts: `phase2_checks.py` exit **0**; `independent_checks.py` exit **0**. The former checked the two mutants, the exact quantization comparisons and the in-memory repair; the latter recomputed the analytic constants and near-static budgets.
- Not rerun in this review: both full static modes (read the phase 2 author summaries and exit receipt), spectral sweeps, Lean, builds, lint, formatters. No source edits outside the review directory.
