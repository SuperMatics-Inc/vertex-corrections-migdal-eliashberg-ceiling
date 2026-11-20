# G1-SIGNED finite Lean formalization

**Scope:** finite algebra and rational guards only. Source references below are numbered lines in `../G1-SIGNED/RIGOR.md`; line 12 of `../G1-SIGNED/RESULT.md` identifies the replay checks. The original result file predates this project.

## Formalized
- Lines 34, 42, 44, 46: `forward_guards`, `corner_window`, and `forward_geometry` cover the contraction constants, line bounds, witness boxes, transfer fraction, and exact margin `55/262144 > 1/8192`.
- Lines 40 and 67: `witness_guards`, `angular_guards`, and `angular_absorption` prove the rational witness-temperature, kernel, and absorption constants given the stated hypotheses; `angular_partial_fraction` proves the unequal-argument rational identity, conditional on nonzero denominators. The spectral decimal enclosure itself is an external input.
- Lines 79 and 116: `overlap_jacobian` proves the Jacobian algebra for nonzero arguments; `overlap_guards` proves the weak-coupling endpoint and contraction constants.
- Lines 17 and 26: `nambu_interpolation` evaluates the literal 2×2 Nambu product at ±1, ±2 and extracts its anomalous coefficient; `nambu_normal` extracts its normal imaginary part. `odd_interpolation` gives the general cubic interpolation identity.
- Lines 22–26: `four_anchor` proves the exact polynomial four-anchor average. Its six pair terms are the symmetric average of the three pair terms in equation (1); equating the integrated forms additionally requires the stated change-of-anchor symmetries.

## Not formalized
- The partial-fraction **integral**, angular measure, Jacobian change of variables, continuum reanchoring and infinite-sum identities require analytic integration or convergence arguments; the algebraic parts above do not certify them.
- The exact finite cutoff-10 ordered/folded angular sum, the finite-grid full/static replay totals and folded ME mean/variance replay are not separately enumerated in Lean. The general Nambu and four-anchor identities do not by themselves formalize their boundary and symmetry bookkeeping.
- The replay's mutant rejection tests and the certified external spectral enclosure are not re-established here; its decimal lower-bound source is an input from `../GAMMA/RESULT.md`.
- Branch existence, spectral claims, limits, onset exclusion, and the universal ceiling are not formalized. This is not a proof about exact many-body Tc.

Lean 4.35 with the shared prebuilt Mathlib overlay: `build.ps1` ran `lake build` and exited 0. The `#print axioms` output for every theorem contains only `propext`, `Classical.choice`, and `Quot.sound`. No replay or other checks were run in this project.
