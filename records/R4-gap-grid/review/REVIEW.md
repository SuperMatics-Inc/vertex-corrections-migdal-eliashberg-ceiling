# G1-GAP adversarial review

**Overall: ACCEPT-WITH-FIXES.** The continuum windows and finite Einstein stability regions follow in the stated model; this is not a universal ceiling. The static certificate needs a fail-closed replay and more explicit binary64 pad accounting.

| Claim | Verdict | Review |
|---|---|---|
| Static normal equation, vertex and Ward identity | **ACCEPT** | `g_j/g` converges in ℓ¹ to the zero-transfer point mass. Momentum inversion turns `b=a+c-p` into the triple convolution. The two end insertions give `2(r*r−Q*Q)` and the center gives `r*r+Q*Q` on the even pairing space. Taking minus the imaginary part gives the Ward identity; positivity is conditional, and map convergence is not a global branch theorem. |
| Small x, gap > 0.00807313 | **ACCEPT-WITH-FIXES** | The independent rational upper/lower computations certify a gap > 0.00808027822. Both normal contraction rows and the invariance/width checks are uniform on the coupling homotopy, by `R≤23/5` and monotonicity of the log-cubed envelope up to the endpoint. The proof uses the reviewed G1-SIGNED overlap/split and form bound, not G1-ONSET or the scalar strict-gap theorem. Fix the replay's decimal directions. |
| Large x, top < 0.70468692 | **ACCEPT** | An independently summed unfolded row at `(λ,η,x,h)=(2,1,19/2,9/25)` has image < 0.354373434329 and row ≤ 0.704684220003. The volume `15/32` also checks by interchanging the lens and radial integrals. |
| Finite Einstein 15/45 and 16/45, including the peak-line boxes | **ACCEPT** | Independently found first-guard witnesses at all ten row thresholds. For fixed `h`, every `Ak_j` is nonincreasing with `x` and `η` and nondecreasing with `λ`, and the positive sums/products preserve this order. No eigenvalue monotonicity or DOS fallback is needed for the rectangles. |
| Static bound 0.70431955 and transfer for η ≥ 2^18 | **ACCEPT-WITH-FIXES** | Both replays close 45/45 with a maximum printed 0.70431954 at (λ,x)=(2,1/16), consistent with the stated 0.70431955 bound. The analytic enclosure and transfer budgets match the code. Both spatial paths use the same scalar padding and tail logic, so agreement is not an independent audit of either. See the arithmetic and replay fixes below. |

## Tail, sectors and branch

- **Small x:** The Schur weight `w_n^(-1/2)` leaves `1/w_m` in each normalized row. Off-diagonals use the monotone Einstein envelope; the canceled diagonal keeps the negative mass. The weighted folded difference identity is `Σ_m P_nm w_m = R w_n`, by symmetric transfer cutoff; only O(1/N) boundary terms remain. The high rows and the low-row omitted tails are bounded separately. The mean/variance split includes all momentum dependence and the crossed form costs at most 10L. The width upper bound converts the form margin to a raw pairing gap.
- **Peak boxes:** The maximum convolution row lies at the nearest integer to the fermion reflection center −1/2; the two tails are bounded by an integral of `m^(-p-2)`. The DOS alternative uses decreasing analytic profiles and the same tails. The row bound integrates the absolute kernel over all momenta, with Schur weight `|z|`. For normal continuation, the image and derivative bounds hold at every `K→θK` in the same ball; the width preserves its sign.
- **Static boxes:** The odd height convolution, integration of its piecewise-linear derivative and the cell range bound control continuous separations. The derivative and absolute Schur envelopes act on the whole inversion-even momentum space, not just radial gaps. The bare tail uses `A/w²` past the first admitted positive block and reflects to negative frequencies. Each slab is a product ball with a uniform derivative bound below 1; its endpoint overlap can be proved componentwise from the static residual and derivative because the **whole off-transfer map**, not only its value at the center, is bounded by the residual budget. Both the normal and pairing derivative bounds include the omitted transfers. Independent maxima on the 45 points: residual/x ≤ 5.726332e−8, derivative remainder ≤ 3.406204e−7.

## Exact fix list

1. **`replay.py:70–73`**: Round printed rational bounds outward (ceil for L, floor for gap) and assert the actual published lower target. The current L print is below the exact proxy it formats, though the independent log-series bound confirms that it still exceeds the actual L. The mathematical gap is unaffected.
2. **`replay.py:233–251`**: Emit each anchor's guard kind and h; require that the ten rectangle anchors pass the peak guard. Assert the per-slice counts and the large-x decimal bound rather than just printing them. This binds a successful replay to the published regions, not just to an aggregate including static points.
3. **`static.py:260–269`**: Return failure if any requested grid point does not reach `FULL`. For the reported grid assert 45 unique completions and the upper target before emitting a passing summary. At present an `UNRESOLVED` block closes only the inner loop, so exit 0 alone is not a certificate.
4. **`RIGOR.md:61`, `static.py:127–129,140–150`**: Make the signed post-convolution roundoff accounting explicit: bound the coupling-midpoint and signed residual errors and the error in `3*rr−qq`/`1+g*corr` at their pre-cancellation scales, and identify the pad paying each error. The convolution/dot-product gamma bounds and widened coupling endpoints are suitable, but a nextafter after cancellation by itself does not establish direction. State the binary64/elementary-operation trust base with the certificate.

## Numerical trend and scope

At (λ,η)=(2,1), the 160-cell/12-frequency run gives tops 0.88286478 and 0.90518420 at x=0.04 and 0.01, respectively. The sampled margins shrink from about 0.117 to 0.095, much larger than the proved lower bound 0.00807313. The samples are truncated and x=0.01 is more than ten orders of magnitude above 2^−40. The trend is consistent with the bound and does not show closure; approach to zero is excluded inside the tiny proved window. The intervening finite-transfer region still requires a certificate. A failed sufficient inequality is not evidence of a crossing.

## Replays

- `replay.py`: exit 0, finite counts 15/45 and 16/45; five rejected mutants.
- `static.py --independent` with the full grid and `--steps 80`: exit 0, 45 `FULL`, maximum printed 0.70431954 at (λ,x)=(2,1/16); four rejected spatial mutants.
- Direct `static.py` with the same grid and `--steps 80`: exit 0, 45 `FULL`, same maximum printed 0.70431954 and four rejected spatial mutants.
- `independent_checks.py`: exit 0; exact rational π/log bounds, independent Matsubara sums, volume, rectangle anchors and nonstatic budgets. No imports from the submitted replay.
- Not run: spectral diagnostics (read the reported sector log), Lean, builds, lint or formatters.
