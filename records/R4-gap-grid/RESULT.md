# G1-GAP: stability certified on parts of the gap; universal ceiling open

**Answer:** No counterexample was found. The available bounds do not prove stability for every $T>T_{\rm ceil}$ with $1/9\le\lambda\le2$. In particular, a failed box check is not an eigenvalue above one or evidence that the regular normal branch ends.

Use $x=\pi k_BT/W$ and $\eta=T/T_{\rm ceil}$. All results apply to the flat-DOS spherical band with $A=1$, fixed chemical potential, and only the rainbow plus first crossed skeleton. See [RIGOR.md](RIGOR.md) for normalization, branch construction and bounds.

**PROVED, static-transfer limit.** At fixed $x,\lambda$, letting $\Omega/W\to0$ decouples the Matsubara blocks. With $W=1$ and $g=2\lambda x/\pi$, the exact normal equation is $\Sigma=g\int z+g^2z*z*z$. With $z=r-iQ$, the unsymmetrized static vertex is $g+g^2(3r*r-Q*Q)$. The negative $Q*Q$ term remains. The static Ward identity gives a positivity criterion but not unconditional stability. The limit is $\eta\to\infty$, not the slice $\eta=1$.

**PROVED with certified constants, continuum windows.** For $0<x\le2^{-40}$, the regular normal branch exists and the pairing gap is greater than $0.00807313$ for the full requested coupling interval and all $\eta\ge1$. This uses the absolute G1-SIGNED overlap estimate and retains the canceled normal mass; it does not rely on the under-review G1-ONSET signed proof. For $x\ge19/2$, the pairing spectral top is below $0.70468692$ on the same parameter range by a zero-coupling contraction. These windows are disjoint.

**PROVED, weak adiabatic edge.** The sharper Hölder guard gives pairing upper bound 0.78467306 whenever $\sum_jK_j\le3/25$. It closes all $x>0$ in the static model for $1/9\le\lambda\le3/25$. For the finite Einstein model with $1/9\le\lambda\le7/60$, it closes all $x$ at $\eta\ge30$, and all $T\ge T_{\rm ceil}$ if $\Omega/W\le2^{-43}$ by overlap with the small-$x$ proof.

**CERTIFIED, finite Einstein grid.** Sampled $x=1/16,1/8,1/4,1/2,1,2,4,8,16$, with full Matsubara tails, all allowed momentum sectors, and regularity along the zero-coupling continuation. Each peak-box threshold also covers larger $x,\eta$ and smaller $\lambda$ by monotonicity of the majorant, not of the eigenvalue. The signed-normal tube additions are pointwise:

| $\lambda$ | Peak, $\eta=1$ | Extra, $\eta=1$ | Peak, $\eta=3/2$ | Extra, $\eta=3/2$ |
|---:|---:|:---|---:|:---|
| $1/9$ | 1 | $1/4,1/2$ | 1 | $1/8,1/4,1/2$ |
| $1/4$ | 2 | 1 | 2 | 1 |
| $1/2$ | 4 | 2 | 4 | 2 |
| 1 | 8 | 4 | 8 | $2,4$ |
| 2 | 16 | $4,8$ | 8 | 4 |

Peak boxes certify 15/45 and 16/45. The finite-transfer replay attempts the other 59 nodes and closes 15 (7 and 8) with pairing upper bound 0.50417617 on those additions. Phase-2 subtotal: 22/45 at $\eta=1$, 24/45 at $\eta=3/2$, leaving 44 nodes for the momentum-dependent replay.

**CERTIFIED, momentum-dependent continuation:** A second spatial/frequency implementation with shared slab, pairing-growth and optional moment bounds ([RIGOR.md §7](RIGOR.md#7-finite-einstein-grid-momentum-dependent-normal-enclosures)) confirms 28 added nodes (12 at $\eta=1$, 16 at $\eta=3/2$), for 34/45 and 40/45 combined. The historical receipt union, using the first qualifying receipt per node, gives pairing upper bound 0.72827964 and contraction bound $q_G\le0.70929592$ for the preconditioned normal map $G$, not for the unpreconditioned map.

| $\lambda$ | Phase 3 added $x$ at $\eta=1$ | Phase 3 added $x$ at $\eta=3/2$ |
|---:|:---|:---|
| $1/9$ | $1/16,1/8$ | $1/16$ |
| $1/4$ | $1/8,1/4,1/2$ | $1/16,1/8,1/4,1/2$ |
| $1/2$ | $1/4,1/2,1$ | $1/16,1/8,1/4,1/2,1$ |
| 1 | $1/2,1,2$ | $1/8,1/4,1/2,1$ |
| 2 | $2$ | $1,2$ |

This is finite-grid evidence, not a continuum claim. At $(\lambda,\eta,x)=(2,1,1/16)$ the 1024-cell refinement still failed the sufficient normal tube test at slab 25/72: preconditioned center defect 0.0119764 and minimum tested derivative bound 0.467679. Neither this failure nor a pairing upper bound above one detects instability.

**CERTIFIED, static and near-static grid.** All 45 $(\lambda,x)$ points passed the signed radial enclosure, with a common spectral upper bound of 0.70431955 (worst at $\lambda=2,\ x=1/16$). Both spatial implementations checked the full zero-coupling continuation, all allowed momentum sectors, and the infinite Matsubara tail. The proof includes the full finite Einstein model at those same 45 $(\lambda,x)$ points for every $\eta\ge2^{18}$, by an explicit bound on the nonstatic transfers. It is not a uniform assertion between sampled $x$ or $\lambda$, nor a continuation down to $\eta=1$.

**NUMERICAL, independent sine-transform solve.** With 96 radial cells and 10 positive Matsubara frequencies, maxima on the sampled $\ell=0$ grid were 0.88337148 at $\eta=1$, 0.75984269 at $\eta=3/2$, and 0.64367885 in the static lowest block, all at $(\lambda,x)=(2,0.04)$. At $\lambda=2$, a 160-cell, 12-frequency diagnostic extended to $x=0.01$: there it found 0.90518420 at $\eta=1$, 0.86116053 at $\eta=9/8$, and 0.66096644 in the static lowest block. The largest sampled value in $\ell=2,4,6$ was 0.01059268. No observed crossing triggered counterexample certification. These are truncated, calculated values, not bounds; the smallest $x$ has the weakest spatial resolution.

## Replay and scope

From `G1-GAP/`, use `UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy --with mpmath --with scipy --with sympy python -B` followed by the script and arguments. Logs are in `D:/t/G1Gap/`.

- `replay.py`: exit 0 (`replay-weak.log`), exact fraction versus interval guards, independent geometry and derivative fixtures, five rejected mutants.
- `static.py --lam 1/9 1/4 1/2 1 2 --xs 1/16 1/8 1/4 1/2 1 2 4 8 16 --steps 80`, then the same with `--independent`: both exit 0 (`fix-static-primary.log`, `fix-static-independent.log`). Both checked 45 exact angular fixtures, rejected four spatial mutants, and gave an outward upper bound of 0.70431954.
- `static.py --lam 2 --xs 1/16 --steps 1`: deliberately under-resolved; exit 1 with `UNRESOLVED` and no pass summary (`fix-unresolved.log`). The four corrected replay statuses are in `fix-exit-codes.log`.
- `finite.py --published` and `finite_check.py`: both exit 0 (`phase2-finite-rounding.log`, `phase2-check-rounding.log`) after fixing the inverse-moment enclosure before integer rounding. All 59 targeted nodes attempted; 15 passed both verifiers, and local fixtures reject three formula mutants.
- `radial_moments.py --cells 128 --head 24 --steps 30 --panel 8`: exit 1 (`map128.log`), 19 explicit two-verifier `CLOSED` lines, 25 unresolved. The nonzero exit is intentional when the selected set is not completely certified.
- `radial_averaged.py --lam 1/9 1/4 1/2 --xs 1/16 1/8 1/4 --cells 128 --head 24 --steps 36 --panel 8`: exit 1 (`averaged-low.log`), 8/15, adding $(1/4,1,1/8)$.
- `radial_moments.py --lam 1 2 --eta 1 3/2 --xs 1/2 1 --cells 256 --head 32 --steps 48 --panel 12`: exit 1 (`mid-growth.log`), 4/8, adding $(2,3/2,1)$.
- `refined-low.log`: exit 1, 4/8 with three new closures; `refined-high.log`: exit 1, 2/8 with two new closures; `last-slab.log`: exit 0, 1/1 new. Exact commands in [RIGOR.md §9](RIGOR.md#9-parameter-grids-and-replay-receipts).
- `radial.py --lam 1/4 1/2 1 2 --xs 2 --cells 128 --head 16 --steps 36`: exit 0 (`radial-pass2.log`), 3/3. `radial_check.py`: exit 0 (`radial-mutants-final.log`), three mutants rejected after tightening the block norm.
- `radial_moments.py --probe --lam 2 --eta 1 --xs 1/16 --cells 1024 --head 48 --steps 72 --panel 16`: exit 1 (`hard1024.log`); primary-only probe and no certificate.
- `spectral.py` and `spectral.py --cells 96 --count 10 --xs .04 .0625 .125 .25 .5 1 2 4 8 16`: exit 0 (`spectral-status.log`, `spectral-refine.log`).
- `spectral.py --cells 160 --count 12 --xs .01 .02 .04 .125 --lam 2 --ratio 1 1.125 inf --sectors 0 2 4 6`: exit 0 (`spectral-sector.log`).
- `spectral.py --cells 256 --count 12 --xs .015625 .03125 .0625 .125 --lam 2 --ratio 1 1.01 1.0625 1.5 --sectors 0`: exit 0 (`spectral-ceiling.log`); maximum at sampled $T>T_{\rm ceil}$ was 0.89644635 at $x=1/64,\eta=1.01$.
- `grid_receipts.py D:/t/G1Gap`: exit 0 (`coverage.log`), exact union and open-node lists. Detailed refinement commands and statuses are in [RIGOR.md §9](RIGOR.md#9-parameter-grids-and-replay-receipts).

No git branch or commit. No Lean, build, lint, formatter, cloud, paper retrieval, or external actions. The remaining gap is a regular-branch enclosure and pairing bound near the ceiling outside the certified nodes and windows.
