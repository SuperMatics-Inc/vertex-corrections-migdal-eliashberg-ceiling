# Vertex corrections and the Migdal-Eliashberg ceiling: a geometric dichotomy between nested and curved Fermi surfaces

Davis Rens  
SuperMatics Inc.  
[ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)  
Version 3.0.0

Paper DOI: <https://doi.org/10.5281/zenodo.24415577>  
Supplementary records DOI: <https://doi.org/10.5281/zenodo.24415572>  
Previous version: paper <https://doi.org/10.5281/zenodo.24410274>, records <https://doi.org/10.5281/zenodo.24410273>  
First version: <https://doi.org/10.5281/zenodo.24407163>  
Source: <https://github.com/SuperMatics-Inc/vertex-corrections-migdal-eliashberg-ceiling>

This repository accompanies the preprint *Vertex corrections and the Migdal-Eliashberg ceiling: a geometric dichotomy between nested and curved Fermi surfaces*. The manuscript source and bibliography are in `paper/`. The statement-to-source map is `CLAIMS.md`, and the build record with the hashes of the released files is `BUILD.md`. The proof, replay and certificate records are grouped in `records/` as R0 to R7, matching the labels defined in the paper. `MANIFEST.sha256` pins every file in this repository except the manifest itself. `ADMITTED.sha256` records the accepted bytes of the manuscript, the map and the records as prepared for release.

The manuscript PDF is the reviewed text, so its DOI macros still carry the version 2 identifiers. The records of this version are the two DOIs above, and `CITATION.cff` names them.

## Mathematical scope

For the phonon-only scalar Migdal-Eliashberg (ME) equations on the positive reflected normal branches specified in the paper, the linear onset obeys the sharp ceiling `k_B T_lin < C_inf sqrt(A M)`, where `A` is the relative height of the electronic density of states, `M` is the phonon spectral moment, and `C_inf` rounds to `0.18272624777232952` (R0). The paper asks whether that ceiling survives the first crossed skeleton added to the rainbow self-energy, and the answer of this version is a geometric dichotomy rather than a single uniform bound.

- **PROVED (R6, flat-sheet limit).** On a one-dimensional band, and on the parallel flat sheets of a two-dimensional band, the first crossed skeleton has an `O(1)` limit as `Omega_max/W -> 0` at fixed coupling. Three of the four sheet combinations satisfy linear energy conservation and the fourth is separated by the full bandwidth, so the two crossed channels reduce to the exact Cauchy residues `I_1` and `I_23` and the correction does not vanish as `m -> 0` (`records/R6-flat-sheet/`).
- **NUMERICAL, with Rayleigh lower bounds CERTIFIED (padded binary64) at two fixed parameters (R6).** The linearized pairing operator of that limit is self-adjoint, and after Richardson extrapolation its top eigenvalue is `1.004245`, `1.013678`, `1.014628` at `lambda = 4, 6, 8`, so the linear onset lies `0.9%` to `6.5%` above the ceiling. The crossing itself is numerical: there is no enclosure of the continuum operator. At fixed discretization and fixed numerical normal solution, the Rayleigh quotients bound the matrix of the computed branch from below by `1.0090706553` at `(lambda, m) = (6, 0.008)` and `1.0123025473` at `(8, 0.005)`, with an explicit forward-error bound and a 50-digit interval cross-check.
- **CERTIFIED at two discretization points (R7, new in this version).** At `N_k = 24`, `N_w = 12`, `lambda = 8`, `eta = 1`, `mu_0 = -2/5` and phonon-mass ratio `m = 1/100` or `m = 1/200`, ball arithmetic proves that the normal-state map of the finite two-flat-sheet model has exactly one fixed point in the sup-norm ball of radius `10^-12` around the numerical solution, with defect at most `5.83329104641e-17` or `3.63702630192e-17` and contraction constant at most `3.25923499198e-10` or `4.50507433578e-10`, and that the linearized gap operator is symmetric and exceeds unity over the whole enclosing box, by at least `5.7975125e-3` or `1.993371122e-2`. Ten deliberately broken variants of the certificate exit 1, both bounds are reproduced at 400-bit precision, and a separately written second checker with its own transcription and its own 160-bit interval arithmetic reproduces the fixed point and the crossing for both models (`records/R7-flat-sheet-certificate/`, Theorem 4.8 of the paper). The statement concerns the finite discretized equations: it carries no Matsubara tail, no continuum or adiabatic limit, and no material claim.
- **PROVED (R6, curved Fermi surfaces).** For every `d >= 2`, strict convexity of the energy collars, a uniform lower bound on the Fermi velocity and nondegenerate pocket tangencies give a three-line strip overlap bound in which a coarea estimate replaces the spherical joint density. A layer-cake reduction turns it into a relative crossed correction of order `x (1 + log_+(1/x))^4` in `d = 2` and `x (1 + log_+(1/x))^3` in `d >= 3`, so the crossed form is a relative `o(1)` perturbation of the ME form in the adiabatic window. For a constant transfer `V_j = K_j/N_0` and an even density of states the ceiling itself holds up to that factor (Lemma 1, Lemma 2, Theorem 2 of the record).
- The results of versions 1 and 2 remain: an explicit inflation for a local flat band (R1), stability above an inflated threshold in the spherical window (R3), an onset bound strictly below the ceiling for `A lambda <= 1/26^2` (R3), continuum estimates and finite arithmetic certificates for the spherical Einstein model at `1/9 <= lambda <= 2` (R4), and three grid-refined numerical pairing roots in a finite lattice (R5).

The band of the flat-sheet construction is nested, and its static `2 k_F` charge response is strongly enhanced (`rpa_2kF` between `8.41` and `10.37` over the computed range), so the normal branch is not the physical ground state of that band and the flat-sheet results are statements about the truncated equations.

Every result is a statement about the stated model. Bounds are on the linear onset `T_lin` (or, in the lattice record, on a pairing onset `T_pair`), never on a measured critical temperature, and no material is claimed to be a superconductor. The finite certificates are exact rational, interval, ball or padded binary64 calculations under the elementary-operation accuracy assumptions stated in the records; the analytic arguments are proved on paper.

What remains open in this version: intermediate `x` at `lambda >= 1/9` without spherical symmetry, existence of the normal branch without particle-hole symmetry, end-to-end certification of the flat-sheet crossing at the truncations of the padded-binary64 proposition (the fixed-point defect and the crossing are now certified at two discretization points only), the Matsubara tail, the continuum and `m -> 0` limits, and degenerate tangencies or inflections.

## Records

| Record | Directory | Contents |
|---|---|---|
| R0 | `records/R0-scalar-ceiling/` | Scalar ceiling replay: interval enclosures of the model critical temperatures from two independent certificate implementations, the certified endpoint `0.1827262477723295176168720`, and guarded mutants. |
| R1 | `records/R1-local-skeleton/` | Inflated comparison for the local crossed skeleton, exact rational algebra replay with mutants, and the Lean 4 core. |
| R2 | `records/R2-signed-form/` | Signed momentum-resolved crossed form, exact Nambu interpolation and forward-margin guards with mutants, and Lean 4 sources for the finite algebra. |
| R3 | `records/R3-onset/` | Spherical onset estimates: signed fold totals, variance absorption, the Young ball and overlap bounds, and guarded mutants. |
| R4 | `records/R4-gap-grid/` | Spherical Einstein gap and grid: static and radial fixtures, the finite-grid enclosures, the receipt audit of the 28-node momentum-dependent set, and recorded logs. |
| R5 | `records/R5-lattice/` | Source of the finite-DOS lattice calculation behind the three numerical pairing roots. This record is numerical, has no enclosure, and its multi-hour sweeps are not rerun here. |
| R6 | `records/R6-flat-sheet/` | Flat-sheet scaling limit and crossing, the finite-`m` Rayleigh lower bounds with their forward-error bound, the curved-surface coarea and layer-cake lemmas with Theorem 2, the CDW caveat, the replay of the CERTIFIED constants with four rejected formula mutants, and the captured log of the flat-sheet solver. |
| R7 | `records/R7-flat-sheet-certificate/` | Ball-arithmetic certificate for the finite two-flat-sheet model: the transcribed model and the symmetry proposition, the fixed-point and crossing enclosures for both parameter sets, ten refused variants, the 400-bit replays, the second checker with its own interval arithmetic and its captured runs, the ladder diagnostics, and their logs. |

Each `independent_checks.py` and `phase*_checks.py` under `records/*/review/` is an independently written replay that CI runs; those directories hold no review narrative.

## Replay

The replays of R0 to R6 use Python 3.12 with pinned dependencies: `numpy==2.5.5`, `scipy==1.18.1`, `sympy==1.14.0`, `mpmath==1.3.0`, `python-flint`. The certificate of R7 uses Python 3.13 with `mpmath==1.4.1`, `numpy==2.5.5`, `scipy==1.18.1`, `python-flint==0.9.0` and `gmpy2==2.3.2`. From the repository root:

```sh
python -m pip install numpy==2.5.5 scipy==1.18.1 sympy==1.14.0 mpmath==1.3.0 python-flint
(cd records/R0-scalar-ceiling && python -B replay.py)
(cd records/R1-local-skeleton && python -B check_exact.py)
(cd records/R2-signed-form && python -B replay_form.py && python -B replay_bounds.py)
(cd records/R2-signed-form/review && python -B independent_checks.py)
(cd records/R3-onset && python -B replay.py)
(cd records/R3-onset/review && python -B independent_checks.py)
(cd records/R4-gap-grid && python -B replay.py && python -B finite_check.py && python -B radial_check.py && python -B finite.py --published && python -B grid_receipts.py logs)
(cd records/R4-gap-grid/review && python -B independent_checks.py)
(cd records/R6-flat-sheet && python -B replay.py)
(cd records/R6-flat-sheet && for m in layer_cake_factor u_floor residue_sign no_decay; do python -B replay.py --mutate $m || echo "mutant $m refused"; done)
```

Each replay exits 0 when its fixtures pass and prints the guards it rejected. `records/R6-flat-sheet/replay.py` certifies the Lemma 2 and Theorem 2(1) constants and the exact flat-sheet overlap ratio `9/4` in exact rational and 60-digit directed-rounding interval arithmetic, cross-checks the bare crossed residue `I_1` and the circle strip-overlap decay numerically, and must exit nonzero on each of its four formula mutants.

The flat-sheet solver `records/R6-flat-sheet/verify_1d_counterexample.py` is the long computation of R6: a full run takes about 8 minutes and its captured output is `records/R6-flat-sheet/verify_1d.log`. CI runs the solver and, separately, checks the captured log by SHA-256.

The certificate of R7 is short: with the five dependencies above,

```sh
(cd records/R7-flat-sheet-certificate && python -B cert.py)
(cd records/R7-flat-sheet-certificate && python -B independent.py --selftest)
(cd records/R7-flat-sheet-certificate && python -B independent.py --selftest-gap)
```

The default run of `cert.py` decides both parameter sets in one process in about 77 s on the reference machine and prints the certified fixed-point radius, contraction constant and Rayleigh bound of each; each of its ten variants is a separate short run that must exit 1; the 400-bit replays double the precision of both parameter sets; and `independent.py --nk 24 --m 100` and `independent.py --nk 24 --m 200` are the long second-checker runs of the record, at about 2020 s and about 2100 s. The captured logs of all of these are in `records/R7-flat-sheet-certificate/`.

The replays check finite algebraic statements, exact rational identities, outward intervals and guarded mutants. They do not formalize the analytic arguments. The Lean 4 sources under `records/R1-local-skeleton/` and `records/R2-signed-form/lean/` cover the finite algebra only; the recorded toolchain is Lean 4.35 with Mathlib, and Lean is not installed in this repository's CI. The R5 lattice sources are included for completeness and are not run in CI, because their sweeps take hours.

CI (`.github/workflows/ci.yml`) verifies `MANIFEST.sha256` and `ADMITTED.sha256`, runs the replays above plus the four R6 mutants, runs the default R7 certificate and both self-tests of the R7 second checker, applies the public prose safeguards to the manuscript and the shipped prose, and builds the paper with Tectonic 0.17.0 at a fixed source date, failing on any overfull box, undefined reference or unresolved citation. The ten R7 variants, the 400-bit replays and the long R7 second-checker runs are not rerun in CI; their logs are shipped and pinned by hash. The certificate job reads the printed quantities numerically rather than as fixed strings: the certified Rayleigh bounds and branch-guard margins are compared against the manuscript's values to eleven digits, and the printed fixed-point radii against the manuscript's quoted bounds, because the radius is computed from the platform's floating-point candidate and differs in its leading digits between the reference machine (about `5.8e-17`) and the CI runner (about `4.8e-17`) while staying below the quoted bound on both.

## Changes in this version

- New Section 4.4, `A certified crossing at two discretisation points`, with Definition 4.6, Proposition 4.7 and Theorem 4.8, and the `Version 3 note` paragraph that states what is added.
- New record R7 with the certificate program, the second checker, the ten refused variants, the 400-bit replays, the ladder diagnostics and their logs.
- Gap 1 of the open flat-sheet problem is closed at `N_k = 24`, `N_w = 12` for `m = 1/100` and `m = 1/200`; Gaps 2 and 3 remain open.
- No statement of versions 1 and 2 is withdrawn, and `refs.bib` is unchanged.
- The released PDF is the reviewed manuscript; its DOI macros still name the version 2 records, and the version 3 identifiers appear in this file and in `CITATION.cff`.

## Disclosure and review

AI agents generated the proofs, code and manuscript drafts under D. Rens' direction. Separate AI agents performed the independent reviews. This preprint has not been peer reviewed.

## Rights

The manuscript text and bibliography in `paper/` are distributed under CC-BY-4.0 through the paper record. The replay code, records and other software materials in this repository are distributed under Apache-2.0. The paper's license does not change the software license.
