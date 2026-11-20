# Vertex corrections and the Migdal-Eliashberg ceiling: rigorous bounds in a spherical overlap model

Davis Rens  
SuperMatics Inc.  
[ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)  
Version 1.0.0

Paper DOI: <https://doi.org/10.5281/zenodo.24407163>  
Software supplement DOI: <https://doi.org/10.5281/zenodo.24407165>  
Source: <https://github.com/SuperMatics-Inc/vertex-corrections-migdal-eliashberg-ceiling>

This repository accompanies the preprint *Vertex corrections and the Migdal-Eliashberg ceiling: rigorous bounds in a spherical overlap model*. The manuscript source and bibliography are in `paper/`. The statement-to-source map is `CLAIMS.md`. The proof, replay, certificate and review records are grouped in `records/` as R0 to R5, matching the labels defined in the paper's introduction. `MANIFEST.sha256` pins every file in this repository except the manifest itself. `ADMITTED.sha256` records the accepted bytes of the manuscript, the map and the records as prepared for release.

## Mathematical scope

For the phonon-only scalar Migdal-Eliashberg (ME) equations on the positive reflected normal branches specified in the paper, the linear critical temperature obeys the sharp ceiling `k_B T_lin < C_inf sqrt(A M)`, where `A` is the relative height of the electronic density of states, `M` is the phonon spectral moment, and `C_inf` rounds to `0.18272624777232952`. The paper asks whether that ceiling survives the first crossed skeleton added to the rainbow self-energy.

For a local flat band it gives an explicit inflation and recovers the uninflated bound at fixed coupling for sufficiently small `Omega_max/W`. In a projected spherical band it proves stability above an inflated threshold within a small-`x` window. A pointwise contraction for `A R(T) <= 1/9` gives a high-temperature bound under `A lambda < 1/9` and an onset bound strictly below the ceiling under `A lambda <= 1/26^2`. An overlap argument gives the uninflated bound for bounded weak coupling under a sufficient bandwidth-ratio guard. For the spherical Einstein model at `1/9 <= lambda <= 2` the finite arithmetic certificates extend coverage but leave regions above the ceiling unclosed, and a separate lattice calculation gives three grid-refined numerical pairing roots below the ceiling.

Neither a universal bound for the crossed skeleton nor a counterexample is established. The results are mathematical statements about the stated model. They bound a linear or pairing onset temperature in that model, not a measured critical temperature, and they confirm no superconductor. The finite certificates are exact rational and padded binary64 interval calculations under the elementary-operation accuracy assumptions stated in the records; the analytic arguments are proved on paper.

## Records

| Record | Directory | Contents |
|---|---|---|
| R0 | `records/R0-scalar-ceiling/` | Scalar ceiling replay: interval enclosures of the model critical temperatures from two independent certificate implementations, the certified endpoint `0.1827262477723295176168720`, and guarded mutants. |
| R1 | `records/R1-local-skeleton/` | Inflated comparison for the local crossed skeleton, exact rational algebra replay with mutants, Lean 4 core, and the independent review. |
| R2 | `records/R2-signed-form/` | Signed momentum-resolved crossed form, exact Nambu interpolation and forward-margin guards with mutants, Lean 4 sources for the finite algebra, and the independent review. |
| R3 | `records/R3-onset/` | Spherical onset estimates: signed fold totals, variance absorption, the Young ball and overlap bounds, guarded mutants, and the independent review. |
| R4 | `records/R4-gap-grid/` | Spherical Einstein gap and grid: static and radial fixtures, the finite-grid enclosures, the receipt audit of the 28-node momentum-dependent set, recorded logs, and three review rounds. |
| R5 | `records/R5-lattice/` | Source of the finite-DOS lattice calculation behind the three numerical pairing roots. This record is numerical, has no enclosure and no independent review, and its multi-hour sweeps are not rerun here. |

`paper/`, `CLAIMS.md` and `records/` were prepared in the working layout referenced by `CLAIMS.md`, whose provenance comments in `refs.bib` and whose record proofs keep the source record names; the record table above maps each R label to its directory here. `reviews/manuscript/` holds the three review rounds of the manuscript itself.

## Replay

Python 3.12 is used with pinned dependencies: `numpy==2.5.5`, `scipy==1.18.1`, `sympy==1.14.0`, `mpmath==1.3.0`, `python-flint`. From the repository root:

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
```

Each replay exits 0 when its fixtures pass and prints the guards it rejected. The replays check finite algebraic statements, exact rational identities, outward intervals and guarded mutants. They do not formalize the analytic arguments. The Lean 4 sources under `records/R1-local-skeleton/` and `records/R2-signed-form/lean/` cover the finite algebra only; the recorded toolchain is Lean 4.35 with Mathlib, and Lean is not installed in this repository's CI. The R5 lattice sources are included for completeness and are not run in CI, because their sweeps take hours.

CI (`.github/workflows/ci.yml`) verifies `MANIFEST.sha256`, runs the replays above, applies the public prose safeguards to the manuscript and this README, and builds the paper with Tectonic, failing on any overfull box or undefined reference.

## Disclosure and review

AI agents generated the proofs, code, and manuscript drafts under D. Rens' direction. Separate AI agents performed the independent reviews. This preprint has not been peer reviewed.

## Rights

The manuscript text and bibliography in `paper/` are distributed under CC-BY-4.0 through the paper record. The replay code, records and other software materials in this repository are distributed under Apache-2.0. The paper's license does not change the software license.
