# ACCEPT-WITH-FIXES

All six PARTIAL items from `REVIEW2.md` (items 7, 8, 12, 14, 15, 22) and all four delta-review V4 wording items are **RESOLVED** in `MS-VERTEX/main.tex`, `MS-VERTEX/CLAIMS.md`, and `MS-VERTEX/refs.bib`. No mathematical or numerical discrepancy was found across Theorems V1--V4, Proposition `prop:obstruction`, or the lattice section L1; all fast source replays and both Lean 4.35 formalizations (`G1-VERTEXBOUND/Core.lean` and `G1-SIGNED-LEAN/Signed.lean`) pass with exit 0, and a clean Tectonic rebuild in `C:/t/MsVertexDelta` passes with exit 0 and zero overfull boxes. Six small defects remain (one release-convention leak of an internal machine name in the PDF, undefined internal theorem/phase tags in PDF prose, one minor overclaim/caption tension, one wording inaccuracy in the Item 7 addition, and two minor line-anchor slips in `CLAIMS.md`). Fix Findings 1--5 before public release.

---

## 1. Disposition of the six REVIEW2 PARTIAL items and four delta-review items

| REVIEW2 item | Status | Evidence in revised files |
|---|---|---|
| **7** (joint measure: equality vs bound, $K=0$ case) | **RESOLVED** (see minor wording note in Finding 5) | `main.tex:110-118`, `CLAIMS.md:15`: states explicitly that for $K>0$ the measure is a.e. equal to the product of Eq.~\eqref{eq:joint} and $\Theta$ (not merely bounded by it), and separates $K=0$ as an envelope upper bound. |
| **8** (V1 inflation limit in intro) | **RESOLVED** | `main.tex:71`, `CLAIMS.md:16`: states "tends to $1$ as $\Omega_{\max}/W\to0$ at fixed $\lambda$". (Note: `FINISH.md` mislabeled its Item 8 heading as a duplicate of Item 7, but `main.tex:71` has the required text.) |
| **12** ($63/W^3$ fold-cost scope) | **RESOLVED** | `main.tex:348`, `CLAIMS.md:21`: explicitly restricts the fold costs $1,7,56,672$ and $63/W^3$ bound to centered reference lines with fixed nonzero widths at external core momenta, not arbitrary dressed lines or exterior momenta. |
| **14** (`lem:local` hypotheses and folded `4X` convention) | **RESOLVED** | `main.tex:297,309,311`, `CLAIMS.md:17`: `lem:local` explicitly inherits the centered flat-band normal-branch hypotheses of Theorem~\ref{thm:V1} ($\chi_n=0$, positive reflected widths), and line 311 explains that reflection halves both the full-line crossed form and $\|f\|_2^2$, leaving $4X$ unchanged on the half-line. (Note: `FINISH.md` misnumbered this as the verification-description item, but both edits are present in `main.tex`.) |
| **15** (`PROVED` reviewed-proof wording; `CERTIFIED` split; binary64 assumptions) | **RESOLVED** | `main.tex:130,392`, `CLAIMS.md:7`: restores "passed an independent review" in `PROVED`; separates exact rational/interval peak guards from padded IEEE 754 binary64 round-to-nearest signed-tube/radial (including static) enclosures, with stated elementary-operation and BLAS accuracy assumptions. |
| **22** (`CLAIMS.md` disposition and obstruction row) | **RESOLVED** | `CLAIMS.md:18,38`: obstruction row explicitly records bare ($\lambda$) and ME-dressed ($\lambda/(1+\lambda)$) flat-band pointwise limits and bare-only simple-cubic limit; disposition row updated for `REVIEW2.md` and the V4 delta fixes. |
| **Delta V4-1** ($\Omega/W$ relation in coverage-map text) | **RESOLVED** | `main.tex:229`: gives $\Omega/W=x/(\pi\Cinf\eta\sqrt{\lambda})$ and notes that panels are not heating paths at fixed $\Omega/W$. |
| **Delta V4-2** (strictly above $\Tceil$ sampled maximum) | **RESOLVED** | `main.tex:225`, `CLAIMS.md:28`: "Among the temperatures strictly above $\Tceil$ sampled in this run, the largest pairing top was $0.89644635$". |
| **Delta V4-3** (static-limit description in `lem:static` and `thm:V4`(i)) | **RESOLVED** | `main.tex:204,361`, `CLAIMS.md:24`: describes decoupling into per-frequency blocks, normalized spatial convolution of zero-extended lines, and the unsymmetrized vertex $g+g^2(3\,r_n*r_n-Q_n*Q_n)$ between $|z_n(p)|$ and $|z_n(k)|$. |
| **Delta V4-4** (Table 1 and trust-base verification scope) | **RESOLVED** (see Finding 3 on `"phase 3"`) | `main.tex:373,385,392`: distinguishes finite algebra, arithmetic enclosures, formula mutants, and receipt bookkeeping; qualifies `grid_receipts.py` as receipt audit only. |

---

## 2. Audit of `finish.diff`, `finish-claims.diff`, `finish-refs.diff`, and full-paper pass

- **Diff integrity**: `diff -u prefinish/main.tex main.tex`, `diff -u prefinish/CLAIMS.md CLAIMS.md`, and `diff -u prefinish/refs.bib refs.bib` match `finish.diff`, `finish-claims.diff`, and `finish-refs.diff` byte-for-byte.
- **Theorem statements, labels, and scope (V1--V4, L1)**: Checked against `G1-VERTEXBOUND/{RESULT,RIGOR}.md`, `G1-SIGNED/{RESULT,RIGOR}.md`, `G1-SIGNED-LEAN/RESULT.md`, `G1-ONSET/{RESULT,RIGOR}.md`, `G1-GAP/{RESULT,RIGOR}.md`, `G1-LATTICE/RESULT.md`, and all independent reviews (`G1-*-REVIEW/REVIEW*.md`). Every stated constant, hypothesis, rational guard, and scope qualifier matches its reviewed source.
- **Release conventions check**:
  - Author line (`main.tex:35`): `Davis Rens \\ SuperMatics Inc. \\ ORCID 0000-0001-5839-5197` — matches `ADIABATIC/EXTHOL-PAPER/main.tex:16`.
  - Peer-review sentence (`main.tex:446`): `"This preprint has not been peer reviewed."` — present.
  - AI disclosure (`main.tex:446`): `"AI agents generated the proofs, code and drafts under D. Rens' direction. Independent reviews were performed by separate AI agents."` — present.
  - Em dashes: `0` in `main.tex`, `CLAIMS.md`, `refs.bib`, and compiled `main.pdf` (verified via Python/PyMuPDF).
  - Workspace directory/lane names (`G1-VERTEXBOUND`, `G1-SIGNED`, `G1-ONSET`, `G1-GAP`, `G1-LATTICE`, `G3`, `GAMMA`, `D:/dev`): `0` occurrences in `main.pdf`. However, `"Mac"` (`main.tex:402`), undefined internal theorem codes `V1, V2, V3a, V4, B2` (`main.tex:67,107,188,356`), and undefined `"phase 3"` (`main.tex:385,392`) remain in `main.pdf` (Findings 1--3).

---

## 3. Numbered findings

### Finding 1 (P1 — Release convention: internal machine name in PDF text)
- **Location**: `MS-VERTEX/main.tex:402`
- **Issue**: The parenthetical `(Mac audit, exit 0)` appears verbatim in the compiled PDF text, leaking an internal machine name (whereas Table 1 line 383 already replaced `(Mac)` with `(second machine)`).
- **Exact replacement** (`main.tex:402`):
  - Replace: `(Mac audit, exit 0)`
  - With: `(independent discretised-equation audit, exit 0)`

### Finding 2 (P2 — Undefined internal theorem codes `V1`, `V2`, `V3a`, `V4`, `B2` in PDF prose)
- **Location**: `MS-VERTEX/main.tex:67`, `main.tex:107`, `main.tex:188`, `main.tex:356`
- **Issue**: Internal campaign labels `V1`, `V2`, `V3a`, `V4`, and `B2` appear in the manuscript text without ever being defined in `main.tex`. A reader sees `Theorem 3.1`--`Theorem 3.7`, not `V1`--`V4` or `B2`.
- **Exact replacements**:
  - `main.tex:67`: replace `Most bounds use the regular branch; V1 applies to any solution with the specified centering and widths.` with `Most bounds use the regular branch; Theorem~\ref{thm:V1} applies to any solution with the specified centering and widths.`
  - `main.tex:107`: replace `The spectrum is specified separately: V2(iii)--(iv) and V4 are Einstein, while V3a and the uniform B2 statement permit general spectra. V3a also treats a balanced spectator circle.` with `The spectrum is specified separately: Theorem~\ref{thm:V2}(iii)--(iv) and Theorem~\ref{thm:V4} are Einstein, while Theorem~\ref{thm:V3a} and the uniform Theorem~\ref{thm:B2} statement permit general spectra. Theorem~\ref{thm:V3a} also treats a balanced spectator circle.`
  - `main.tex:188`: replace `\begin{theorem}[Overlap of the two regimes (B2)]\label{thm:B2}\hfill\PROVED\\` with `\begin{theorem}[Overlap of the two regimes]\label{thm:B2}\hfill\PROVED\\`
  - `main.tex:356`: replace `The spherical Einstein splice joins V2(iii) to the Young criterion; the general-spectrum splice, including the balanced extension, joins V3a to that criterion.` with `The spherical Einstein splice joins Theorem~\ref{thm:V2}(iii) to the Young criterion; the general-spectrum splice, including the balanced extension, joins Theorem~\ref{thm:V3a} to that criterion.`

### Finding 3 (P2 — Undefined internal campaign tag `"phase 3"` in PDF text)
- **Location**: `MS-VERTEX/main.tex:385`, `main.tex:392`
- **Issue**: `"phase 3"` is never defined in `main.tex`; in the PDF it reads as an unexplained internal project phase. Using the manuscript's own terminology from `Proposition~\ref{prop:V4grid}` ("28-node momentum-dependent") preserves the exact qualification requested in `REVIEW2.md` while remaining self-contained.
- **Exact replacements**:
  - `main.tex:385`: replace `the last script audits phase 3 receipts only` with `the last script audits the 28-node momentum-dependent receipts only`
  - `main.tex:392`: replace `The phase 3 receipt audit pins the exact 28-node set, the worst pairing top $0.72827964$ and the worst preconditioned normal contraction bound $q\le0.70929592$.` with `The 28-node momentum-dependent receipt audit pins the exact 28-node set, the worst pairing top $0.72827964$ and the worst preconditioned normal contraction bound $q\le0.70929592$.`

### Finding 4 (P2 — Scope overclaim in Data Availability and minor caption tension in Table 1)
- **Location**: `MS-VERTEX/main.tex:387`, `main.tex:446`
- **Issue**: (a) `main.tex:446` claims `"Every number in this paper is computed from the model parameters stated in the text, and the listed scripts reproduce it."` Table 1 only lists R1--R4 scripts, not the R0 van Hove comparison or R5 lattice scripts. (b) Table 1's caption (`main.tex:387`) says `"No verification artifact was rerun in preparing the manuscript."`, while line 392 reports a separate spot rerun of `radial_averaged.py` at `(1,1,1/2)`.
- **Exact replacements**:
  - `main.tex:387`: replace `No verification artifact was rerun in preparing the manuscript.` with `No verification artifact was rerun by the manuscript author, other than the single spot check noted in the trust-base paragraph below.`
  - `main.tex:446`: replace `Every number in this paper is computed from the model parameters stated in the text, and the listed scripts reproduce it.` with `Every number in this paper is computed from the model parameters stated in the text, and the scripts in records R0 to R5 reproduce it.`

### Finding 5 (P3 — Wording inaccuracy in Item 7 addition on $\Theta$)
- **Location**: `MS-VERTEX/main.tex:118`
- **Issue**: Line 118 says `"$\Theta$ is the exact density of the angular variable at fixed radii"`. In `G1-ONSET/RIGOR.md:32-35`, $\Theta = |[|R_a-K|, R_a+K] \cap [|R_b-R_c|, R_b+R_c]|/(2K)$ is the normalized admissible-distance overlap factor obtained after integrating out the angular variables at fixed radii, not the density of an angular variable.
- **Exact replacement** (`main.tex:118`):
  - Replace: `$\Theta$ is the exact density of the angular variable at fixed radii, and the fold estimates below use this exact density.`
  - With: `$\Theta$ is the exact normalized angular-overlap factor at fixed radii, and the fold estimates below use this exact joint density.`

### Finding 6 (P3 — Minor line-anchor slip in `CLAIMS.md`)
- **Location**: `MS-VERTEX/CLAIMS.md:14`, `CLAIMS.md:16`
- **Issue**: `CLAIMS.md:16` lists `139–152` for `thm:V1, eq:V1`, overlapping `prop:obstruction` (which starts at line 149); `thm:V1` ends at line 147. In `CLAIMS.md:14`, `def:skeleton` starts at line 94 (line 89 is `eq:moments`).
- **Exact replacements**:
  - `CLAIMS.md:14`: replace `| \`def:skeleton\`, \`eq:skeleton\`, \`lem:signed\` | 89–100, 284–294 |` with `| \`eq:moments\`, \`def:skeleton\`, \`eq:skeleton\`, \`lem:signed\` | 89–100, 284–294 |`
  - `CLAIMS.md:16`: replace `| \`thm:V1\`, \`eq:V1\` | 69–75, 139–152, 325–333 |` with `| \`thm:V1\`, \`eq:V1\` | 69–75, 139–147, 325–333 |`

---

## 4. Journal-grade assessment (e.g., *Phys. Rev. B*)

The underlying mathematical results—especially the sharp local flat-DOS ceiling (Theorem~\ref{thm:V1}), the momentum-resolved signed crossed form and negative-form witnesses (Theorem~\ref{thm:V2}), the spherical weak-coupling/Young-ball windows (Theorems~\ref{thm:V3a}--\ref{thm:B2}), and the computer-assisted static/finite-frequency coverage (Theorem~\ref{thm:V4}, Propositions~\ref{prop:V4grid}--\ref{prop:V4static})—are sound, clearly scoped, and of publishable technical caliber for a *Phys. Rev. B* or *J. Math. Phys.* methods/bounds paper, **but the manuscript in its current form is not yet ready for *Phys. Rev. B* submission**. Reaching PRB grade would require: (i) making the proofs self-contained in appendices or a cleanly structured Supplemental Material rather than deferring core steps to 15 inline citations of workspace `RIGOR.md` files and `R0`--`R5` record codes; (ii) recasting into standard REVTeX4-2 format with fuller physical context tying the signed-vertex cancellations and 1D lattice numerics (and their CDW/RPA caveats) to the nonadiabatic Eliashberg literature; (iii) fixing Findings 1--5 above so no internal campaign tags remain; and (iv) packaging the Lean and Python verification artifacts into a single versioned public archive (assigning the currently placeholder supplement DOI) with relative-path receipt verification so a referee can reproduce Table 1 without hardcoded local directories such as `D:/t/G1Gap`.

---

## 5. Verification run and exit codes

All checks below were executed read-only against `D:/dev/scratch/ceiling-20261118/` (Python scripts via `python -B` or offline `uv run --no-project --with numpy --with mpmath --with scipy --with sympy python -B`; Lean via `leanprover--lean4---v4.35.0/bin/lean.exe -j1` with shared Mathlib overlay).

| Check | Command / target | Exit code | Wall time | Key observed output |
|---|---|---:|---:|---|
| Tectonic rebuild (copy) | `tectonic.exe --keep-logs main.tex` in `C:/t/MsVertexDelta` | **0** | 4.8s | 15 pages, 45,304 chars, 0 overfull boxes, 0 undefined refs/cites, 1 underfull `\hbox` in `main.bbl:8`; 0 em dashes; 0 margin overflows |
| Diff exactness | `diff -u prefinish/<f> <f>` vs `finish*.diff` (3 files) | **0** | <1s | `finish.diff`, `finish-claims.diff`, `finish-refs.diff` all exact |
| R1 exact algebra + mutants | `G1-VERTEXBOUND/check_exact.py` | **0** | 0s | Signed form, row moments, rational guards, 3 mutants; ratio `13427/12985` |
| R1 Lean 4.35 | `G1-VERTEXBOUND/Core.lean` | **0** | 77s | 8 declarations, axioms `[propext, Classical.choice, Quot.sound]` only |
| R2 signed form + mutants | `G1-SIGNED/replay_form.py` | **0** | 7s | Exact Nambu interpolation = signed form; 6 mutants rejected |
| R2 rational guards + mutants | `G1-SIGNED/replay_bounds.py` | **0** | 0s | Margin `55/262144 > 1/8192`; 5 mutants rejected |
| R2 independent checks | `G1-SIGNED-REVIEW/independent_checks.py` | **0** | 3s | Independent polynomial, static limit, and rational guards agree |
| R2 Lean 4.35 | `G1-SIGNED-LEAN/Signed.lean` | **0** | 12s | 13 declarations, axioms `[propext, Classical.choice, Quot.sound]` only |
| R3 geometry/folds/guards + mutants | `G1-ONSET/replay.py` | **0** | 0s | 625 fixtures, 6 mutants rejected |
| R3 independent checks | `G1-ONSET-REVIEW/independent_checks.py` | **0** | 0s | Young `r=1/9`, fold totals `1331,1573,1936,2652`, inflation `52000*B^4` |
| R4 continuum + peak guards + mutants | `G1-GAP/replay.py` | **0** | 10s | `certified=54/135`, `x>=19/2` upper `<=0.70468692`, 5 mutants rejected |
| R4 radial fixtures + mutants | `G1-GAP/radial_check.py` | **0** | 2s | `RADIAL FIXTURES PASS: 3/3 mutants rejected` |
| R4 finite rational rounding check | `G1-GAP/finite_check.py` | **0** | 1s | `FINITE CHECK inverse-moment rational rounding fixture PASS` |
| R4 finite published grid + mutants | `G1-GAP/finite.py --published` | **0** | 17s | 3 mutants rejected; `NEW closed=15/59`, `worst pairing upper<=0.50417617` |
| R4 under-resolved static control | `G1-GAP/static.py --lam 2 --xs 1/16 --steps 1` | **1** (expected) | 6s | Rejects under-resolved slab with `ArithmeticError` |
| R4 independent checks | `G1-GAP-REVIEW/independent_checks.py` | **0** | 2s | 7 peak anchors and near-static remainders pass |
| R4 receipt audit | `G1-GAP/grid_receipts.py logs/G1Gap` | **0** | 0s | `eta=1 closed=32/45 (radial=12)`, `eta=3/2 closed=40/45 (radial=16)`, `pairing<=0.72827964`, `q<=0.70929592` |

### Checks not completed / skipped
- `G1-GAP/static.py --lam 1/9 1/4 1/2 1 2 --xs 1/16 ... 16 --steps 80` (full 45-node static grid) and `--independent`: timed out at 300s (reached $\lambda=1, x=1/8$; a separate 120s spot run at the worst node $\lambda=2, x=1/16$ confirmed `w=1` upper bound `0.70431954 <= 0.70431955` before timing out on higher $w$ blocks).
- `G1-GAP-REVIEW/phase3_checks.py`: exited 1 (`FileNotFoundError`) because it hardcodes `D:/t/G1Gap`, whereas the archived logs live in `D:/dev/scratch/ceiling-20261118/logs/G1Gap` (verified directly via `grid_receipts.py` above).
- Multi-hour 28-node `radial_averaged.py` / `radial_moments.py` runs and R5 `G1-LATTICE` numerical sweeps: not rerun (exceed the <10 min fast-replay budget).
- Temporary copy `C:/t/MsVertexDelta` was deleted after verification.
