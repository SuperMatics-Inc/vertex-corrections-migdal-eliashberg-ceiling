# BUILD.md, version 3.0.0

Build record for the released version 3 of *Vertex corrections and the Migdal-Eliashberg ceiling: a geometric dichotomy between nested and curved Fermi surfaces*, tagged `v3.0.0` at the deposited tree. Every hash below is SHA-256 of the file bytes as released. `MANIFEST.sha256` pins the whole tree, `ADMITTED.sha256` the same bytes in the form accepted for release, and CI verifies both.

## 1. Manuscript and paper record

| File | Bytes | sha256 |
|---|---|---|
| `paper/main.tex` | 118647 | `9a4944297579fbdfdaa5590e59a0f7d06b504f7c729de60057fd7405c4401e41` |
| `paper/refs.bib` | 7931 | `a1cbb440feaa4dadd979061789e552dc8bc33762dabc0b1ceea494a910077e07` |
| `paper/main.pdf` | 299526 | `ee9f15012993dd2b445c85cc8e0f0cc3868d7ee240fbc013a752509af3ea60d6` |

Pages 31. The paper record `10.5281/zenodo.24415577` carries `main.pdf`; the sources travel in this repository.

The deposited manuscript is the reviewed text, byte for byte: no editorial byte was changed after review, and `refs.bib` is unchanged from version 2. The title block and the data-availability paragraph of the manuscript still name the version 2 records, `10.5281/zenodo.24410274` and `10.5281/zenodo.24410273`, because the PDF was not rebuilt for the release and the released bytes are the reviewed bytes. The current records of this version are `10.5281/zenodo.24415577` for the paper and `10.5281/zenodo.24415572` for the supplement; that is a metadata-only difference between the manuscript's DOI literals and the record identifiers, and it is recorded as such in the release ledger. The repository, the release notes, `CITATION.cff` and `README.md` carry the version 3 DOIs.

## 2. How the PDF was built

Toolchain: Tectonic 0.17.0, XeTeX engine, from the TeX Live bundle the tool downloads.

```sh
rm -f main.aux main.bbl main.log main.blg main.out
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

The build reports 31 pages, zero overfull boxes, zero undefined references and zero unresolved citations; the single recorded box warning is one underfull `\hbox` (badness 2662) on the closing line of the abstract, which does not affect the text. A second run of the identical command produced a byte-identical PDF, so the hash above is reproducible on this Tectonic build. `SOURCE_DATE_EPOCH` fixes the PDF timestamps.

The PDF is not byte-reproducible across Tectonic versions or machines, so CI does not compare its hash. CI rebuilds the manuscript from `paper/main.tex`, fails on any overfull box, undefined reference or unresolved citation, and checks the text of the rebuilt PDF: the DOI links it carries, the repository URL, the byline, the peer-review sentence, the AI disclosure, and the absence of the banned tokens listed there.

## 3. New record R7

24 files, all byte-identical to the reviewed certificate record, which is the source of every number of Theorem 4.8 of the manuscript.

| File | Bytes | sha256 |
|---|---|---|
| `records/R7-flat-sheet-certificate/CERT.md` | 33327 | `5e761b02d9cd847332d2239e1f43a10232386ada984d630b8ce2960595a9fd6d` |
| `records/R7-flat-sheet-certificate/RESULT.md` | 6235 | `9084d0aae887551d5454d213168d042a81adc9b439280059379027ec53c69d18` |
| `records/R7-flat-sheet-certificate/cert.py` | 39945 | `7014e3c3190d0f43824c2abfeb618dc38a07c87efdbb5da0a08d218076e1230b` |
| `records/R7-flat-sheet-certificate/independent.py` | 52098 | `098cf0f7655b3c5286bd0dfe0cbac0c9b377abde20ec39db32e578681b897880` |
| `records/R7-flat-sheet-certificate/cert_default.log` | 2540 | `e6cf2a2c069c2fca381294a2d4286940401db5f5ab1aecb6960e83cf007ca3cb` |
| `records/R7-flat-sheet-certificate/cert_main.log` | 1189 | `bb22b80a9dd0d27132fde64d69d7dc62ed44ad8a3da47eb497ef8919919b55f7` |
| `records/R7-flat-sheet-certificate/cert_m200.log` | 1187 | `7ba513d544ba68a0b72b69dcc01d9dd1559c4c0a751f26ab9ed891c2fae0fda0` |
| `records/R7-flat-sheet-certificate/cert_prec400.log` | 1189 | `a735ea0373798a14ac413d1823f169764b249b279652fc40459f492524197af7` |
| `records/R7-flat-sheet-certificate/cert_prec400_m200.log` | 1187 | `ef5947de2c093754560c91b0ca412858a6df4f34ed14d0360b0a54c17e39231b` |
| `records/R7-flat-sheet-certificate/independent_nk8.log` | 1358 | `a3ff22c21e1d6c2d155c9b95f4eac0cf6c78f2d3d63e70e967e4c73ba1e840cf` |
| `records/R7-flat-sheet-certificate/independent_nk24.log` | 1364 | `b9c66f3825d376f99fd4a0205a1522ee0cbde8f1b31bf28a9c2328b00f57cb29` |
| `records/R7-flat-sheet-certificate/independent_nk24_m200.log` | 1045 | `8ac2340be927d11c569906115af093c9c673d4c10e5299c5af8cfa96a1276d1b` |
| `records/R7-flat-sheet-certificate/independent_nk24_m200_rerun.log` | 1369 | `69e8a8dc372b7c0158d3614db59d220ee682968885680daef7bcf1dad23a9e64` |
| `records/R7-flat-sheet-certificate/ladder.log` | 4322 | `61f17628caae1745d491c041cd0e723a7d053f7a443cb24e8815b2edce7040d1` |
| `records/R7-flat-sheet-certificate/mutants/diag_drop.log` | 1191 | `27e4d23b51899d23d0ff3e24f4745aa2c9a0807eb8eefb94eb19aac2631cb99a` |
| `records/R7-flat-sheet-certificate/mutants/drop_lb_t3.log` | 1191 | `1ef255b27fbafd9b364675edb86170c18bff63257b1f08f4ee9dd28e42c30b89` |
| `records/R7-flat-sheet-certificate/mutants/float_leak.log` | 1233 | `1739797daa60c8e24adf4442746c71ac9427ff9aa6474d55d2d3fba7effccd5f` |
| `records/R7-flat-sheet-certificate/mutants/m_threshold.log` | 1194 | `f416bca1aafaabeeaf36179637b433366b87df687cec50624412868305e34b92` |
| `records/R7-flat-sheet-certificate/mutants/no_vertex.log` | 1281 | `29714a395c9daf874a55e2452c09cec424835eaca694e5674b4ae335ef22a7d9` |
| `records/R7-flat-sheet-certificate/mutants/point_not_ball.log` | 1228 | `851d7a9118e3e936d24dddfff12d609c6402a1df8a8377b1e08401c27ebd9f6b` |
| `records/R7-flat-sheet-certificate/mutants/ray_thresh.log` | 1192 | `17c53ae1cfa3a6c5a5d6067791c3c0efbfa43cba40cf17033abfdcf4023c94ea` |
| `records/R7-flat-sheet-certificate/mutants/rb_scale.log` | 1174 | `0ed7c1b889696b633dc476a50951683f52f07f4741489edacc38e227ce2f9222` |
| `records/R7-flat-sheet-certificate/mutants/res_scale.log` | 1189 | `bcbb604bd7fbcdcdc7cf450b2e0c9d539b91500e26a95747b25adcbfd3286329` |
| `records/R7-flat-sheet-certificate/mutants/rho_zero.log` | 1188 | `55b0421899cf7ae0087d75cbc009ed520dfa2ce73536ec4f4ebaa6a03a6cf102` |

Two properties of the record to note before reading it. `CERT.md` is the working record of the computation and its runtime section names the memory guard that killed one launched run and the local temporary directory of that machine; the files are shipped unmodified because the reviewed numbers and their line references are pinned to these bytes. The record's review documents, superseded drafts and byte code are not part of it and are not in this tree.

The certificate of the record runs in CI as its own job, on Python 3.13 with `mpmath==1.4.1`, `numpy==2.5.5`, `scipy==1.18.1`, `python-flint==0.9.0` and `gmpy2==2.3.2`. That job runs `python -B cert.py` at its default settings (both parameter sets, one process, about 77 s on the reference machine, 81 s when replayed for this release) and checks the printed certified bounds, then runs both self-tests of the second checker. The ten mutant runs (38 s to 52 s each on the reference machine), the 400-bit replays and the long second-checker runs at `N_k=8` and `N_k=24` are not rerun in CI; their captured logs are shipped above and their hashes are pinned both by the manifests and, for `cert_default.log` and `independent_nk24_m200_rerun.log`, by literal lines of the workflow.

## 4. Diff against version 2

`paper/main.tex` differs from the released version 2 text in 21 hunks, 123 added and 17 removed lines:

1. The leading comment records the labels R0 to R7.
2. The abstract adds the certification sentence and extends the "on the negative side" and "earlier results remain" sentences so that the certified finite-model statement is separated from the numerical and open continuum statements. No version 2 sentence is deleted.
3. The introduction summary gains one bullet for the certified finite-model crossing, and the "What is not claimed" paragraph gains the finite-model caveat.
4. A new paragraph, `Version 3 note`, states what changed against version 2.
5. A new Section 4.4, `A certified crossing at two discretisation points`, holds Definition 4.6 (the finite model), Proposition 4.7 (real symmetry of the gap operator, with proof) and Theorem 4.8 (parts (N), (K), (S) and (V)), the method paragraph, the precision and mutant paragraph with its table, the adjacent-discretisation paragraph and the scope paragraph. The preceding "three gaps" paragraph now records Gap 1 as closed at these discretisation points.
6. Open problem `prob:flatcert` states the same, with Gaps 2 and 3 untouched, and warns that the inverse at the truncations of Proposition `prop:flatcert` acts in a space of order `10^5`.
7. The supplementary record table gains the row `R7 finite model certificate`, the reproduction section gains the subsection `Theorem 4.8` and the `Reproduction of R7` paragraph, and the Lean and trust-base paragraphs are extended.
8. The limitations sentence of the conclusions, the "what remains open" item on end-to-end certification and the closing support paragraph are updated to R0 to R7.

Nothing was withdrawn: every statement, label and number of version 2 is present with the same label. `refs.bib` is unchanged, so the bibliography is that of version 2.

## 5. Supplement archive

`vertex-corrections-migdal-eliashberg-ceiling-replay.zip` in the supplementary record `10.5281/zenodo.24415572` is `git archive` of tag `v3.0.0`, so it holds the records R0 to R7 without review records and without superseded drafts. Its byte count and SHA-256 are recorded in the release receipt rather than here, because a file cannot contain the hash of the archive of the tree that contains it. CI verifies both manifests and runs every replay and the certificate job on exactly this tree before the archive is taken.

## 6. Documents of the release

| File | Bytes | sha256 |
|---|---|---|
| `CLAIMS.md` | 11718 | `3bcc295190f12fd88de595c203c13a4e849a491cc7e1a8f95c420d10e2872729` |
| `README.md` | 14072 | `8e568af76f8da5e4d48f58c977b3249b30552a1c1dce21050058f4ea57128ca7` |

`MANIFEST.sha256` pins every tracked file except the two manifests, and `ADMITTED.sha256` carries the same 155 entries in binary-mode form. `BUILD.md` is not listed here, because it cannot contain its own hash; the manifest pins it.

`SUPPLEMENT-ADDITIONS.md`, which lists the files of record R7 for the supplementary upload, is a custody note of the release and is not part of the tree; its hashes are in the release receipt.
