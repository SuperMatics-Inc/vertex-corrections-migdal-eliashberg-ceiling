# BUILD.md, version 2.0.0

Build record for the released version 2.0.0 of *Vertex corrections and the Migdal-Eliashberg ceiling: a geometric dichotomy between nested and curved Fermi surfaces*. Every hash below is SHA-256 of the file bytes as released. `MANIFEST.sha256` pins the whole tree, `ADMITTED.sha256` the same bytes in the form accepted for release, and CI verifies both.

## 1. Manuscript and paper record

| File | Bytes | sha256 |
|---|---|---|
| `paper/main.tex` | 96929 | `ad5ebc1fedbd44683055683379107011f3d8a1bf05147a112d0efbcc967243aa` |
| `paper/refs.bib` | 7931 | `a1cbb440feaa4dadd979061789e552dc8bc33762dabc0b1ceea494a910077e07` |
| `paper/main.pdf` | 263168 | `0de7aad87ac5bc05f7c01a7c3b7eea3f0c150f0b900a1c0c1f1b148226a690f0` |

Pages 26. The paper record `10.5281/zenodo.24410274` carries `main.pdf`; the sources travel in this repository.

The reviewed manuscript is the file whose `\MSPAPERDOI` and `\MSSUPPDOI` macros still hold the first version's DOIs, sha256 `ae60f330b9b46f1cd66fc57488754207b35a6b7320111d9b071e8b91ddd85ae8`. `paper/main.tex` differs from it on exactly two lines: the two DOI macros now hold `10.5281/zenodo.24410274` and `10.5281/zenodo.24410273`. Nothing else in the manuscript changed, and the release check re-derives the reviewed file from the released one by restoring those two literals and comparing the hashes.

## 2. How the PDF was built

Toolchain: Tectonic 0.17.0, XeTeX engine, from the TeX Live bundle the tool downloads.

```sh
rm -f main.aux main.bbl main.log main.blg main.out
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

The build was run on `paper/main.tex` in a directory holding `main.tex`, `refs.bib` and nothing else from the release. It reports 26 pages, zero overfull boxes, zero undefined references and zero unresolved citations; the single recorded box warning is an underfull `\hbox` in the Chubukov reference entry (badness 2662), which does not affect the text. `SOURCE_DATE_EPOCH` fixes the PDF timestamps.

The PDF is not byte-reproducible across Tectonic versions or machines, so CI does not compare its hash. CI rebuilds the manuscript from `paper/main.tex`, fails on any overfull box, undefined reference or unresolved citation, and checks the text of the rebuilt PDF: the two DOIs, the repository URL, the byline, the peer-review sentence and the AI disclosure, and the absence of the banned tokens listed there.

## 3. New record R6

| File | Bytes | sha256 |
|---|---|---|
| `records/R6-flat-sheet/RIGOR.md` | 37855 | `dc935bc9473801595a427f36a26cae6d694b209e0eeaed04a7365e8d2c005392` |
| `records/R6-flat-sheet/RESULT.md` | 10982 | `c15e3380c1e9987f04b802108dc2417b8586acea06c47a6584423a357588ecf3` |
| `records/R6-flat-sheet/replay.py` | 10496 | `08a16090313453b516a08514057803fd16eaf8cac642f9e7f30863167e6fa776` |
| `records/R6-flat-sheet/verify_1d_counterexample.py` | 20476 | `5711ed271c2275fbd16a0416aa15869145cff3415f0acad665e8d9ca29c3300c` |
| `records/R6-flat-sheet/verify_1d.log` | 5252 | `a82fcf97f95f6f6f7707a4890e786cd6c8ea39f0b07fbe1036fc2df73e4e0a7d` |

`replay.py`, the solver and the captured log are byte-identical to the reviewed computation. `RIGOR.md` and `RESULT.md` were reissued for the upload with the working record names, filesystem paths and working-file list replaced by the record labels `R2` and `R3` and by repository-relative commands; their line counts are unchanged at 163 and 61, so every line reference in `CLAIMS.md` and in the manuscript remains valid.

## 4. Supplement archive

`vertex-corrections-migdal-eliashberg-ceiling-replay.zip` in the supplementary record `10.5281/zenodo.24410273` is `git archive` of tag `v2.0.0`, so it holds the tree of records R0 to R6 without review records and without superseded drafts. Its hash is recorded in the release receipt of the Zenodo upload, not here, because a file cannot contain the hash of the archive of the tree that contains it.

## 5. Documents of the release

| File | Bytes | sha256 |
|---|---|---|
| `CLAIMS.md` | 12452 | `97dccf98085a1ed15a4a91d9c5724ced3b740d706ab0a1621084d28dc05db847` |
| `README.md` | 10618 | `92cc993bdc114d22f9163124c3bfa34c4403c51c6536b6a253daad0c25566d2d` |

`SUPPLEMENT-ADDITIONS.md`, which lists the files of record R6 for the supplementary upload, is a custody note of the release; its hashes are in the release receipt.
