## Resubmission

This is a resubmission of tppcm 0.1.0. It replaces the submission of
2 October 2026, which is still waiting for manual review (its auto-check
passed with 1 NOTE: new submission and possibly misspelled words). The
version number is unchanged because the package has never been on CRAN.

Changes since that submission:

* `inst/CITATION` uses `meta$Version` instead of
  `utils::packageVersion("tppcm")`, so it can be read when the package is
  not installed (requested by Uwe Ligges on 2 October 2026).
* New interface: submodels are specified by a design formula for the step
  discriminations (`design = ~ item * step` etc.) in `tppcm()`,
  `xxirt_tppcm()`, `tppcmtree()` and `score_test()`.
* New functions for differential step functioning: `dsf_design()` and
  `dsf()`, and a per-step test in `dif_test(by = "step")`.
* `inst/WORDLIST` lists the author names and acronyms in DESCRIPTION that
  the spell check flagged.

See NEWS.md for the full list.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

Local:

* macOS, R 4.5: `R CMD check --as-cran`: 0 errors | 0 warnings.
  The PDF manual was not built locally (inconsolata.sty is missing); it was
  built on GitHub Actions and win-builder (below).

GitHub Actions (r-lib/actions, `R CMD check --as-cran` including the PDF
manual, DATE_GHA):

GHA_TABLE

win-builder: WINBUILDER_RESULT

`urlchecker::url_check()` and `spelling::spell_check_package()`: no broken
URLs, no misspellings.

## Notes for submission

* Two vignettes in English and two in Traditional Chinese (UTF-8). The
  Chinese ones duplicate the English ones with identical code.
* Examples that fit models are wrapped in \donttest{} because they call TAM;
  with --run-donttest every example file runs in under 1 second.
* Suggested packages (mirt, psychotools, psychotree, partykit, strucchange,
  numDeriv) are used only in examples, vignettes and tests, with
  requireNamespace() guards in the functions.
