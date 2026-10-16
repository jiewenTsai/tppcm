## Resubmission

This is a resubmission of tppcm 0.1.0. The submission of 16 October 2026
answered the review by Konstanze Lauseker on 15 October 2026 but was
rejected by the incoming pretest because the overall check time on Windows
was 14 min (> 10 min). The version number is unchanged because the package
has never been on CRAN.

Changes since the submission of 16 October 2026 (check time):

* The two Traditional Chinese vignettes, which repeated the code of the two
  English vignettes, are no longer in the package; they stay in the GitHub
  repository. This halves the time for re-building the vignettes.
* The slowest tests (the multigroup refits in `test-dsf_design.R` and two
  blocks in `test-design.R` that fit many designs) are skipped on CRAN with
  `skip_on_cran()`; they run on GitHub Actions.

Changes in response to the review of 15 October 2026 (unchanged since the
submission of 16 October 2026):

* DESCRIPTION no longer starts with "Tools for".
* All acronyms in the Description are spelled out (generalized partial
  credit model, partial credit model, root mean square error of
  approximation, standardized root mean squared residual, comparative fit
  index, Tucker-Lewis index); TPPCM is defined at its first use.
* Yu (1991) and Glas (1998) are linked in the form authors (year) <https:...>
  (the dissertation has no DOI; the link is its handle in the University of
  Illinois repository, and Glas (1998) links to the Statistica Sinica page).
* The vignettes save the options they change in the setup chunk
  (`old_options <- options(width = 90)`) and restore them in a final hidden
  chunk (`options(old_options)`). The `par()` calls already restored the old
  settings (`op <- par(...); ...; par(op)`). Examples do not change
  `options()`, `par()` or the working directory.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

Local:

* macOS, R 4.5: `R CMD check --as-cran`: 0 errors | 0 warnings.
  The PDF manual was not built locally (inconsolata.sty is missing); it was
  built on GitHub Actions and win-builder (below).

GitHub Actions (r-lib/actions, `R CMD check --as-cran` including the PDF
manual, 2026-10-16):

| Platform | R | Status |
|---|---|---|
| macOS (aarch64) | 4.6.1 | OK |
| Windows Server (x86_64, ucrt) | 4.6.1 | OK |
| Ubuntu 24.04 (x86_64) | 4.6.1 | OK |
| Ubuntu 24.04 (x86_64) | 4.5 (oldrel-1) | OK |
| Ubuntu 24.04 (x86_64) | R-devel | OK |

"OK" means 0 errors and 0 warnings; the PDF manual was built and checked on
every platform. The one NOTE on the runners comes from the HTML manual
check (HTML Tidy or V8 not available on the runner, so that validation was
skipped). All 374 testthat expectations pass on every platform.

Check time on the GitHub Windows runner: re-building the vignettes 187 s
(373 s before the Chinese vignettes were removed); the test suite 21 s with
`NOT_CRAN=false` (99 s with all tests).

win-builder (R-devel 2026-10-05 r90641 ucrt, 2026-10-06; run on the same
sources with the development version number 0.1.0.9000): 0 errors |
0 warnings | 1 note, including the PDF manual.

## Notes for submission

* Examples that fit models are wrapped in \donttest{} because they call TAM;
  with --run-donttest all examples together take about 35 s on the GitHub
  Windows runner.
* Suggested packages (mirt, psychotools, psychotree, partykit, strucchange,
  numDeriv) are used only in examples, vignettes and tests, with
  requireNamespace() guards in the functions.
