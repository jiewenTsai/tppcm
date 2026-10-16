## Resubmission

This is a resubmission of tppcm 0.1.0, in response to the review by
Konstanze Lauseker on 15 October 2026. The version number is unchanged
because the package has never been on CRAN.

Changes:

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
manual, 2026-10-06):

| Platform | R | Status |
|---|---|---|
| macOS (aarch64) | 4.6.1 | OK |
| Windows Server (x86_64, ucrt) | 4.6.1 | OK |
| Ubuntu 24.04 (x86_64) | 4.6.1 | OK |
| Ubuntu 24.04 (x86_64) | 4.5 (oldrel-1) | OK |
| Ubuntu 24.04 (x86_64) | R-devel (2026-10-03 r90638) | OK |

"OK" means 0 errors and 0 warnings; the PDF manual was built and checked on
every platform. The one NOTE on the runners comes from the HTML manual
check (HTML Tidy or V8 not available on the runner, so that validation was
skipped). All 374
testthat expectations pass on every platform.

win-builder (R-devel 2026-10-05 r90641 ucrt, 2026-10-06; run on the same
sources with the development version number 0.1.0.9000): 0 errors |
0 warnings | 1 note, including the PDF manual.

## Notes for submission

* Two vignettes in English and two in Traditional Chinese (UTF-8). The
  Chinese ones duplicate the English ones with identical code.
* Examples that fit models are wrapped in \donttest{} because they call TAM;
  with --run-donttest every example file runs in under 1 second.
* Suggested packages (mirt, psychotools, psychotree, partykit, strucchange,
  numDeriv) are used only in examples, vignettes and tests, with
  requireNamespace() guards in the functions.
