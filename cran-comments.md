## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

Local:

* macOS 27, R 4.5, TAM 4.3-25, sirt 4.2-133:
  `R CMD check --as-cran --run-donttest`: 0 errors | 0 warnings | 3 notes
  (new submission; "unable to verify current time" because the machine
  had no network access; HTML validation skipped because the local HTML
  Tidy is too old. The last two are local and not package issues.)

GitHub Actions (r-lib/actions check-standard, `--as-cran`, 2026-10-01):

| Platform | R | Status |
|---|---|---|
| macOS (aarch64) | 4.6.1 | OK |
| Windows Server 2022 (x86_64, ucrt) | 4.6.1 | OK |
| Ubuntu 24.04 (x86_64) | 4.6.1 | OK |
| Ubuntu 24.04 (x86_64) | 4.5.3 (oldrel-1) | OK |
| Ubuntu 24.04 (x86_64) | R-devel (2026-09-30 r90605) | OK |

All 230 testthat expectations pass on every platform.

R-hub v2 (2026-10-01; all three jobs completed successfully):

| Platform | R |
|---|---|
| linux (Ubuntu 24.04, x86_64) | R-devel |
| windows | R-devel |
| macos (x86_64) | R-devel |

win-builder (R-devel): PENDING (to be run by the maintainer with
`devtools::check_win_devel()`)

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
