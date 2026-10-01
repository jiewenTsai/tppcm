## R CMD check results

Local `R CMD check --as-cran --run-donttest` (macOS 27, R 4.5, TAM 4.3-25,
sirt 4.2-133, 2026-10-02): 0 errors | 0 warnings | 3 notes

* New submission.
* "unable to verify current time": the check machine had no network access
  (not a package issue).
* HTML validation skipped: the local HTML Tidy is too old (not a package issue).

`urlchecker::url_check()` and `spelling::spell_check_package()`: see the
session log; no broken URLs.

Cross-platform checks (win-builder R-devel, macOS builder, R-hub) have not
been run yet; they require the maintainer's mailbox and a GitHub repository.

## Notes for submission

* Two vignettes in English and two in Traditional Chinese (UTF-8). The
  Chinese ones duplicate the English ones with identical code.
* Examples that fit models are wrapped in \donttest{} because they call TAM;
  with --run-donttest every example file runs in under 1 second.
* Suggested packages (mirt, psychotools, psychotree, partykit, strucchange,
  numDeriv) are used only in examples, vignettes and tests, with
  requireNamespace() guards in the functions.
