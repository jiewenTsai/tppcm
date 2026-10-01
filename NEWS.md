# tppcm 0.1.0

First release. Changes after the internal code review (statistician and
applied-user reviews, October 2026):

* Multigroup fits: the free latent means and variances of groups 2, ..., G are
  part of the parameter vector (scores, Louis information, standard errors);
  they were treated as fixed before, which understated item standard errors.
  `irt_pars()` reports their standard errors.
* `dif_test()` prints a reminder about impact (ability differences between the
  groups of the covariate) for single-group fits; `tppcmtree()` gained
  `impact = ` (group means and variances estimated in every node) and keeps
  rows with missing responses.
* `score_test()` checks that `against` contains the fitted model, and gives a
  clear message when the product form is requested with unequal numbers of
  categories.
* `irt_pars()` blanks the difficulty of steps whose discrimination is near 0
  and lists them; `as.data.frame()` method; `level` is validated.
* Model labels name the design (`GPCM design (tam.mml.3pl)` etc.).
* `xxirt_tppcm()` objects print a short summary; `coef()` for `get_parts()`
  objects; `dif_test()` keeps the item order; `tppcmtree()` prints the splits
  only; input checks and messages in several functions.
* `m2()` documents the MCAR assumption and the difference to mirt's C2 under
  misfit; `mi()` documents the EPC as a one-step approximation.

* `tppcm()` (alias `step_design()`): step design for `TAM::tam.mml.3pl()`;
  saturated TPPCM, GPCM, step model, PCM and any linear restriction through
  `index`; items may have different numbers of categories.
* `xxirt_tppcm()`, `xxirt_vcov()`: product form `a_il = alpha_i gamma_l` and
  the bound `a_il >= 0` with `sirt::xxirt()`.
* `irt_pars()`: compact report with fit header, correct standard errors
  (observed information), long format and parameters of the fitted model.
* `get_parts()`: components of a fit (`what = ` as in `lme4::getME()`),
  with `vcov()`, `logLik()` and `nobs()` methods.
* `step_info()`: Fisher information split into step contributions.
* Tests: `wald_test()`, `score_test()`, `item_test()`, `mi()` (modification
  indices and EPCs), `m2()` (M2, C2, M2*, RMSEA2, SRMSR, CFI, TLI; missing
  data under MCAR).
* DIF: `dif_test()` (score-based tests per step, item or block) and
  `tppcmtree()` (model-based recursive partitioning).
* `sim_tppcm()` for simulating data.
* Vignettes: tutorial and recipes, in English and Traditional Chinese, using
  data sets from TAM, CDM, mirt and psychotools.
