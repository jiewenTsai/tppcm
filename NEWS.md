# tppcm (development version)

## Breaking changes

* Submodels are specified by a design formula for the step discriminations
  (CBD), read as an item x step layout on the log scale,
  `log a_il = mu + alpha_i + gamma_l + (alpha gamma)_il`: `~ 1` (PCM),
  `~ item` (GPCM), `~ step` (discriminations common to all items),
  `~ item + step` (Rank-1, `a_il = alpha_i * gamma_l`) and `~ item * step`
  (saturated TPPCM, default; `~ item:step` is the same model).
  `tppcm(dat, design = )`, `xxirt_tppcm(dat, design = )`,
  `dsf_design(design = )`, `tppcmtree(design = )` and
  `score_test(against = )` take formulas. The `model =` arguments, the model
  names (`"pcm"`, `"gpcm"`, `"step"`, `"rank1"`, `"tppcm"`) and
  `step_design()` were removed without deprecation. `index =` is unchanged.
  `tppcm()` refuses `~ item + step` (multiplicative, not a linear TAM design)
  and points to `xxirt_tppcm()`.
* Tables of tests share one convention: `dif_test()`, `dsf()`, `item_test()`,
  `wald_test()` and `mi()` take `adjust = "holm"` (any `p.adjust()` method)
  and report adjusted p-values in `p_adj`. `dif_test()` returns `p_adj`
  instead of `p_holm` and a new `df` column. Numerical results are unchanged.
* `dif_test()`: the default is now `parm = c("disc", "diff")` (both blocks;
  `"all"` is still accepted). `tppcmtree()` uses the same `parm` values.

## New features: differential step functioning (DSF)

* `dif_test(fit, g, by = "step")` on a multigroup fit
  (`tam.mml.3pl(..., group = g, est.variance = FALSE)`) tests every step
  jointly in its discrimination and location (2 df per group contrast), an LM
  test in the sense of Glas (1998). It is computed in the slope-intercept
  form, so it does not depend on the parameterization and steps with
  discrimination near 0 keep both degrees of freedom. This is the recommended
  DSF screen.
* `dsf_design()` builds multigroup TPPCMs for `TAM::tam.mml.3pl()` with
  lavaan-style `group.equal` (`NULL` configural, `"cbd"` metric,
  `c("cbd", "ibd")` scalar; `"cbl"` for `"ibd"`) and `group.partial` (single
  steps released), including the identification of the latent means and
  variances (`variance.fixed`, `beta.fixed`).
* `dsf(fit, group.equal = c("cbd", "ibd"), group.partial = NULL)`: likelihood
  ratio tests that release every restricted step parameter in turn, with the
  latent means and variances of the groups estimated. CBD tests keep the
  location of the step free, so they do not depend on the origin of the
  scale. `items =` also accepts single steps (`"item:step"`), so the units
  flagged by `dif_test(by = "step")` can be passed directly.

## Improvements

* `dif_test()` on multigroup fits with a factor covariate: the LM test is the
  score test of the model with the tested parameters split by the levels of
  the covariate. The previous Brownian-bridge LM ignored that the latent
  means and variances of groups 2, ..., G are group specific and was
  conservative (about 2% instead of 5%, not improving with N). The default
  `info = "robust"` is the generalized score test of Boos (1992). In
  simulations of per-step tests (two groups, 300 or 600 persons per group,
  with and without local dependence) it was the only variant that kept the
  nominal level in every condition and, at equal size, was as powerful as the
  likelihood ratio test; `info = "louis"` agrees closely with the LR test;
  `info = "opg"` is liberal.
* Model labels and printed headers show the design formula (e.g.
  `GPCM, design ~ item (tam.mml.3pl)`); `get_parts(fit, "formula")` returns
  it. Every table of tests starts with a one-line header stating what was
  tested, against what, and how `p_adj` was adjusted.
* `irt_pars()$group` and `get_parts(, "groups")` use TAM's group labels
  instead of the internal codes 1, ..., G.
* Error messages are consistent (no call, actionable hints).

## Bug fixes

* `dif_test()` and `tppcmtree(orthogonal = TRUE)` no longer fail with "system
  is computationally singular" when case-wise scores are (nearly) linearly
  dependent, e.g. when a step discrimination is estimated at 0. Null
  directions are dropped and `npar` reports the degrees of freedom used.
* `dif_test()` no longer lists the latent means and variances of a multigroup
  fit as (untested) units, and splitting all slopes or all locations no
  longer triggers a spurious "degenerate scores" message.
* `print()` of tables keeps group names such as `"36"` (was `"X36"`).
* Recipes: `irt_pars()$SE` (returns `NULL`) corrected to `$se`.

## Documentation

* `?tppcm` and `?tppcm-package` explain the log-scale item x step view of the
  five models (the interaction means that items differ in the shape of their
  step discriminations), the analogy to TAM/ConQuest location designs
  (`~ item + step` = RSM, `~ item + item:step` = PCM), and the parameter
  names (`disc`/`slope` = CBD, `diff` = CBL, `int` = IBD).
* `?dif_test` leads with the per-step DSF screen and documents why the 1-df
  `disc` rows of `by = "param"` are not clean CBD tests (location DSF leaks
  into them).
* The four vignettes have a design-formula section and a DSF workflow
  (screen with `dif_test(by = "step")`, decompose with `dsf()`,
  configural/metric/scalar with `dsf_design()`), with simulation evidence.

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
