#' tppcm: Step discriminations for the two-parameter partial credit model
#'
#' The two-parameter partial credit model (TPPCM; Yu, 1991) gives every step
#' \eqn{l} of item \eqn{i} of Masters' partial credit model its own
#' discrimination \eqn{a_{il}}. With categories
#' \eqn{k = 0, \dots, K-1},
#' \deqn{P(X_i = k \mid \theta) \propto \exp\Big(\sum_{l \le k} (a_{il}\theta - d_{il})\Big).}
#' The saturated TPPCM is a reparameterization of the nominal response model,
#' so the package does not ship an estimator of its own. It builds on:
#'
#' * **TAM** for estimation: [TAM::tam.mml.3pl()] with the step design from
#'   [step_design()] (the slope parameters are then the step discriminations),
#'   and [TAM::tam.mml.2pl()] / [TAM::tam.mml()] for the GPCM and the PCM.
#' * **sirt** for what TAM cannot express exactly: the product form
#'   \eqn{a_{il} = \alpha_i \gamma_l} and the bound \eqn{a_{il} \ge 0}
#'   ([xxirt_tppcm()]).
#'
#' The functions below take the fitted TAM (or xxirt) object directly; the
#' responses are read from the fit. The underlying mathematical objects are
#' available with [get_parts()].
#'
#' * [irt_pars()] with `long = TRUE`: step discriminations and difficulties with standard
#'   errors from the full observed information (TAM's `se.gammaslope` ignores
#'   the covariance with the thresholds and is too small);
#' * [wald_test()]: per-item test of equal step discriminations;
#' * [irt_pars()]: a compact, mirt-style table of item and group parameters
#'   with standard errors and a fit header;
#' * [m2()]: limited-information fit (M2, RMSEA2, SRMSR, CFI, TLI) for any
#'   step design;
#' * [score_test()], [mi()], [item_test()]: score tests, modification indices
#'   and expected parameter changes, computed from the restricted model only;
#' * [dif_test()], [tppcmtree()]: score-based DIF tests and DIF trees for
#'   step discriminations and step difficulties separately.
#'
#' * [step_info()]: step, item and test information.
#'
#' The vignettes `vignette("tppcm-tutorial")` and `vignette("tppcm-recipes")`
#' (also in Traditional Chinese: `"tppcm-tutorial-zh"`, `"tppcm-recipes-zh"`)
#' show every function on data sets from TAM, CDM, mirt and psychotools.
#'
#' @references Yu, M.-N. (1991). *A two-parameter partial credit model*
#'   (Doctoral dissertation). University of Illinois at Urbana-Champaign.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats pchisq qnorm lm coef dnorm rnorm runif setNames logLik nobs vcov
#' @importFrom TAM tam.mml.3pl
## usethis namespace: end
NULL
