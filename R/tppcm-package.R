#' tppcm: Boundary-level diagnostics for the two-parameter partial credit model
#'
#' The two-parameter partial credit model (TPPCM; Yu, 1991) gives every step
#' \eqn{l} of item \eqn{i} of Masters' partial credit model its own
#' discrimination \eqn{a_{il}}. With categories
#' \eqn{k = 0, \dots, K-1},
#' \deqn{P(X_i = k \mid \theta) \propto \exp\Big(\sum_{l \le k} a_{il}(\theta - b_{il})\Big).}
#' The saturated TPPCM is a reparameterization of the nominal response model,
#' so the package does not ship an estimator of its own: it provides
#' post-estimation, boundary-level diagnostics for fits from TAM and sirt.
#'
#' Estimation is left to
#'
#' * **TAM**: [TAM::tam.mml.3pl()] with the step design from [tppcm()] (the
#'   slope parameters are then the step discriminations), and
#'   [TAM::tam.mml.2pl()] / [TAM::tam.mml()] for the GPCM and the PCM;
#' * **sirt**: the product form \eqn{a_{il} = \alpha_i \gamma_l} and the
#'   bound \eqn{a_{il} \ge 0}, which TAM cannot express ([xxirt_tppcm()]).
#'
#' The functions take the fitted TAM (or xxirt) object directly; the
#' responses are read from the fit.
#'
#' @section Parameter names:
#' Every step (category boundary) has a discrimination and a location. The
#' package uses two sets of names, one for reporting single-group fits and
#' one, from the category boundary literature (Preston et al., 2011), for
#' differential step functioning (DSF) across groups:
#'
#' | Quantity | Kernel | [irt_pars()], [dif_test()] (`param = "irt"`) | `param = "si"` / `IRTpars = FALSE` | [dsf()], [dsf_design()] |
#' |---|---|---|---|---|
#' | step discrimination, CBD (category boundary discrimination) | \eqn{a_{il}} | `disc` | `slope` | `"cbd"` |
#' | step difficulty, CBL (category boundary location) | \eqn{b_{il}} | `diff` | | `"cbl"` |
#' | step intercept, IBD (intercept boundary difference) | \eqn{-a_{il} b_{il}} | | `int` | `"ibd"` |
#'
#' The step kernel is `disc * (theta - diff)` = `slope * theta + int`. The CBD
#' equals the difference of adjacent category slopes of the nominal response
#' model, and the IBD the difference of adjacent category intercepts. In
#' multigroup comparisons the CBD plays the role of a factor loading and the
#' IBD that of an intercept (lavaan's `group.equal = c("loadings",
#' "intercepts")`); with equal CBDs, equal IBDs and equal CBLs are the same
#' restriction.
#'
#' @section Submodels (an item x step design on the log scale):
#' The submodels are written as formulas for the step discriminations of the
#' items, read as a two-way item x step layout on the log scale,
#' \eqn{\log a_{il} = \mu + \alpha_i + \gamma_l + (\alpha \gamma)_{il}}
#' (argument `design` of [tppcm()], [xxirt_tppcm()], [score_test()],
#' [dsf_design()] and [tppcmtree()]):
#'
#' | `design` | Model | \eqn{a_{il}} | Shape of the steps |
#' |---|---|---|---|
#' | `~ 1` | PCM (common slope) | \eqn{a} | flat, same level for all items |
#' | `~ item` | GPCM | \eqn{\alpha_i} | flat, item-specific level |
#' | `~ step` | step model | \eqn{\gamma_l} | same shape and level for all items |
#' | `~ item + step` | Rank-1 (product form; [xxirt_tppcm()]) | \eqn{\alpha_i\gamma_l} | same shape, item-specific level |
#' | `~ item * step` | saturated TPPCM (default) | \eqn{a_{il}} | item-specific shapes (item x step interaction) |
#'
#' The models are nested from top to bottom, except that `~ item` and
#' `~ step` do not contain each other. The score test of the GPCM against
#' `~ item + step` asks whether there is a step effect, that of Rank-1
#' against `~ item * step` whether there is an item x step interaction (items
#' differ in the ratios \eqn{a_{il}/a_{il'}}). The log scale presumes
#' \eqn{a_{il} > 0}, as in Yu (1991); TAM does not impose it. In the same
#' language, TAM's and ConQuest's location formulas `~ item + step` (rating
#' scale model) and `~ item + item:step` (partial credit model) describe the
#' step locations; `design = ~ step` is the rating-scale idea applied to the
#' discriminations.
#'
#' @section Functions:
#' * Design: [tppcm()] (`design = ~ item * step` and its submodels) for [TAM::tam.mml.3pl()];
#'   [xxirt_tppcm()] for the product form and bounds; [sim_tppcm()] for
#'   simulation.
#' * Parameters: [irt_pars()], a compact table of item and group parameters
#'   with standard errors from the full observed information (TAM's
#'   `se.gammaslope` ignores the covariance with the thresholds) and a fit
#'   header; `long = TRUE` for a long table with confidence intervals.
#' * Fit and information: [m2()] (M2, RMSEA2, SRMSR, CFI, TLI for any step
#'   design); [step_info()] (step, item and test information).
#' * Tests between models: [wald_test()] (equal step discriminations within
#'   items); [score_test()], [item_test()], [mi()] (score tests, modification
#'   indices and expected parameter changes, computed from the restricted
#'   model only).
#' * Differential step functioning: [dif_test()] (score tests from one fit;
#'   `by = "step"` is the per-boundary screen on a multigroup fit); [dsf()]
#'   (likelihood ratio tests of single CBDs and IBDs); [dsf_design()]
#'   (configural, metric, scalar and partial multigroup models);
#'   [tppcmtree()] (DIF trees for unknown groupings).
#' * Components: [get_parts()] (scores, information, Jacobian, posterior, ...).
#'
#' Tables of tests report `p` and an adjusted `p_adj` (argument `adjust`,
#' Holm by default); the header of the printed table says what was tested
#' and how `p_adj` was adjusted.
#'
#' The vignettes `vignette("tppcm-tutorial")` and `vignette("tppcm-recipes")`
#' (also in Traditional Chinese: `"tppcm-tutorial-zh"`, `"tppcm-recipes-zh"`)
#' show every function on data sets from TAM, CDM, mirt and psychotools.
#'
#' @references Preston, K., Reise, S., Cai, L., & Hays, R. D. (2011). Using
#'   the nominal response model to evaluate response category discrimination
#'   in the PROMIS emotional distress item pools. *Educational and
#'   Psychological Measurement, 71*, 523-550.
#'
#'   Yu, M.-N. (1991). *A two-parameter partial credit model*
#'   (Doctoral dissertation). University of Illinois at Urbana-Champaign.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats pchisq qnorm lm coef dnorm rnorm runif setNames logLik nobs vcov
#' @importFrom TAM tam.mml.3pl
## usethis namespace: end
NULL
