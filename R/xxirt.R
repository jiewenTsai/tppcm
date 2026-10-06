#' Fit TPPCM-family models with sirt::xxirt()
#'
#' Use this for what TAM cannot estimate exactly: Rank-1, the multiplicative
#' product form (`design = ~ item + step`: \eqn{a_{il} = \alpha_i\gamma_l},
#' \eqn{\gamma_1 = 1}, \eqn{\gamma} shared across items; see [tppcm()] for
#' the item x step view of the designs) and the bound \eqn{a_{il} \ge 0}. The
#' saturated model (`~ item * step`) and the GPCM (`~ item`) are also
#' available, e.g. for simulations; fit `~ 1` and `~ step` with
#' [TAM::tam.mml.3pl()] and [tppcm()].
#'
#' The latent trait is \eqn{N(0, 1)}. The default nodes are those of
#' TAM (`seq(-6, 6, len = 21)`), so log-likelihoods can be compared with TAM
#' fits of the same data (e.g. with `CDM::IRT.compareModels()`; use its
#' `LRtest` part, not `anova()`).
#'
#' @param dat Item responses coded `0, ..., K-1`, same `K` for all items.
#' @param design `~ item * step` (default, saturated), `~ item + step`
#'   (Rank-1) or `~ item` (GPCM); see [tppcm()].
#' @param Theta Quadrature nodes (matrix with one column).
#' @param lower0 Impose \eqn{a_{il} \ge 0} (and \eqn{\alpha, \gamma \ge 0}).
#' @param ... Passed to [sirt::xxirt()] (e.g. `conv`, `maxit`).
#' Item parameters are stored in slope-intercept form: step `l` has kernel
#' `a_l * theta + int_l` (`alpha * g_l * theta + int_l` for the product form).
#' Use [irt_pars()] for discriminations and difficulties.
#' @return An `xxirt` object (with the additional class `xxirt_tppcm`, whose
#'   `print()` method shows a short summary instead of the whole list). Use
#'   [xxirt_vcov()] for its covariance matrix and [irt_pars()], [score_test()]
#'   and the other tools of this package.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(500, disc = outer(c(0.8, 1, 1.2, 1.5), c(1, 1.5, 0.7)),
#'                       diff = matrix(c(-1, 0, 1), 4, 3, byrow = TRUE))
#' fit <- xxirt_tppcm(dat, design = ~ item + step)       # Rank-1
#' fit
#' irt_pars(fit, restricted = TRUE)                       # alpha_i, gamma_l
#' score_test(fit)                                        # item x step interaction?
#' }
#' @export
xxirt_tppcm <- function(dat, design = ~ item * step,
                        Theta = matrix(seq(-6, 6, length.out = 21), ncol = 1),
                        lower0 = TRUE, ...) {
  if (!is.null(list(...)$model))      # would otherwise reach sirt::xxirt()
    stop("xxirt_tppcm() has no 'model' argument; use design = ~ item + step (Rank-1), ~ item * step ",
         "or ~ item", call. = FALSE)
  model <- .parse_design(design)
  if (!model %in% c("tppcm", "rank1", "gpcm"))
    stop("xxirt_tppcm() fits design = ~ item * step, ~ item + step or ~ item; fit the linear design ",
         .design_codes[[model]], " with tam.mml.3pl(dat, E = tppcm(dat, design = ",
         .design_codes[[model]], "), est.variance = FALSE)", call. = FALSE)
  dat <- as.data.frame(dat)
  ncat <- apply(dat, 2, max, na.rm = TRUE) + 1
  if (length(unique(ncat)) > 1) stop("all items must have the same number of categories", call. = FALSE)
  K <- unname(ncat[1]); L <- K - 1; I <- ncol(dat); items <- colnames(dat)
  Theta <- as.matrix(Theta)
  lo <- if (lower0) 0 else -Inf
  int0 <- seq(1, -1, length.out = L)
  P_sat <- function(par, Theta, ncat) .tppcm_P(par[1:L], -par[L + 1:L], Theta)   # slope*theta + int
  P_r1  <- function(par, Theta, ncat) .tppcm_P(par[1] * par[1 + 1:L], -par[1 + L + 1:L], Theta)
  items_def <- list(
    sirt::xxirt_createDiscItem(
      name = "TPPCM", P = P_sat, est = rep(TRUE, 2 * L),
      par = stats::setNames(c(rep(1, L), int0), c(paste0("a", 1:L), paste0("int", 1:L))),
      lower = c(rep(lo, L), rep(-Inf, L))),
    sirt::xxirt_createDiscItem(
      name = "RANK1TPPCM", P = P_r1, est = rep(TRUE, 1 + 2 * L),
      par = stats::setNames(c(1, rep(1, L), int0), c("alpha", paste0("g", 1:L), paste0("int", 1:L))),
      lower = c(lo, rep(lo, L), rep(-Inf, L))))
  theta_def <- sirt::xxirt_createThetaDistribution(
    par = c(mu = 0, sigma = 1), est = c(FALSE, FALSE),
    P = function(par, Theta, G) {
      p <- stats::dnorm(Theta[, 1], mean = par[1], sd = par[2])
      matrix(p / sum(p), ncol = 1)
    })
  if (model == "rank1") {
    pt <- sirt::xxirt_createParTable(dat, itemtype = rep("RANK1TPPCM", I), customItems = items_def)
    pt <- sirt::xxirt_modifyParTable(pt, parname = "g1", value = 1, est = FALSE)
    idx0 <- max(pt$parindex, na.rm = TRUE)
    for (l in 2:L) pt <- sirt::xxirt_modifyParTable(pt, parname = paste0("g", l), item = items,
                                                    parindex = idx0 + l - 1)
  } else {
    pt <- sirt::xxirt_createParTable(dat, itemtype = rep("TPPCM", I), customItems = items_def)
    if (model == "gpcm") {
      idx0 <- max(pt$parindex, na.rm = TRUE)
      for (i in 1:I) for (nm in paste0("a", 1:L))
        pt <- sirt::xxirt_modifyParTable(pt, parname = nm, item = items[i], parindex = idx0 + i)
    }
  }
  fit <- sirt::xxirt(dat = dat, Theta = Theta, partable = pt, customItems = items_def,
                     customTheta = theta_def, verbose = FALSE, ...)
  fit$tppcm_model <- model
  fit$tppcm_design <- .design_formula(model)
  class(fit) <- c("xxirt_tppcm", class(fit))
  fit
}

#' @export
print.xxirt_tppcm <- function(x, ...) {
  cat(sprintf("xxirt fit of the %s (design %s): %d items, %d persons, logLik = %.2f, %d parameters\n",
              .design_names[[x$tppcm_model]], .design_codes[[x$tppcm_model]], ncol(x$dat), nrow(x$dat), as.numeric(x$loglik),
              length(x$par_items)))
  cat("irt_pars(x) for the item parameters, xxirt_vcov(x) for the covariance matrix,\n",
      "summary(x) for sirt's summary\n", sep = "")
  invisible(x)
}

#' Covariance matrix of an xxirt fit with shared or fixed parameters
#'
#' `vcov()` for `xxirt` objects (sirt 4.2-133) fails with
#' "the condition has length > 1" when parameters are shared across rows of the
#' parameter table or fixed. This recomputes the numerical Hessian without the
#' shortcut that causes the error.
#'
#' @param x An `xxirt` object.
#' @return Inverse of the negative Hessian of the log-likelihood.
#' @examples
#' \donttest{
#' data(data.Students, package = "CDM")
#' sc <- as.matrix(data.Students[, c("sc1", "sc2", "sc3", "sc4")])
#' fit <- xxirt_tppcm(sc, design = ~ item + step)
#' fit
#' round(sqrt(diag(xxirt_vcov(fit))), 3)
#' }
#' @export
xxirt_vcov <- function(x) {
  solve(-sirt::xxirt_hessian(x, use_shortcut = FALSE))
}
