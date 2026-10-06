#' Components of a fitted model
#'
#' Extracts the mathematical components behind the package functions from a
#' fitted model, in the spirit of `lme4::getME()`, `lavaan::lavInspect()` and
#' `mirt::extract.mirt()`. Without `what` it returns all components as one
#' object that every function of the package accepts in place of the fit, so
#' they are computed only once. Nothing is re-estimated: everything follows
#' from the fitted item response functions, the posterior and the responses
#' stored in the fit.
#'
#' Supported fits (the model is detected automatically):
#'
#' | Fit | `model` | `formula` |
#' |---|---|---|
#' | [TAM::tam.mml.3pl()] with [tppcm()] | `"E"` (linear step design; fixed steps detected) | detected from `E`: `~ 1`, `~ item`, `~ step`, `~ item * step`, or `NULL` (custom `index`) |
#' | [TAM::tam.mml.2pl()], `irtmodel = "2PL"` | `"tppcm"` (saturated TPPCM = NRM) | `~ item * step` |
#' | [TAM::tam.mml.2pl()], `irtmodel = "GPCM"` | `"gpcm"` | `~ item` |
#' | [TAM::tam.mml()] | `"pcm"` (free latent variance) | `~ 1` |
#' | [xxirt_tppcm()] | `"tppcm"`, `"gpcm"`, `"rank1"` | its `design` |
#'
#' Components (`what`):
#'
#' | `what` | Content |
#' |---|---|
#' | `"disc"`, `"diff"`, `"int"` | items x steps matrices: discriminations, difficulties (`disc * (theta - diff)`), intercepts (`slope * theta + int`) |
#' | `"est"`, `"se"` | slope parameters of the fitted model (e.g. item slopes, \eqn{\alpha, \gamma}) and their standard errors |
#' | `"vcov"` | covariance of the fitted model's parameters |
#' | `"par"`, `"partable"` | saturated parameters (per item \eqn{a_1..a_L, d_1..d_L}, \eqn{d = a \cdot diff}) and their description |
#' | `"vcov.sat"`, `"se.sat"` | implied covariance and standard errors of the saturated parameters |
#' | `"information"`, `"information.opg"` | observed (Louis) and outer-product information of the saturated parameters |
#' | `"gradient"`, `"scores"` | score vector and case-wise scores of the saturated parameters |
#' | `"estfun"`, `"estfun.si"` | case-wise scores of the fitted slopes and the step difficulties (`disc`, `diff`) or step intercepts (`slope`, `int`); used by [dif_test()] and [tppcmtree()] |
#' | `"jacobian"`, `"curvature"` | Jacobian of the fitted model in the saturated space and its curvature term |
#' | `"posterior"`, `"nodes"`, `"prior"`, `"probs"` | posterior, quadrature nodes, latent distribution at the nodes, item response probabilities |
#' | `"data"`, `"nobs"`, `"loglik"`, `"npar"` | responses, number of persons, log-likelihood, number of parameters |
#' | `"model"`, `"formula"` | type of fit, and its design formula (see [tppcm()]; `NULL` for a custom design) |
#' | `"groups"`, `"design"` | latent means and variances, design array `E` |
#'
#' The fitted model is a map \eqn{h(\vartheta)} into the saturated
#' parameters with Jacobian \eqn{J}; its observed information is
#' \eqn{J'IJ - \sum_k g_k \nabla^2 h_k}, where \eqn{I} is the Louis
#' information and \eqn{g} the score of the saturated model. The second term
#' is non-zero only for the product form. The latent trait distribution is
#' held fixed at its estimate.
#'
#' Any unidimensional fit whose item response functions are of the
#' divide-by-total form \eqn{\log P_k / P_{k-1} = a_k\theta - d_k} is accepted,
#' including binary items (one step), several groups, person weights
#' (`pweights`; [m2()] and [dif_test()] then refuse) and latent regression
#' (with a warning; regression coefficients are treated as fixed, [m2()]
#' refuses). Fits with guessing parameters, many-facet models
#' (`tam.mml.mfr()`), JML fits and multidimensional fits are refused.
#'
#' For TAM fits the location parameters follow TAM's design matrix `A`, so
#' restricted thresholds (e.g. `tam.mml(irtmodel = "RSM")`) and fixed `xsi`
#' are taken into account. Multidimensional fits are not supported. With
#' several groups (`group = ` in TAM), the free latent means and variances of
#' groups 2, ..., G are part of the parameter vector: their scores and
#' information are included, the item standard errors account for them, and
#' `what = "groups"` reports them with standard errors. Treating them as fixed
#' would understate the item standard errors by up to 40% on real data.
#'
#' @param x A fitted model (see table) or the object returned by
#'   `get_parts()`.
#' @param what `NULL` (all components) or one or more names from the table;
#'   several names give a named list.
#' @param model Fitted model type; detected from `x` when `NULL`. One of
#'   `"E"`, `"tppcm"`, `"gpcm"`, `"pcm"`, `"pcm_fixed"`, `"rank1"`.
#' @param fixed For `tam.mml.3pl` fits: indices of fixed slope parameters
#'   (default: those with `se.gammaslope == 0`).
#' @return Without `what`, an object of class `tppcm_parts` with methods for
#'   `print()`, `coef()` (the slope parameters of the fitted model, as
#'   `what = "est"`), `vcov()`, `logLik()` and `nobs()`; otherwise the
#'   component(s). The internal elements of the object are not the components
#'   of the table: use `get_parts(x, what)`, not `x$what`.
#' @examples
#' \donttest{
#' data(data.Students, package = "CDM")
#' sc <- data.Students[, c("sc1", "sc2", "sc3", "sc4")]
#' fit <- tam.mml.3pl(sc, E = tppcm(sc), est.variance = FALSE, verbose = FALSE)
#' get_parts(fit, "disc")                 # step discriminations
#' get_parts(fit, c("loglik", "npar"))
#' x <- get_parts(fit)                    # all components, computed once
#' dim(get_parts(x, "information"))
#' irt_pars(x); wald_test(x)            # reuse without recomputing
#' }
#' @export
get_parts <- function(x, what = NULL, model = NULL, fixed = NULL) {
  x <- if (inherits(x, "tppcm_parts")) x else .parts(x, model, fixed)
  if (is.null(what)) return(x)
  out <- lapply(what, function(w) .part(x, w))
  if (length(what) == 1) out[[1]] else stats::setNames(out, what)
}

.part_names <- c("disc", "diff", "int", "est", "se", "vcov", "par", "partable", "vcov.sat",
                 "se.sat", "information", "information.opg", "gradient", "scores", "estfun",
                 "estfun.si", "jacobian",
                 "curvature", "posterior", "nodes", "prior", "probs", "data", "nobs", "loglik",
                 "npar", "model", "formula", "groups", "design")

.part <- function(x, what) {
  if (!what %in% .part_names) stop("unknown component '", what, "'; see ?get_parts", call. = FALSE)
  w <- what
  r <- x$restricted
  switch(w,
    disc = x$a, diff = x$diff, int = x$int,
    est = stats::setNames(r$est, r$par), se = stats::setNames(r$se, r$par),
    vcov = x$vcov, par = x$par, partable = x$parindex, vcov.sat = x$vcov_sat,
    se.sat = sqrt(pmax(diag(x$vcov_sat), 0)),
    information = x$info, information.opg = x$info_opg, gradient = x$grad,
    scores = x$scores, estfun = .step_scores(x, "irt"), estfun.si = .step_scores(x, "si"),
    jacobian = x$J, curvature = x$curv,
    posterior = x$post, nodes = x$theta, prior = x$weights, probs = x$probs,
    data = x$data, nobs = x$nobs, loglik = x$loglik, npar = .npar(x),
    model = x$model, formula = .design_formula(.fit_design(x)), groups = x$groups, design = x$E)
}

.npar <- function(x) ncol(x$J)      # includes the free group means and variances

#' @export
vcov.tppcm_parts <- function(object, ...) object$vcov

#' @export
coef.tppcm_parts <- function(object, ...) .part(object, "est")

#' @export
logLik.tppcm_parts <- function(object, ...) {
  structure(object$loglik, df = .npar(object), nobs = object$nobs, class = "logLik")
}

#' @export
nobs.tppcm_parts <- function(object, ...) object$nobs

.parts <- function(fit, model = NULL, fixed = NULL) {
  if (is.null(model)) model <- .detect_model(fit)
  dat <- as.matrix(if (inherits(fit, "xxirt")) fit$dat else fit$resp)
  items <- .item_names(dat)
  ncat <- apply(dat, 2, max, na.rm = TRUE) + 1
  post  <- CDM::IRT.posterior(fit)
  probs <- CDM::IRT.irfprob(fit)
  if (NCOL(attr(post, "theta")) > 1) stop("only unidimensional models are supported", call. = FALSE)
  theta <- attr(post, "theta")[, 1]
  par   <- .par_from_probs(probs, theta, ncat, items)
  wts   <- if (inherits(fit, "xxirt")) fit$weights else fit$pweights
  if (!is.null(wts) && all(abs(wts - 1) < 1e-12)) wts <- NULL
  n_grp <- if (!inherits(fit, "xxirt") && !is.null(fit$group)) length(unique(fit$group)) else 1
  latreg <- !inherits(fit, "xxirt") && !is.null(fit$Y) && NCOL(fit$Y) > n_grp
  if (latreg) warning("latent regression: the regression coefficients are treated as fixed ",
                      "in the standard errors and tests", call. = FALSE)
  groups <- .group_pars(fit)
  grp <- NULL
  if (n_grp > 1 && !latreg) {                      # free group means and variances (groups 2..G)
    lev <- sort(unique(fit$group))
    grp <- list(g = match(fit$group, lev), MEAN = groups$MEAN, VAR = groups$VAR, free = seq_len(n_grp)[-1])
  }
  core  <- .score_core(dat, post, probs, theta, ncat, weights = wts, grp = grp)
  pidx  <- core$parindex
  ng    <- length(core$gnames)

  if (any(apply(dat, 2, min, na.rm = TRUE) > 0))
    warning("category 0 is not observed for some items; were the data coded 1, ..., K?", call. = FALSE)
  E <- NULL
  if (model == "E" && !is.null(fit$variance) && abs(as.numeric(fit$variance)[1] - 1) > 1e-6)
    warning("the latent variance was estimated together with free slopes, so the scale is not ",
            "identified; refit with tam.mml.3pl(..., est.variance = FALSE)", call. = FALSE)
  if (model == "E") {
    E <- fit$E
    if (is.null(fixed)) fixed <- which(fit$se.gammaslope == 0)
  }
  loc  <- .loc_design(fit, pidx)
  J    <- .jacobian(model, pidx, par, E = E, fixed = fixed, loc = loc)
  curv <- .curvature(model, J, pidx, core$grad[seq_len(nrow(J))])
  if (ng) {                                        # group parameters: identity block
    J <- rbind(cbind(J, matrix(0, nrow(J), ng, dimnames = list(NULL, core$gnames))),
               cbind(matrix(0, ng, ncol(J)), diag(ng)))
    rownames(J) <- c(pidx$name, core$gnames)
    curv <- rbind(cbind(curv, matrix(0, nrow(curv), ng)), matrix(0, ng, ncol(curv) + ng))
    dimnames(curv) <- list(colnames(J), colnames(J))
  }
  V    <- solve(crossprod(J, core$info %*% J) - curv)
  Vsat <- J %*% V %*% t(J)
  p    <- nrow(pidx)
  if (ng) {
    groups$SE_MEAN <- groups$SE_VAR <- NA_real_
    se_g <- sqrt(diag(V))[core$gnames]
    groups$SE_MEAN[grp$free] <- se_g[paste0("grp:mean", grp$free)]
    groups$SE_VAR[grp$free]  <- se_g[paste0("grp:var", grp$free)]
  }

  ia <- which(pidx$type == "a"); id <- which(pidx$type == "d")
  L  <- max(ncat) - 1
  to_mat <- function(v, vals = par[v]) {
    M <- matrix(NA_real_, length(items), L, dimnames = list(items, paste0("step", 1:L)))
    M[cbind(pidx$itemnr[v], pidx$step[v])] <- vals
    M
  }
  sc <- .slope_cols(J, pidx)
  est <- if (model == "rank1") {
    A <- to_mat(ia); c(A[, 1], colMeans(A / A[, 1])[-1])
  } else if (length(sc)) qr.solve(J[ia, sc, drop = FALSE], par[ia]) else numeric(0)

  structure(list(model = model, items = items, ncat = ncat,
                 a = to_mat(ia), diff = to_mat(id, par[id] / par[ia]), int = to_mat(id, -par[id]),
                 par = par, parindex = pidx,
                 theta = theta, post = post, probs = probs,
                 grad = core$grad, scores = core$scores, info = core$info,
                 info_opg = core$info_opg, nobs = core$nobs, gnames = core$gnames,
                 J = J, curv = curv, vcov = V, vcov_sat = Vsat[seq_len(p), seq_len(p)],
                 restricted = data.frame(par = sc, est = as.numeric(est),
                                         se = sqrt(diag(V))[sc], row.names = NULL),
                 E = E, fixed = fixed, n_loc = sum(grepl("^loc:", colnames(J))),
                 groups = groups, grp = grp, data = dat,
                 weights = .prior_weights(fit, post, theta),
                 weighted = !is.null(wts), latreg = latreg,
                 loglik = tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NA_real_)),
            class = "tppcm_parts")
}

# Location design of a TAM fit: d_il as linear functions of the free xsi
.loc_design <- function(fit, pidx) {
  if (!inherits(fit, c("tam.mml", "tam.mml.3pl")) || is.null(fit$A)) return(NULL)
  A <- fit$A; A[is.na(A)] <- 0
  M <- matrix(0, nrow(pidx), dim(A)[3])
  for (r in which(pidx$type == "d")) {
    i <- pidx$itemnr[r]; l <- pidx$step[r]
    M[r, ] <- A[i, l + 1, ] - A[i, l, ]
  }
  if (!is.null(fit$xsi.fixed)) M <- M[, -fit$xsi.fixed[, 1], drop = FALSE]
  M <- M[, colSums(abs(M)) > 0, drop = FALSE]
  q <- qr(M)
  M <- M[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  if (ncol(M) == sum(pidx$type == "d")) return(NULL)   # all thresholds free
  colnames(M) <- paste0("loc:", seq_len(ncol(M)))
  M
}

# Latent distribution at the nodes used in the likelihood. TAM's "prob.theta"
# (pi.k) is the mean posterior, not the normal prior of the fitted model.
.prior_weights <- function(fit, post, theta) {
  if (inherits(fit, c("tam.mml", "tam.mml.3pl"))) {
    g <- .group_pars(fit)
    w <- stats::dnorm(theta, g$MEAN[1], sqrt(g$VAR[1]))
  } else w <- as.numeric(as.matrix(attr(post, "prob.theta"))[, 1])
  w / sum(w)
}

# Columns of J that belong to slope parameters
.slope_cols <- function(J, pidx) {
  setdiff(colnames(J), c(pidx$name[pidx$type == "d"], grep("^(loc|grp):", colnames(J), value = TRUE)))
}

# Latent means and variances per group
.group_pars <- function(fit) {
  if (inherits(fit, "xxirt")) {
    return(data.frame(group = 1, MEAN = unname(fit$customTheta$par[1]),
                      VAR = unname(fit$customTheta$par[2])^2))
  }
  if (is.null(fit$beta)) return(NULL)
  g <- if (is.null(fit$group)) rep(1, nrow(as.matrix(fit$resp))) else fit$group
  lev <- sort(unique(g))
  v <- as.numeric(fit$variance)
  var_g <- if (length(v) == length(g)) v[match(lev, g)] else rep(v, length.out = length(lev))
  lab <- if (!is.null(fit$groups) && length(fit$groups) == length(lev)) fit$groups else lev   # TAM recodes group to 1..G
  data.frame(group = lab, MEAN = as.numeric(as.matrix(fit$beta)[seq_along(lev), 1]), VAR = var_g)
}

#' @export
print.tppcm_parts <- function(x, digits = 3, ...) {
  cat("Components of a fitted model:", .label(x), "|", length(x$items), "items,",
      max(x$ncat), "categories,", x$nobs, "persons\n")
  cat("get_parts(x, what = ) with what in:\n")
  cat(strwrap(paste(.part_names, collapse = ", "), indent = 2, exdent = 2), sep = "\n")
  invisible(x)
}

.label <- function(x) {
  code <- .fit_design(x)
  lab <- if (code %in% names(.design_codes)) sprintf("%s, design %s", .design_names[[code]], .design_codes[[code]])
         else switch(code, custom = "custom step design", fixed = "fixed-slope design",
                     pcm_fixed = "PCM (fixed slopes and variance)")
  if (x$model == "pcm") lab <- paste(lab, "(TAM: slope 1, free variance)")
  if (x$model == "E") {
    if (isTRUE(attr(code, "fixed"))) lab <- paste(lab, "with fixed steps")
    lab <- paste(lab, "(tam.mml.3pl)")
  }
  if (isTRUE(x$n_loc > 0)) lab <- sprintf("%s, %d free location parameters", lab, x$n_loc)
  lab
}

# Design of the fitted model: "pcm", "gpcm", "step", "rank1", "tppcm" (the
# design formulas of tppcm()), "custom", "fixed" or "pcm_fixed"; for
# tam.mml.3pl fits detected from the step design, with attribute "fixed" when
# some steps are fixed
.fit_design <- function(x) {
  if (x$model != "E") return(if (x$model == "pcm") "pcm" else x$model)
  pidx <- x$parindex; ia <- which(pidx$type == "a")
  A <- x$J[ia, .slope_cols(x$J, pidx), drop = FALSE]
  fixed <- any(rowSums(A != 0) == 0)
  if (!ncol(A)) return("fixed")
  per_item <- apply(A, 2, function(v) length(unique(pidx$itemnr[ia][v != 0])))
  per_step <- apply(A, 2, function(v) length(unique(pidx$step[ia][v != 0])))
  n_item <- length(unique(pidx$itemnr)); n_step <- max(pidx$step)
  code <- if (all(colSums(A != 0) == 1) && ncol(A) == sum(rowSums(A != 0) > 0)) "tppcm"
          else if (ncol(A) == 1) "pcm"
          else if (all(per_item == 1) && ncol(A) == n_item) "gpcm"
          else if (all(per_step == 1) && ncol(A) == n_step) "step"
          else "custom"
  structure(code, fixed = fixed)
}
