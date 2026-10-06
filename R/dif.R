#' Score (LM) tests of differential step and item functioning
#'
#' Tests whether item parameters differ between groups, or vary with a person
#' covariate, using only the fitted model (no refitting). The main use is the
#' per-step DSF screen on a multigroup fit:
#' ```
#' fit <- tam.mml.3pl(dat, E = tppcm(dat), group = g, est.variance = FALSE)
#' dif_test(fit, g, by = "step")
#' ```
#' which tests, for every category boundary, its discrimination (CBD) and
#' location jointly across the groups: the Lagrange multiplier (score) test
#' for DIF of Glas (1998), with 2 df per group contrast. The units of flagged
#' steps (`"item:step"`) can be passed to `dsf(fit, items = )` as they are, to
#' tell CBD from location DSF by likelihood ratio tests. See [tppcm-package]
#' for the parameter names (CBD, CBL, IBD).
#'
#' **Multigroup fits (recommended).** Fit the model with the grouping
#' (`tam.mml.3pl(..., group = g, est.variance = FALSE)`): the latent means and
#' variances of groups 2, ..., G are then estimated and only item-level
#' differences remain in the scores. For a factor covariate the test is the
#' score test of the model in which the tested parameters are split by the
#' levels of the covariate, with `npar * (levels - 1)` degrees of freedom.
#' By default (`info = "robust"`) it is the generalized score test of Boos
#' (1992): the observed (Louis) information of each level as the Hessian and
#' the outer product of the case-wise scores in the middle, which keeps its
#' level when the model does not fit exactly (as with most rating-scale data;
#' Falk & Monroe, 2018; Guastadisegni et al., 2021). In simulations of
#' per-step tests (two groups, 300 or 600 persons per group, with and without
#' local dependence) the robust test kept its level in every condition; it is
#' slightly conservative with 300 per group (about .03), and with the size
#' equalized it was as powerful as the likelihood ratio test. `info =
#' "louis"` (model-based) agrees closely with the likelihood ratio test
#' without refitting and shares its slight excess of rejections in small
#' samples; `info = "opg"` is liberal and is for comparison only.
#'
#' **Which unit.** `by = "step"` is the screen: the joint test of a step has a
#' clean null (the boundary is the same in all groups) and does not depend on
#' the parameterization (CBD and CBL, or CBD and IBD, give the same
#' statistic); it is computed in the slope-intercept form, so a step with
#' discrimination near 0 keeps both degrees of freedom. With `by = "param"`
#' the 1-df `disc` row of a step tests its CBD with its location (`diff`, the
#' CBL) held equal across groups (with `param = "si"`: its intercept, the
#' IBD), so location DSF leaks into it: in a simulation with 600 persons per
#' group, a CBL shift of 0.5 in one step was flagged by its `disc` row in 16.5%
#' of the replications. Use [dsf()] to decompose a flagged step instead.
#' `by = "item"` and `"block"` test all selected parameters of an item or of
#' the test jointly.
#'
#' **Impact is not DIF.** A single-group fit assumes one latent distribution
#' for everybody. When the covariate is related to the latent trait (groups
#' that differ in ability, a numeric covariate correlated with it), the scores
#' of all parameters move with the covariate and most items are flagged. Fit
#' the multigroup model (`group = `) or the latent regression (`Y = `) with
#' that covariate first and test on that fit.
#'
#' **Single-group fits and other covariates.** Without groups in the fit (or
#' for a numeric or ordered covariate) the case-wise scores
#' (`get_parts(fit, "estfun")`) are ordered by the covariate and aggregated
#' with a functional from \pkg{strucchange} (Merkle & Zeileis, 2013; Strobl,
#' Kopf & Zeileis, 2015). Each tested unit is replaced by its efficient scores
#' (residuals from the regression on all other parameters, outer-product
#' information), so DIF in the other parameters does not leak into the test.
#'
#' @param x A fitted model (see [get_parts()]) or a `get_parts` object.
#' @param covariate Person covariate in the order of the fitted responses
#'   (no missing values), usually the grouping of a multigroup fit. A factor
#'   (or character, logical) gives the LM test, an ordered factor the
#'   weighted maximum LM test (`"maxLMo"`), a numeric covariate the
#'   Cramer-von Mises type test (`"CvM"`).
#' @param parm Parameters to test: `c("disc", "diff")` (default, both
#'   blocks), `"disc"` (or `"slope"`; slope block), `"diff"` (or `"int"`;
#'   location block), or column names of `get_parts(fit, "estfun")`
#'   (`get_parts(fit, "estfun.si")` with `param = "si"`). `"all"` is accepted
#'   for `c("disc", "diff")`.
#' @param by Testing unit: every parameter (`"param"`), every step (`"step"`:
#'   discrimination and location jointly, units named `"item:step"`; needs
#'   both blocks and a TPPCM step design), every item (`"item"`) or all
#'   selected parameters together (`"block"`).
#' @param functional `NULL` (chosen from the covariate type), one of `"LM"`,
#'   `"maxLMo"`, `"DM"`, `"CvM"`, `"maxLM"`, or a \pkg{strucchange}
#'   functional. Not used by the split LM test of multigroup fits.
#' @param param Parameterization of the tested locations: `"irt"` (`disc`,
#'   `diff`: kernel `disc * (theta - diff)`) or `"si"` (`slope`, `int`: kernel
#'   `slope * theta + int`). `by = "step"` does not depend on it.
#' @param info Information used by the split LM test of multigroup fits:
#'   `"robust"` (default), `"louis"` or `"opg"`; see Details. Ignored
#'   otherwise.
#' @param adjust Multiple-testing adjustment of `p_adj` over all rows, a method
#'   of [stats::p.adjust()] (`"none"` for none).
#' @return Data frame (class `tppcm_table`) in the order of the items, with
#'   `unit` (parameter name, `"item:step"`, item name or `"all"`), `npar`
#'   (number of parameters tested), `statistic`, `df` (chi-square degrees of
#'   freedom of the LM tests, `npar * (levels - 1)`; `NA` for the other
#'   functionals), `p` and `p_adj`.
#' @seealso [dsf()] for likelihood ratio tests of single step parameters,
#'   [dsf_design()] for configural, metric and scalar models,
#'   [tppcmtree()] for unknown groupings.
#' @references Boos, D. D. (1992). On generalized score tests. *The American
#'   Statistician, 46*, 327-333.
#'
#'   Falk, C. F., & Monroe, S. (2018). On Lagrange multiplier tests in
#'   multidimensional item response theory: Information matrices and model
#'   misspecification. *Educational and Psychological Measurement, 78*,
#'   653-678.
#'
#'   Glas, C. A. W. (1998). Detection of differential item functioning using
#'   Lagrange multiplier tests. *Statistica Sinica, 8*, 647-667.
#'
#'   Guastadisegni, L., Cagnone, S., Moustaki, I., & Vasdekis, V. (2021). Use
#'   of the Lagrange multiplier test for assessing measurement invariance
#'   under model misspecification. *Educational and Psychological
#'   Measurement*. \doi{10.1177/00131644211020355}
#'
#'   Merkle, E. C., & Zeileis, A. (2013). Tests of measurement invariance
#'   without subgroups: A generalization of classical methods.
#'   *Psychometrika, 78*, 59-82.
#'
#'   Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method for
#'   detecting differential item functioning in the Rasch model.
#'   *Psychometrika, 80*, 289-316.
#' @examples
#' \donttest{
#' set.seed(1)
#' a <- matrix(c(1, 1.5, 1.2), 5, 3, byrow = TRUE)
#' b <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
#' a2 <- a; a2[3, 2] <- 2.5                       # DSF in the CBD of item 3, step 2
#' dat <- rbind(sim_tppcm(600, a, b),
#'              sim_tppcm(600, a2, b, theta = rnorm(600, 0.3, 0.8)))   # impact
#' g <- rep(c("ref", "foc"), each = 600)
#' fit <- tam.mml.3pl(dat, E = tppcm(dat), group = g, est.variance = FALSE,
#'                    verbose = FALSE)
#' st <- dif_test(fit, g, by = "step")            # screen: every step, 2 df
#' st
#' st$unit[st$p_adj < 0.05]                       # flagged steps, e.g. for dsf()
#' dif_test(fit, g, by = "item")                  # every item jointly
#' }
#' @export
dif_test <- function(x, covariate, parm = c("disc", "diff"), by = c("param", "step", "item", "block"),
                     functional = NULL, param = c("irt", "si"),
                     info = c("robust", "louis", "opg"), adjust = "holm") {
  if (!requireNamespace("strucchange", quietly = TRUE))
    stop("dif_test() needs the strucchange package: install.packages(\"strucchange\")", call. = FALSE)
  by <- match.arg(by); param <- match.arg(param); info <- match.arg(info)
  adjust <- .check_adjust(adjust)
  x <- get_parts(x)
  if (isTRUE(x$n_loc > 0))
    stop("dif_test() needs free step locations; RSM-type fits (restricted thresholds) are not supported",
         call. = FALSE)
  if (isTRUE(x$weighted)) stop("dif_test() does not support person weights (pweights)", call. = FALSE)
  blk_sel <- .parm_block(parm)
  if (by == "step") {
    if (!identical(blk_sel, "all"))
      stop('by = "step" tests the discrimination and location of a step jointly; use parm = c("disc", "diff")',
           call. = FALSE)
    param <- "si"                                     # joint test: invariant, no degenerate locations
  }
  S <- .step_scores(x, param)
  if (length(covariate) != nrow(S))
    stop("'covariate' must have one value per person of the fit (", nrow(S), "), not ", length(covariate),
         call. = FALSE)
  if (anyNA(covariate)) stop("'covariate' has missing values; drop those persons and refit", call. = FALSE)
  if (is.character(covariate) || is.logical(covariate)) covariate <- factor(covariate)
  if (is.factor(covariate) && (is.null(x$groups) || nrow(x$groups) == 1) && !isTRUE(x$latreg))
    message("single-group fit: a difference in ability between the groups of 'covariate' (impact) ",
            "is reported as DIF in the difficulties; refit with group = covariate if that is possible (see ?dif_test)")
  blk <- attr(S, "block")
  S <- sweep(S, 2, colMeans(S))                       # scores sum to zero at the MLE
  tested <- switch(if (is.null(blk_sel)) "names" else blk_sel,
    slope = which(blk == "slope"), location = which(blk == "location"),
    all = which(blk %in% c("slope", "location")),        # item parameters, not the group means/variances
    names = match(parm, colnames(S)))
  if (anyNA(tested))
    stop("unknown parameter(s) in 'parm': ", paste(parm[is.na(tested)], collapse = ", "),
         "; see colnames(get_parts(fit, \"", if (param == "si") "estfun.si" else "estfun", "\"))",
         call. = FALSE)
  nm <- colnames(S)[tested]
  it <- sub("_.*$", "", nm)
  if (by == "step") {
    st <- sub("^.*_(slope|int)([0-9]+)$", "\\2", nm)
    if (any(st == nm))
      stop('by = "step" needs one discrimination per step (a TPPCM step design, e.g. E = tppcm(dat))',
           call. = FALSE)
    su <- paste0(it, ":", st)
  }
  units <- switch(by,
    param = stats::setNames(as.list(tested), nm),
    step  = split(tested, factor(su, levels = unique(su[order(match(it, unique(it)), as.integer(st))]))),
    item  = split(tested, factor(it, levels = unique(it))),   # keep the item order
    block = list(all = tested))
  ref <- stats::median(sqrt(colSums(S^2)))           # scale for judging zero score columns
  fun_name <- if (is.character(functional) || is.null(functional)) .dif_functional_name(covariate, functional) else "user"
  # In a multigroup fit the latent means/variances of groups 2..G act only in
  # their own group. The Brownian-bridge functionals assume parameters common to
  # all persons, which makes the LM test conservative when the covariate is (or
  # overlaps with) the grouping; use the score test of the model in which the
  # tested parameters are split by the levels of the covariate.
  split_lm <- fun_name == "LM" && !is.null(x$groups) && nrow(x$groups) > 1
  n_lev <- if (is.factor(covariate)) nlevels(droplevels(covariate)) else NA_integer_
  if (split_lm) {
    zf <- droplevels(as.factor(covariate))
    Lz <- .louis_by_level(x, zf, param)
  }
  out <- do.call(rbind, lapply(names(units), function(u) {
    t <- units[[u]]
    if (split_lm) return(cbind(unit = u, .lm_split(Lz, S, t, zf, info)))
    St <- .full_rank_scores(.efficient(S, t), ref)
    if (ncol(St) == 0)
      return(data.frame(unit = u, npar = 0L, statistic = NA_real_, df = NA_integer_, p = NA_real_))
    gp <- strucchange::gefp(St, fit = NULL, scores = function(m) m, order.by = covariate)
    fn <- if (fun_name == "user") functional else .dif_functional(fun_name, gp)
    st <- strucchange::sctest(gp, functional = fn)
    data.frame(unit = u, npar = ncol(St), statistic = unname(st$statistic),
               df = if (fun_name == "LM") ncol(St) * (n_lev - 1L) else NA_integer_,
               p = unname(st$p.value))
  }))
  # expected testable directions: in a split test, splitting all slopes (all
  # locations) leaves the latent variances (means) of groups 2..G unidentified
  expected <- lengths(units)
  if (split_lm) expected <- expected - vapply(units, function(t)
    sum(all(which(blk == "slope") %in% t), all(which(blk == "location") %in% t)), 0)
  short <- out$npar < expected
  if (any(short))
    message("degenerate scores (e.g. the difficulty of a step whose discrimination is 0) left ",
            "fewer testable directions in: ", paste(out$unit[short], collapse = ", "),
            "; 'npar' and 'df' give what was tested",
            if (param == "irt") "; param = \"si\" avoids the problem for flat steps" else "")
  out <- .add_p_adj(out, adjust)
  rownames(out) <- NULL
  sl <- if (param == "irt") c("disc", "diff") else c("slope", "int")
  what <- switch(if (is.null(blk_sel)) "names" else blk_sel,
                 slope = sl[1], location = sl[2], all = paste(sl, collapse = " + "),
                 names = "selected parameters")
  unit <- switch(by, param = sprintf("parameter (%s)", what),
                 step = "step (discrimination and location jointly)",
                 item = sprintf("item (%s)", what), block = sprintf("all %s jointly", what))
  title <- if (split_lm)
    sprintf("LM tests of DIF between the %d groups of the covariate (%s information) | unit: %s | %s",
            nlevels(zf), info, unit, .adj_label(adjust))
  else sprintf("Score-based DIF tests along the covariate (%s functional) | unit: %s | %s",
               fun_name, unit, .adj_label(adjust))
  .as_table(out, title)
}

# Score test of H0: the parameters in columns t of S are equal across the levels
# of the factor z, against the model with one copy of them per level. The
# information of the augmented model is built from the observed (Louis)
# information of each level, mapped to the step parameterization with the same
# linear map as the scores; LM = U' I_aug^- U with U the level-wise score sums.
.lm_split <- function(Lz, S, t, z, info = c("robust", "louis", "opg")) {
  info <- match.arg(info)
  o <- setdiff(seq_len(ncol(S)), t); k <- length(t); G <- length(Lz)
  m <- k * G + length(o); io <- k * G + seq_along(o)
  Saug <- cbind(do.call(cbind, lapply(levels(z), function(l) S[, t, drop = FALSE] * (z == l))),
                S[, o, drop = FALSE])
  B <- crossprod(Saug)                                # OPG of the augmented model
  H <- if (info == "opg") B else {                    # observed (Louis) information
    Iaug <- matrix(0, m, m)
    for (l in seq_len(G)) {
      il <- (l - 1) * k + seq_len(k)
      Iaug[il, il] <- Lz[[l]][t, t]
      Iaug[il, io] <- Lz[[l]][t, o]; Iaug[io, il] <- Lz[[l]][o, t]
      Iaug[io, io] <- Iaug[io, io] + Lz[[l]][o, o]
    }
    Iaug
  }
  U <- colSums(Saug)
  Hr <- .psd_inv(H); Hi <- Hr$inv
  H0 <- if (info == "opg") crossprod(S) else Reduce(`+`, Lz)
  df <- Hr$rank - .psd_inv(H0)$rank                   # testable restrictions
  none <- data.frame(npar = 0L, statistic = NA_real_, df = NA_integer_, p = NA_real_)
  if (df <= 0) return(none)
  R <- cbind(kronecker(cbind(diag(G - 1), -1), diag(k)), matrix(0, k * (G - 1), length(o)))
  RHU <- R %*% Hi %*% U
  M <- if (info == "louis") R %*% Hi %*% t(R) else R %*% Hi %*% B %*% Hi %*% t(R)
  Mi <- .psd_inv(M)
  df <- min(df, Mi$rank)
  if (df <= 0) return(none)
  stat <- max(drop(crossprod(RHU, Mi$inv %*% RHU)), 0)
  data.frame(npar = df %/% (G - 1), statistic = stat, df = as.integer(df),
             p = stats::pchisq(stat, df, lower.tail = FALSE))
}

# Observed (Louis) information of each level of z in the step parameterization
.louis_by_level <- function(x, z, param) {
  P <- ncol(x$scores)
  e <- x; e$scores <- diag(P); colnames(e$scores) <- colnames(x$scores)
  Tm <- .step_scores(e, param)                       # saturated -> step map (rows: saturated)
  lapply(levels(z), function(l) {
    sel <- which(z == l)
    g <- x$grp; if (!is.null(g)) g$g <- g$g[sel]
    Ls <- .score_core(x$data[sel, , drop = FALSE], x$post[sel, , drop = FALSE], x$probs,
                      x$theta, x$ncat, grp = g)$info
    crossprod(Tm, Ls %*% Tm)
  })
}

# Generalized inverse and numerical rank of a symmetric information matrix,
# after scaling to unit diagonal so that very different parameter scales
# (e.g. a flat step in the IRT parameterization) do not hide null directions.
.psd_inv <- function(A, tol = sqrt(.Machine$double.eps)) {
  d <- sqrt(pmax(diag(A), 0)); ok <- d > tol * max(stats::median(d), .Machine$double.xmin)
  inv <- matrix(0, nrow(A), ncol(A))
  if (!any(ok)) return(list(inv = inv, rank = 0L))
  B <- A[ok, ok, drop = FALSE] / tcrossprod(d[ok])
  ev <- eigen((B + t(B)) / 2, symmetric = TRUE)
  keep <- ev$values > tol * max(ev$values[1], 0)
  Bi <- ev$vectors[, keep, drop = FALSE] %*% (t(ev$vectors[, keep, drop = FALSE]) / ev$values[keep])
  inv[ok, ok] <- Bi / tcrossprod(d[ok])
  list(inv = inv, rank = sum(keep))
}

# Efficient scores of the columns t given all other columns (OPG information)
.efficient <- function(S, t) {
  St <- S[, t, drop = FALSE]
  if (length(t) == ncol(S)) return(St)
  .resid_span(St, S[, -t, drop = FALSE])
}

# Residuals of Y after projection onto the numerically identified column space
# of X. Equal to Y - X (X'X)^{-1} X'Y when X has full column rank; when it does
# not (e.g. the score of a step difficulty whose discrimination is 0 is a zero
# column), the dependent directions are dropped instead of failing in solve().
.resid_span <- function(Y, X, tol = sqrt(.Machine$double.eps)) {
  X <- .unit_columns(X, tol)
  if (ncol(X) == 0) return(Y)
  s <- svd(X, nv = 0)
  U <- s$u[, s$d > tol * s$d[1], drop = FALSE]
  Y - U %*% crossprod(U, Y)
}

# Drop (numerically) zero columns, judged against the median column norm so that
# one huge column (e.g. the slope score of a flat step in the IRT
# parameterization) does not make all others look negligible, and scale the rest
# to unit length.
.unit_columns <- function(X, tol, ref = NULL) {
  cn <- sqrt(colSums(X^2))
  if (is.null(ref)) ref <- stats::median(cn)
  keep <- cn > tol * max(ref, .Machine$double.xmin)
  sweep(X[, keep, drop = FALSE], 2, cn[keep], "/")
}

# Rotate the tested scores onto their numerically identified directions. The
# LM-type statistics are invariant to nonsingular transformations of the score
# columns, so dropping null directions only reduces the degrees of freedom.
.full_rank_scores <- function(St, ref = NULL, tol = sqrt(.Machine$double.eps)) {
  Z <- .unit_columns(St, tol, ref)
  if (ncol(Z) == 0) return(Z)
  s <- svd(Z)
  r <- s$d > tol * s$d[1]
  if (all(r) && ncol(Z) == ncol(St)) return(St)
  Z %*% s$v[, r, drop = FALSE]
}

# Block selected by 'parm': "slope", "location", "all" (both), or NULL (column names)
.parm_block <- function(parm) {
  key <- c(disc = "slope", slope = "slope", diff = "location", int = "location", all = "all")
  if (!is.character(parm) || !all(parm %in% names(key))) return(NULL)
  b <- unique(key[parm])
  if ("all" %in% b || all(c("slope", "location") %in% b)) "all" else b
}

.dif_functional_name <- function(z, functional) {
  if (!is.null(functional)) return(match.arg(functional, c("LM", "maxLMo", "DM", "CvM", "maxLM")))
  if (is.ordered(z)) "maxLMo" else if (is.factor(z)) "LM" else "CvM"
}

.dif_functional <- function(name, gp) {
  switch(name,
         LM     = strucchange::catL2BB(gp),
         maxLMo = strucchange::ordwmax(gp),
         DM     = strucchange::maxBB,
         CvM    = strucchange::meanL2BB,
         maxLM  = strucchange::supLM(0.1))
}

