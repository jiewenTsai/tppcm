# User-facing names of saturated parameters: I1_a2 -> I1_disc2, I1_d2 -> I1_d2
.disc_name <- function(nm) sub("_a([0-9]+)$", "_disc\\1", nm)

# Score statistic for extra directions D (columns in the saturated space)
.score_stat <- function(x, D, info) {
  I  <- if (info == "louis") x$info else x$info_opg
  J0 <- x$J; g <- x$grad
  D  <- .pad_rows(D, nrow(J0))
  IJJ <- crossprod(J0, I %*% J0) - x$curv
  IDJ <- crossprod(D, I %*% J0)
  M   <- crossprod(D, I %*% D) - IDJ %*% solve(IJJ, t(IDJ))
  gD  <- crossprod(D, g) - IDJ %*% solve(IJJ, crossprod(J0, g))
  Minv <- MASS::ginv(M)
  stat <- as.numeric(crossprod(gD, Minv %*% gD))
  df <- qr(cbind(J0, D))$rank - qr(J0)$rank
  list(stat = stat, df = df, p = stats::pchisq(stat, df, lower.tail = FALSE),
       epc = stats::setNames(as.numeric(Minv %*% gD), colnames(D)))
}

#' Score test against a less restricted TPPCM
#'
#' Tests the restricted model of `x` against `against` without fitting the
#' alternative. With extra directions \eqn{D},
#' \deqn{S = g_D' \big(D'ID - D'IJ (J'IJ)^{-1} J'ID\big)^{-1} g_D,}
#' where \eqn{g_D} is the score in the directions \eqn{D} adjusted for the
#' restricted parameters and \eqn{J'IJ} includes the curvature term.
#' Against the saturated model \eqn{S = g'I^{-1}g}.
#'
#' @param x A fitted model (see [get_parts()]) or a `get_parts` object.
#' @param against Alternative model as a design formula (see [tppcm()]); must
#'   contain the model of `x`: `~ item * step` (default, saturated),
#'   `~ item + step` (Rank-1, \eqn{\alpha_i\gamma_l}), `~ item` (GPCM,
#'   \eqn{\alpha_i}) or `~ step` (\eqn{\gamma_l}, shared by items). The
#'   models form a hierarchy: `~ 1` is contained in `~ item` and `~ step`,
#'   both are contained in `~ item + step`, and that in `~ item * step`. Thus
#'   the GPCM against `~ item + step` tests for a step effect, Rank-1 against
#'   `~ item * step` for an item x step interaction.
#' @param D Optional matrix of extra directions in the saturated space
#'   (rows as in `x$parindex`); overrides `against`.
#' @param info `"louis"` (observed information) or `"opg"`.
#' @return Object of class `tppcm_test` with `statistic`, `df`, `p.value`, `epc`.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(800, disc = outer(c(0.8, 1, 1.2, 1.5), c(1, 1.5, 0.7)),
#'                  diff = matrix(c(-1, 0, 1), 4, 3, byrow = TRUE))
#' fit <- tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
#' score_test(fit)                   # GPCM vs saturated TPPCM
#' score_test(fit, against = ~ item + step)   # GPCM vs Rank-1: step effect?
#' item_test(fit)
#' mi(fit)
#' }
#' @export
score_test <- function(x, against = ~ item * step, D = NULL, info = c("louis", "opg")) {
  info <- match.arg(info)
  x <- get_parts(x)
  if (is.null(D)) {
    against <- .parse_design(against, "against")
    if (against == "rank1" && length(unique(x$ncat)) > 1)
      stop("Rank-1 (against = ~ item + step) needs the same number of categories in every item",
           call. = FALSE)
    Jalt <- .pad_rows(.jacobian(against, x$parindex, x$par), nrow(x$J))
    rownames(Jalt) <- rownames(x$J)
    if (qr(cbind(Jalt, x$J[, setdiff(colnames(x$J), x$gnames), drop = FALSE]))$rank > qr(Jalt)$rank)
      stop("the model ", .design_codes[[against]], " does not contain the fitted model (", .label(x),
           "); 'against' must be a larger model (~ 1 < ~ item, ~ step < ~ item + step < ~ item * step)",
           call. = FALSE)
    D <- tryCatch(.extra_dirs(x$J, Jalt), error = function(e)
      stop("the fitted model (", .label(x), ") already contains the model ", .design_codes[[against]],
           "; nothing to test", call. = FALSE))
  } else against <- "user-defined directions"
  st <- .score_stat(x, D, info)
  structure(list(statistic = st$stat, df = st$df, p.value = st$p, epc = st$epc,
                 model = x$model, label = .label(x), against = against, info = info),
            class = "tppcm_test")
}

# "the saturated TPPCM (~ item * step)", ...
.against_label <- stats::setNames(paste0("the ", .design_names, " (", .design_codes, ")"), names(.design_codes))

#' @export
print.tppcm_test <- function(x, digits = 4, ...) {
  against <- if (x$against %in% names(.against_label)) .against_label[[x$against]] else x$against
  cat(sprintf("Score test (%s information) of the fitted %s against %s\n",
              if (identical(x$info, "opg")) "OPG" else "Louis",
              if (is.null(x$label)) x$model else x$label, against))
  cat(sprintf("  S = %.*f, df = %d, p = %s\n", digits, x$statistic, x$df,
              format.pval(x$p.value, digits = digits)))
  invisible(x)
}

#' Modification indices for single step discriminations
#'
#' For each step discrimination that is restricted in `x`, the 1-df score
#' test for freeing it (MI) and its expected parameter change (EPC). The EPC
#' is the first-order (one Newton step) approximation to the change of the
#' estimate when the step is freed; when the MI is large it can overshoot the
#' actual change noticeably.
#'
#' For an item with two steps, freeing one step is the same 1-df test as
#' freeing the other (the two rows have equal MI and EPCs of opposite sign).
#'
#' @inheritParams score_test
#' @param pars Names of step discriminations (default: all, e.g. `"I1_disc2"`).
#' @param adjust Multiple-testing adjustment of `p_adj` over all rows, a method
#'   of [stats::p.adjust()] (`"none"` for none).
#' @return Data frame (class `tppcm_table`) with `par`, `MI`, `df`, `p`,
#'   `p_adj`, `EPC` (expected change of the step discrimination), sorted by
#'   `MI`.
#' @examples
#' \donttest{
#' data(data.Students, package = "CDM")
#' sc <- data.Students[, c("sc1", "sc2", "sc3", "sc4")]
#' fit <- tam.mml.2pl(sc, irtmodel = "GPCM", verbose = FALSE)
#' mi(fit)                                  # which steps deviate from the GPCM
#' mi(fit, pars = c("sc1_disc3", "sc2_disc3"))
#' }
#' @export
mi <- function(x, pars = NULL, info = c("louis", "opg"), adjust = "holm") {
  info <- match.arg(info); adjust <- .check_adjust(adjust)
  x <- get_parts(x)
  pidx <- x$parindex
  disc_names <- .disc_name(pidx$name)
  if (is.null(pars)) pars <- disc_names[pidx$type == "a"]
  idx <- match(pars, disc_names)
  idx[is.na(idx)] <- match(pars[is.na(idx)], pidx$name)            # internal names also accepted
  if (anyNA(idx)) stop("unknown step discrimination(s): ", paste(pars[is.na(idx)], collapse = ", "),
                       "; use names such as \"", disc_names[pidx$type == "a"][1], "\"", call. = FALSE)
  Q <- qr.Q(qr(x$J))
  out <- do.call(rbind, lapply(idx, function(j) {
    nm <- disc_names[j]
    D <- .pad_rows(diag(nrow(pidx))[, j, drop = FALSE], nrow(x$J))
    colnames(D) <- nm
    if (sqrt(sum((D - Q %*% crossprod(Q, D))^2)) < 1e-7) return(NULL)   # already free
    st <- .score_stat(x, D, info)
    data.frame(par = nm, MI = st$stat, df = st$df, p = st$p, EPC = st$epc)
  }))
  if (is.null(out)) return(.as_table(NULL, "Modification indices for step discriminations",
                                      empty = "all step discriminations are already free in the fitted model"))
  out <- .add_p_adj(out, adjust)[, c("par", "MI", "df", "p", "p_adj", "EPC")]
  out <- out[order(-out$MI), ]
  rownames(out) <- NULL
  .as_table(out, sprintf("Modification indices (1-df score tests) for freeing single step discriminations of the fitted %s | %s",
                         .label(x), .adj_label(adjust)))
}

#' Item-level score tests
#'
#' For each item, the score test for freeing all of its step discriminations
#' relative to the model of `x` (e.g. whether the item deviates from the GPCM).
#'
#' @inheritParams mi
#' @return Data frame (class `tppcm_table`) with `item`, `statistic`, `df`,
#'   `p`, `p_adj`.
#' @examples
#' \donttest{
#' data(data.Students, package = "CDM")
#' sc <- data.Students[, c("sc1", "sc2", "sc3", "sc4")]
#' fit <- tam.mml.2pl(sc, irtmodel = "GPCM", verbose = FALSE)
#' item_test(fit)                           # which items deviate from the GPCM
#' }
#' @export
item_test <- function(x, info = c("louis", "opg"), adjust = "holm") {
  info <- match.arg(info); adjust <- .check_adjust(adjust)
  x <- get_parts(x)
  pidx <- x$parindex
  out <- do.call(rbind, lapply(unique(pidx$item), function(it) {
    sel <- which(pidx$item == it & pidx$type == "a")
    U <- .pad_rows(diag(nrow(pidx))[, sel, drop = FALSE], nrow(x$J)); colnames(U) <- pidx$name[sel]
    D <- tryCatch(.extra_dirs(x$J, U), error = function(e) NULL)
    if (is.null(D)) return(NULL)
    st <- .score_stat(x, D, info)
    data.frame(item = it, statistic = st$stat, df = st$df, p = st$p)
  }))
  .as_table(.add_p_adj(out, adjust),
            sprintf("Score tests for freeing the step discriminations of each item in the fitted %s | %s",
                    .label(x), .adj_label(adjust)),
            empty = "all step discriminations are already free in the fitted model")
}

#' Wald test of equal step discriminations within items
#'
#' For each item, tests \eqn{a_{i1} = \dots = a_{iL}} (the item follows the
#' GPCM) with `CDM::WaldTest()`. Items whose steps are already constrained
#' equal in `x` are skipped.
#'
#' @param x A fitted model (see [get_parts()]) or a `get_parts` object in
#'   which the tested steps are free.
#' @inheritParams mi
#' @return Data frame (class `tppcm_table`) with `item`, `statistic`, `df`,
#'   `p`, `p_adj`.
#' @examples
#' \donttest{
#' data(data.Students, package = "CDM")
#' sc <- data.Students[, c("sc1", "sc2", "sc3", "sc4")]
#' fit <- tam.mml.3pl(sc, E = tppcm(sc), est.variance = FALSE, verbose = FALSE)
#' wald_test(fit)                           # a_i1 = a_i2 = a_i3 for each item?
#' }
#' @export
wald_test <- function(x, adjust = "holm") {
  adjust <- .check_adjust(adjust)
  x <- get_parts(x)
  pidx <- x$parindex
  out <- do.call(rbind, lapply(unique(pidx$item), function(it) {
    ia <- which(pidx$item == it & pidx$type == "a")
    if (length(ia) < 2) return(NULL)
    R <- matrix(0, length(ia) - 1, length(x$par))
    for (r in seq_len(nrow(R))) R[r, ia[r:(r + 1)]] <- c(-1, 1)
    if (qr(R %*% x$vcov_sat %*% t(R))$rank < nrow(R)) return(NULL)
    w <- CDM::WaldTest(delta = x$par, vcov = x$vcov_sat, R = R, nobs = x$nobs)
    data.frame(item = it, statistic = w$X2, df = w$df, p = w$p)
  }))
  .as_table(.add_p_adj(out, adjust),
            sprintf("Wald tests of equal step discriminations within items (H0: the item follows the GPCM) | %s",
                    .adj_label(adjust)),
            empty = "no item has free step discriminations; fit the saturated TPPCM, e.g. E = tppcm(dat)")
}
