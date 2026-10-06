#' Item and group parameters in a compact table
#'
#' A clean report of a fitted TPPCM-family model, in the spirit of
#' `coef(mod, IRTpars = TRUE, simplify = TRUE)` in mirt: one row per item,
#' one column per parameter, standard errors in a matching table, and the
#' latent mean and variance.
#'
#' With `IRTpars = TRUE` the step kernel is `disc * (theta - diff)`: `disc`
#' is the step discrimination and `diff` the step difficulty (Masters' step
#' difficulty; TAM's `beta + tau`; mirt's `b` with `IRTpars = TRUE`;
#' `psychotools::threshpar()`). With `IRTpars = FALSE` it is
#' `slope * theta + int`, where the step intercepts `int` cumulate to the
#' category intercepts of mirt / psychotools (`d_k`), which equal `-AXsi_` in TAM.
#' In the terms of the category boundary literature `disc` (= `slope`) is
#' the CBD, `diff` the CBL and `int` the IBD (see [tppcm-package]).
#' Standard errors come from the observed information of the fitted model
#' (delta method for `diff`), not from TAM's `se.gammaslope`/`se.xsi`.
#'
#' @param x A fitted model (see [get_parts()]) or a `get_parts` object.
#' @param IRTpars Report discriminations and difficulties (default) or slopes
#'   and intercepts.
#' @param se Include standard errors.
#' @param fit Include model fit: log-likelihood, number of parameters, AIC,
#'   BIC and the limited-information statistics of [m2()] (M2, RMSEA2 with
#'   confidence interval, SRMSR, CFI, TLI). M2 needs a single group;
#'   otherwise only the information criteria are shown.
#' @param long Return a long table (one row per parameter) with standard
#'   errors and confidence limits instead of the compact report.
#' @param level Confidence level for `long = TRUE` or `restricted = TRUE`.
#' @param restricted Return the parameters of the fitted model instead, e.g.
#'   the item slopes of the GPCM or \eqn{\alpha_i, \gamma_l} of the product
#'   form, with standard errors and confidence limits.
#' @details
#' A step discrimination that is numerically 0 (an estimate on the bound of
#' `xxirt_tppcm()`) or practically flat (`|disc| < 0.1`) makes the step
#' difficulty `diff = -int / disc` meaningless; for such steps `diff` and its
#' SE are left blank and the step is listed below the table. The
#' slope-intercept form (`IRTpars = FALSE`) is well behaved for them.
#'
#' With `fit = TRUE` the header includes [m2()], which takes a few seconds for
#' long tests (e.g. 24 items); use `fit = FALSE` to skip it.
#'
#' @return Object of class `irt_pars`: a list with `items` (estimates, items x
#'   parameters), `se` (standard errors or `NULL`), `group` (latent mean and
#'   variance per group), `model`, `label`, `IRTpars`, `fit` (the header) and
#'   `at_bound`, `near_zero` (names of step discriminations at or near 0).
#'   `as.data.frame()` gives the long table. With `long = TRUE` or
#'   `restricted = TRUE`, a data frame (class `tppcm_table`).
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(800, disc = outer(c(0.8, 1, 1.2, 1.5), c(1, 1.5, 0.7)),
#'                  diff = matrix(c(-1, 0, 1), 4, 3, byrow = TRUE))
#' irt_pars(tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE))
#' fit <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
#' irt_pars(fit)
#' irt_pars(fit, long = TRUE)       # long format with confidence limits
#' fit$se.gammaslope                # TAM's standard errors: too small
#' }
#' @export
irt_pars <- function(x, IRTpars = TRUE, se = TRUE, fit = TRUE, long = FALSE,
                     level = 0.95, restricted = FALSE) {
  x <- get_parts(x)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("'level' must be a number between 0 and 1, e.g. 0.95", call. = FALSE)
  z <- stats::qnorm(1 - (1 - level) / 2)
  if (restricted) {
    r <- x$restricted
    out <- data.frame(est = r$est, se = r$se, lower = r$est - z * r$se,
                      upper = r$est + z * r$se, row.names = r$par)
    return(.as_table(out, sprintf("Parameters of the fitted model (%s), %g%% intervals",
                                  .label(x), 100 * level)))
  }
  pidx <- x$parindex; L <- ncol(x$a)
  ia <- which(pidx$type == "a"); id <- which(pidx$type == "d")
  a <- x$par[ia]; d <- x$par[id]
  V <- x$vcov_sat
  va <- diag(V)[ia]; vd <- diag(V)[id]
  pcm <- x$model == "pcm"
  at_bound <- !pcm & abs(a) < 1e-3                               # numerically 0 (bound of xxirt)
  near0 <- !pcm & !at_bound & abs(a) < 0.1                       # practically flat step
  if (IRTpars && pcm) {
    # TAM's PCM: slope fixed at 1, latent variance free. With the trait scale held
    # at its estimate the variance acts as a common slope c (VAR = c^2 VAR_hat),
    # and the location is b = d with SE(b) = SE(d).
    loc <- d; vloc <- vd
  } else if (IRTpars) {
    loc <- d / a
    cov_ad <- V[cbind(ia, id)]
    vloc <- vd / a^2 + d^2 * va / a^4 - 2 * d * cov_ad / a^3     # delta method
    loc[at_bound | near0] <- NA; vloc[at_bound | near0] <- NA    # difficulty undefined
  } else {
    loc <- -d; vloc <- vd                                        # int = -d
  }
  lab_s <- if (IRTpars) "disc" else "slope"; lab <- if (IRTpars) "diff" else "int"
  wide <- function(va_, vl_) {
    M <- matrix(NA_real_, length(x$items), 2 * L,
                dimnames = list(x$items, c(paste0(lab_s, 1:L), paste0(lab, 1:L))))
    M[cbind(pidx$itemnr[ia], pidx$step[ia])] <- va_
    M[cbind(pidx$itemnr[id], L + pidx$step[id])] <- vl_
    M
  }
  a[at_bound] <- 0; va[at_bound] <- NA
  items <- wide(a, loc)
  ses <- if (se) wide(sqrt(pmax(va, 0)), sqrt(pmax(vloc, 0))) else NULL
  if (se) {
    ses[abs(ses) < 1e-12] <- NA                    # fixed parameters
    if (pcm) ses[, seq_len(L)] <- NA               # slopes fixed at 1 in TAM's PCM
  }
  group <- x$groups
  if (!is.null(group)) {
    if (is.null(group$SE_VAR)) group$SE_VAR <- NA_real_
    if (pcm && nrow(group) == 1)
      group$SE_VAR <- 2 * group$VAR * x$restricted$se[x$restricted$par == "slope"]
  }
  if (pcm) items[, seq_len(L)][!is.na(items[, seq_len(L)])] <- 1
  if (long) {
    est <- c(if (pcm) rep(1, length(ia)) else a, loc)
    sev <- c(if (pcm) rep(NA, length(ia)) else sqrt(pmax(va, 0)), sqrt(pmax(vloc, 0)))
    sev[!is.na(sev) & sev < 1e-12] <- NA
    type <- rep(c(lab_s, lab), c(length(ia), length(id)))
    item <- c(pidx$item[ia], pidx$item[id]); step <- c(pidx$step[ia], pidx$step[id])
    o <- order(match(item, x$items), type != lab_s, step)
    out <- data.frame(item = item, step = step, par = type, est = est, se = sev,
                      lower = est - z * sev, upper = est + z * sev)[o, ]
    rownames(out) <- paste0(out$item, "_", out$par, out$step)
    return(.as_table(out, sprintf("%s, %g%% intervals", .label(x), 100 * level)))
  }
  fit_tab <- if (fit) .fit_summary(x) else NULL
  structure(list(items = items, se = ses, group = group, model = x$model,
                 label = .label(x), IRTpars = IRTpars, fit = fit_tab,
                 at_bound = .disc_name(pidx$name[ia][at_bound]),
                 near_zero = .disc_name(pidx$name[ia][near0]), parts = x, level = level),
            class = "irt_pars")
}

#' @export
as.data.frame.irt_pars <- function(x, ...) {
  as.data.frame(unclass(irt_pars(x$parts, IRTpars = x$IRTpars, long = TRUE, level = x$level)))
}

#' @export
print.irt_pars <- function(x, digits = 3, ...) {
  fmt <- function(M) {
    out <- formatC(round(M, digits), format = "f", digits = digits)
    out[is.na(M)] <- ""
    dimnames(out) <- dimnames(M)
    out
  }
  cat(sprintf("Model: %s | %s\n", x$label,
              if (x$IRTpars) "kernel disc * (theta - diff)" else "kernel slope * theta + int"))
  if (!is.null(x$fit)) .print_fit(x$fit, digits)
  cat("\n")
  cat("$items\n"); print(noquote(fmt(x$items)), right = TRUE)
  if (!is.null(x$se)) { cat("\n$se\n"); print(noquote(fmt(x$se)), right = TRUE) }
  if (length(x$at_bound)) {
    cat("\nOn the bound 0 (no SE; tests involving them are not valid):",
        paste(x$at_bound, collapse = ", "), "\n")
  }
  if (length(x$near_zero) && x$IRTpars) {
    cat("\nNear 0 (|disc| < 0.1), difficulty not reported; see IRTpars = FALSE:",
        paste(x$near_zero, collapse = ", "), "\n")
  }
  if (!is.null(x$group)) {
    G <- x$group
    if (all(is.na(G$SE_VAR))) G$SE_VAR <- NULL
    if (!is.null(G$SE_MEAN) && all(is.na(G$SE_MEAN))) G$SE_MEAN <- NULL
    cat("\n$group\n"); print(G, digits = digits, row.names = FALSE)
  }
  invisible(x)
}

# Information criteria and M2-based fit for the header of irt_pars()
.fit_summary <- function(x) {
  N <- nrow(x$data)
  npar <- .npar(x)
  ll <- x$loglik
  out <- list(N = N, logLik = ll, npar = npar, AIC = -2 * ll + 2 * npar,
              BIC = -2 * ll + log(N) * npar)
  out$m2 <- tryCatch(m2(x), error = function(e) conditionMessage(e))
  out
}

.print_fit <- function(f, digits = 3) {
  cat(sprintf("Fit: N = %d, logLik = %.2f, npar = %d, AIC = %.2f, BIC = %.2f\n",
              f$N, f$logLik, f$npar, f$AIC, f$BIC))
  if (is.character(f$m2)) { cat("     (no M2: ", f$m2, ")\n", sep = ""); return(invisible()) }
  y <- as.data.frame(unclass(f$m2))
  fm <- function(v) sub("^0\\.", ".", formatC(v, format = "f", digits = digits))
  cat(sprintf("     %s(%d) = %.2f, p = %s, RMSEA = %s [%s, %s], SRMSR = %s, CFI = %s, TLI = %s\n",
              y$type, y$df, y$M2, fm(y$p), fm(y$RMSEA), fm(y$RMSEA_lower), fm(y$RMSEA_upper),
              fm(y$SRMSR), fm(y$CFI), fm(y$TLI)))
}

# Data frame with a title and rounded, readable printing
.as_table <- function(df, title, empty = NULL) {
  if (is.null(df) || !nrow(df)) {
    if (!is.null(empty)) message(empty)
    df <- data.frame()
  }
  structure(df, class = c("tppcm_table", "data.frame"), title = title)
}

#' @export
print.tppcm_table <- function(x, digits = 3, ...) {
  if (!is.null(attr(x, "title"))) cat(attr(x, "title"), "\n")
  if (!nrow(x)) { cat("(nothing to test)\n"); return(invisible(x)) }
  y <- as.data.frame(unclass(x), stringsAsFactors = FALSE, optional = TRUE)
  rownames(y) <- rownames(x)
  for (nm in names(y)) {
    v <- y[[nm]]
    if (!is.numeric(v) || all(v == round(v), na.rm = TRUE)) next
    y[[nm]] <- if (grepl("^p($|_)", nm)) format.pval(v, digits = digits, eps = 1e-4)
               else formatC(v, format = "f", digits = digits)
    y[[nm]][is.na(v)] <- ""
  }
  print(y, right = TRUE, row.names = !identical(rownames(y), as.character(seq_len(nrow(y)))))
  invisible(x)
}
