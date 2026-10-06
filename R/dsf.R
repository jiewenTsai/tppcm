#' Likelihood ratio tests of differential step functioning
#'
#' Partial invariance search at the level of category boundaries, in the way
#' of `group.equal` and `group.partial` in lavaan. The baseline model is the
#' multigroup TPPCM with the step parameters in `group.equal` equal across
#' groups (except the steps in `group.partial`); each restricted parameter is
#' then released one at a time and compared with the baseline by a likelihood
#' ratio test (df = number of groups - 1). The latent means and variances of
#' groups 2, ..., G are estimated as far as the baseline identifies them, and
#' all other steps serve as anchors. Use it to decompose the steps flagged by
#' the joint screen `dif_test(fit, g, by = "step")`:
#' `dsf(fit, items = c("I2:2", "I5:1"))`.
#'
#' * `group.equal = c("cbd", "ibd")` (default, scalar baseline): every step
#'   gives two tests. The IBD (location) test releases the IBD of the step.
#'   The CBD test releases the CBD of the step with its location free, i.e. it
#'   compares the model with CBD and IBD of the step released against the
#'   model with only the IBD released, so that it does not depend on the
#'   origin of the scale and location DSF does not leak into it. The two LR
#'   statistics of a step add up to the joint 2-df LR test of the step.
#' * `group.equal = "cbd"` (metric baseline, all IBDs free): one CBD test per
#'   step, as in testing loadings at the metric stage.
#'
#' `"cbl"` may be written for `"ibd"`: the location rows then report the CBL
#' (CBL = -IBD / CBD) instead of the IBD; with equal CBDs the two restrictions
#' are the same. As in invariance testing, interpret the location of a step
#' only when its CBD is invariant. A significant test can also reflect DSF in
#' other steps, which are assumed invariant (anchor contamination); with many
#' tests use `p_adj`. See [tppcm-package] for CBD, CBL and IBD.
#'
#' Only the data, the groups, the within-group step design (`design` or
#' `index` of [tppcm()]) and the TAM control settings are taken from `fit`;
#' the baseline is refitted with [dsf_design()]. Every test refits the model
#' once or twice, so restrict `items` to the steps of interest for long tests.
#'
#' @param fit A multigroup [TAM::tam.mml.3pl()] fit with a step design `E`
#'   from [tppcm()] and `group = `, e.g.
#'   `tam.mml.3pl(dat, E = tppcm(dat), group = g, est.variance = FALSE)`.
#' @param group.equal `c("cbd", "ibd")` (default) or `"cbd"`; `"cbl"` may be
#'   used for `"ibd"`.
#' @param group.partial Steps already released in the baseline, as in
#'   [dsf_design()]; they are not tested.
#' @param items Items (names) or single steps (`"item:step"`, as in
#'   `group.partial` and the units of `dif_test(by = "step")`) to test;
#'   default all.
#' @param steps Steps to test within `items` (numbers); default all.
#' @param adjust Multiple-testing adjustment of `p_adj` over all rows, a method
#'   of [stats::p.adjust()] (`"none"` for none).
#' @param control TAM `control` list for the refits; default the fit's
#'   control with `conv <= 1e-6`, `convD <= 1e-5` and `maxiter >= 3000`.
#' @param verbose Print progress.
#' @return Data frame (class `tppcm_table`) with one row per test: `item`,
#'   `step`, `par` (`"cbd"`, `"ibd"` or `"cbl"`), the estimate in every group
#'   from the model with the parameter released, `LR`, `df`, `p`, `p_adj`.
#'   Attributes `"impact"` (latent means and variances of the baseline) and
#'   `"loglik0"`.
#' @seealso [dif_test()] for the score-based screen without refitting;
#'   [dsf_design()] for configural, metric, scalar and partial models.
#' @references Stark, S., Chernyshenko, O. S., & Drasgow, F. (2006). Detecting
#'   differential item functioning with confirmatory factor analysis and item
#'   response theory: Toward a unified strategy. *Journal of Applied
#'   Psychology, 91*, 1292-1306.
#' @examples
#' \donttest{
#' set.seed(1)
#' a <- matrix(c(1, 1.5, 1.2), 4, 3, byrow = TRUE)
#' b <- matrix(c(-1, 0, 1), 4, 3, byrow = TRUE)
#' a2 <- a; a2[2, 2] <- 2.5
#' dat <- rbind(sim_tppcm(500, disc = a, diff = b),
#'              sim_tppcm(500, disc = a2, diff = b, theta = rnorm(500, 0.3, 0.8)))
#' g <- rep(c("ref", "foc"), each = 500)
#' fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), group = g, est.variance = FALSE,
#'                         verbose = FALSE)
#' dif_test(fit, g, by = "step")                 # screen: joint 2-df test per step
#' dsf(fit, items = c("I2:2", "I3:2"), verbose = FALSE)          # decompose single steps
#' dsf(fit, group.equal = "cbd", items = "I2", verbose = FALSE)  # metric baseline
#' }
#' @export
dsf <- function(fit, group.equal = c("cbd", "ibd"), group.partial = NULL, items = NULL,
                steps = NULL, adjust = "holm", control = NULL, verbose = TRUE) {
  adjust <- .check_adjust(adjust)
  ge <- .dsf_types(group.equal, "group.equal")
  if (!ge[["cbd"]]) stop('group.equal must contain "cbd": c("cbd", "ibd") or "cbd"', call. = FALSE)
  loc <- if ("cbl" %in% group.equal) "cbl" else "ibd"
  if (!inherits(fit, "tam.mml.3pl"))
    stop("fit must be a TAM::tam.mml.3pl() fit with E = tppcm(dat) and group = g, not an object of class '",
         class(fit)[1], "'", call. = FALSE)
  if (is.null(fit$group) || length(unique(fit$group)) < 2)
    stop("fit has a single group; refit with tam.mml.3pl(..., group = g, est.variance = FALSE)", call. = FALSE)
  if (!is.null(fit$Y) && NCOL(fit$Y) > length(unique(fit$group)))
    stop("latent regression fits are not supported", call. = FALSE)
  dat <- as.matrix(fit$resp); itm <- .item_names(dat)
  group <- if (!is.null(fit$groups)) fit$groups[fit$group] else fit$group
  G <- length(unique(group)); gl <- as.character(sort(unique(group)))
  K <- apply(dat, 2, max, na.rm = TRUE) + 1
  index <- .index_from_E(fit$E, K)
  if (is.null(control)) {
    control <- fit$control
    control$conv <- min(control$conv, 1e-6); control$convD <- min(control$convD, 1e-5)
    control$maxiter <- max(control$maxiter, 3000)
  }
  control$progress <- FALSE
  has <- !is.na(index)
  part0 <- .dsf_partial(group.partial, itm, has, ge)
  if (loc == "cbl") part0$ibd <- part0$ibd | (part0$cbd & ge[["ibd"]])
  design <- function(cbd = NULL, ibd = NULL) {
    p <- part0
    if (!is.null(cbd)) p$cbd[cbd[1], cbd[2]] <- TRUE
    if (!is.null(ibd)) p$ibd[ibd[1], ibd[2]] <- TRUE
    dsf_design(dat, group, group.equal = c("cbd", "ibd")[ge], group.partial = p,
               index = index, K = K)
  }
  refit <- function(d) suppressWarnings(TAM::tam.mml.3pl(
    d$resp, A = d$A, E = d$E, group = group, variance.fixed = d$variance.fixed,
    beta.fixed = d$beta.fixed, pweights = fit$pweights, control = control, verbose = FALSE))
  ll <- function(m) -m$ic$deviance / 2
  m0 <- refit(design())
  sel <- if (is.null(items)) has else .dsf_steps(items, itm, has, "items")
  if (!is.null(steps)) sel[, setdiff(seq_len(ncol(sel)), steps)] <- FALSE
  units <- which(sel, arr.ind = TRUE)
  units <- data.frame(i = units[, 1], l = units[, 2])
  units <- units[order(units$i, units$l), , drop = FALSE]
  if (!nrow(units)) stop("no steps to test; check 'items' and 'steps'", call. = FALSE)
  step_pars <- function(m, d, i, l) {
    vn <- paste0(itm[i], "_", gl)
    a <- unname(m$gammaslope[d$index[vn, l]])
    int <- -unname(m$xsi$xsi[d$index_int[vn, l]])
    list(cbd = a, ibd = int, cbl = -int / a)
  }
  row <- function(i, l, par, est, lr) {
    r <- data.frame(item = itm[i], step = l, par = par)
    for (g in seq_len(G)) r[[gl[g]]] <- est[g]
    r$LR <- max(0, lr); r$df <- G - 1L
    r
  }
  rows <- list()
  for (u in seq_len(nrow(units))) {
    i <- units$i[u]; l <- units$l[u]
    if (verbose) message(sprintf("[%d/%d] %s:%d", u, nrow(units), itm[i], l))
    ibd_eq <- ge[["ibd"]] && !part0$ibd[i, l]
    m_loc <- NULL
    if (ibd_eq) {                                  # location test (and location-free model for the CBD test)
      d_loc <- design(ibd = c(i, l)); m_loc <- refit(d_loc)
      p <- step_pars(m_loc, d_loc, i, l)
      rows[[length(rows) + 1]] <- row(i, l, loc, p[[loc]], 2 * (ll(m_loc) - ll(m0)))
    }
    if (!part0$cbd[i, l]) {                        # CBD test with the location of the step free
      d_cbd <- design(cbd = c(i, l), ibd = c(i, l)); m_cbd <- refit(d_cbd)
      ll_null <- if (ibd_eq) ll(m_loc) else ll(m0)
      p <- step_pars(m_cbd, d_cbd, i, l)
      rows[[length(rows) + 1]] <- row(i, l, "cbd", p$cbd, 2 * (ll(m_cbd) - ll_null))
    }
  }
  if (!length(rows))
    stop("no restricted steps to test: the selected steps are all released in group.partial", call. = FALSE)
  out <- do.call(rbind, rows)
  out <- out[order(match(out$item, itm), out$step, out$par != "cbd"), ]
  rownames(out) <- NULL
  out$p <- stats::pchisq(out$LR, out$df, lower.tail = FALSE)
  out$p_adj <- stats::p.adjust(out$p, adjust)
  grp <- data.frame(group = gl, MEAN = .grp_means(m0), VAR = as.numeric(m0$variance[, 1, 1]))
  base <- if (ge[["ibd"]]) "scalar baseline (CBD and IBD equal)" else "metric baseline (CBD equal)"
  part <- if (any(part0$cbd | part0$ibd)) ", partial" else ""
  title <- sprintf("LR tests of DSF, one step parameter released at a time vs the %s%s | reference group: %s | %s",
                   base, part, gl[1], .adj_label(adjust))
  structure(out, class = c("tppcm_table", "data.frame"), title = title, impact = grp,
            loglik0 = ll(m0))
}

# group means of a multigroup TAM fit (reference group 0)
.grp_means <- function(m) {
  G <- length(unique(m$group))
  b <- as.numeric(m$beta[, 1])
  if (length(b) == G) b - b[1] else c(0, b[-1])
}

# items x steps parameter index of a tppcm() step design array E
.index_from_E <- function(E, K) {
  I <- dim(E)[1]; L <- max(K) - 1
  idx <- matrix(NA_integer_, I, L)
  for (i in seq_len(I)) for (l in seq_len(K[i] - 1)) {
    inc <- E[i, l + 1, 1, ] - E[i, l, 1, ]
    p <- which(abs(inc) > 1e-12)
    if (length(p) != 1 || abs(inc[p] - 1) > 1e-12)
      stop("fit$E is not a step design from tppcm(); refit with E = tppcm(dat, ...)", call. = FALSE)
    idx[i, l] <- p
  }
  idx
}
