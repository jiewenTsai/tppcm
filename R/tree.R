# Case-wise scores of the fitted slopes and the step locations, in the IRT form
# disc * (theta - diff) ("irt") or the slope-intercept form slope * theta + int
# ("si"). Columns: free slopes of the design, then one location per step.
# Attribute "block": "slope" / "location".
.step_scores <- function(x, param = c("irt", "si")) {
  param <- match.arg(param)
  if (isTRUE(x$n_loc > 0)) stop("case-wise scores need free step locations (e.g. not an RSM fit)", call. = FALSE)
  pidx <- x$parindex
  ia <- which(pidx$type == "a"); id <- which(pidx$type == "d")
  Sa <- x$scores[, ia, drop = FALSE]; Sd <- x$scores[, id, drop = FALSE]
  a <- x$par[ia]; d <- x$par[id]
  A <- x$J[ia, .slope_cols(x$J, pidx), drop = FALSE]   # slope design
  N <- nrow(Sa)
  sl <- if (param == "irt") "disc" else "slope"
  lc <- if (param == "irt") "diff" else "int"
  slope_names <- apply(A, 2, function(col) {
    rows <- which(col != 0)
    if (length(rows) == 1) return(paste0(pidx$item[ia][rows], "_", sl, pidx$step[ia][rows]))
    it <- unique(pidx$item[ia][rows])
    if (length(it) == 1) paste0(it, "_", sl) else sl                 # item-level / common
  })
  if (param == "irt") {                               # d = a * diff
    S_slope <- (Sa + Sd * rep(d / a, each = N)) %*% A
    S_loc <- Sd * rep(a, each = N)
  } else {                                            # d = -int
    S_slope <- Sa %*% A
    S_loc <- -Sd
  }
  S_grp <- x$scores[, x$gnames, drop = FALSE]        # free group means and variances, if any
  out <- cbind(S_slope, S_loc, S_grp)
  colnames(out) <- c(unname(slope_names), paste0(pidx$item[id], "_", lc, pidx$step[id]), x$gnames)
  attr(out, "block") <- c(rep("slope", ncol(S_slope)), rep("location", ncol(S_loc)), rep("group", ncol(S_grp)))
  out
}

# Parameter estimates matching the columns of step_scores()
.step_coef <- function(x, param) {
  pidx <- x$parindex
  ia <- which(pidx$type == "a"); id <- which(pidx$type == "d")
  loc <- if (param == "irt") x$par[id] / x$par[ia] else -x$par[id]
  grp <- if (length(x$gnames)) {
    free <- as.integer(sub("grp:(mean|var)", "", x$gnames))
    ifelse(grepl("mean", x$gnames), x$groups$MEAN[free], x$groups$VAR[free])
  } else numeric(0)
  c(x$restricted$est, loc, grp)
}

#' TPPCM trees: model-based recursive partitioning of step discriminations
#'
#' Fits a TPPCM-family model with [TAM::tam.mml.3pl()] in every node and
#' splits the sample along covariates where the item parameters are unstable
#' (Zeileis, Hothorn & Hornik, 2008), using [partykit::mob()]. With
#' `parm = "disc"` (or `"slope"`) only instability of the slopes drives the
#' splits; with `parm = "diff"` (or `"int"`) only instability of the step
#' difficulties (intercepts); with `"all"` both.
#'
#' Node models: `"tppcm"` (step discriminations), `"gpcm"` (item
#' discriminations; compare [psychotree::gpcmtree()]) and `"pcm"` (one common
#' slope; compare [psychotree::pctree()]). With the trait fixed at
#' \eqn{N(0,1)}, the common slope of the PCM node plays the role of the latent
#' standard deviation, so for `"pcm"` use `parm = "diff"` to look for DIF
#' rather than for differences in trait variance.
#'
#' **Impact is not DIF.** The latent trait is \eqn{N(0,1)} in every node, so a
#' covariate related to ability (groups that differ in their mean) is split on
#' even without DIF: the difficulties of all items shift together. `impact`
#' names a grouping variable whose latent means and variances are estimated in
#' every node (`tam.mml.3pl(group = )`, as `impact` in
#' [psychotools::gpcmodel()]); splits on that variable then reflect item
#' differences only. Conditional-maximum-likelihood trees such as
#' [psychotree::pctree()] do not have this problem for the PCM.
#'
#' Rows with missing item responses are kept (TAM uses the observed responses);
#' the partitioning covariates must be complete. The node model needs every
#' category of every item in every node; otherwise the fit stops with a
#' message. Increase `minsize` or use a less demanding node model.
#'
#' @param formula `resp ~ z1 + z2 + ...`, where `resp` is a numeric matrix of
#'   item responses stored as one column of `data` (`data$resp <- as.matrix(items)`)
#'   and `z` are partitioning covariates.
#' @param data Data frame.
#' @param impact `NULL`, or the name of a factor in `data` (or a factor of
#'   length `nrow(data)`) whose latent means and variances are estimated in
#'   every node.
#' @param model Node model: `"tppcm"`, `"gpcm"` or `"pcm"`.
#' @param parm Parameters used in the instability tests: `"all"`, the slope
#'   block (`"disc"` or `"slope"`) or the location block (`"diff"` or `"int"`).
#' @param param Parameterization of the scores: `"irt"` (`disc`, `diff`) or
#'   `"si"` (`slope`, `int`); see `get_parts(fit, "estfun")`.
#' @param orthogonal For a slope (location) block test: replace the tested
#'   scores by efficient scores, i.e. residuals from their regression on the
#'   other block (\eqn{s_a - s_b \hat I_{bb}^{-1} \hat I_{ba}} with the
#'   outer-product information). Without this, `mob()` selects components of
#'   the jointly decorrelated process and instability in one block leaks into
#'   the test of the other (e.g. location DIF triggers slope splits).
#' @param ... Passed to [partykit::mob_control()] (e.g. `alpha`, `minsize`,
#'   `maxdepth`, `bonferroni`).
#' @return A `modelparty` object (class `tppcmtree`). `coef()` gives the
#'   node-wise parameters; `print()` shows the splits and the number of
#'   persons per node.
#' @examples
#' \donttest{
#' data("VerbalAggression", package = "psychotools")
#' tr <- tppcmtree(resp ~ gender + anger, data = VerbalAggression,
#'                 model = "pcm", parm = "diff", minsize = 50)
#' tr
#' head(t(coef(tr)))
#' }
#' @export
tppcmtree <- function(formula, data, model = c("tppcm", "gpcm", "pcm"),
                      parm = c("all", "disc", "diff", "slope", "int"),
                      param = c("irt", "si"), orthogonal = TRUE, impact = NULL, ...) {
  if (!requireNamespace("partykit", quietly = TRUE)) stop("tppcmtree() needs the partykit package", call. = FALSE)
  model <- match.arg(model); parm <- match.arg(parm); param <- match.arg(param)
  vars <- all.vars(formula); resp_name <- vars[1]
  y <- data[[resp_name]]
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y) || ncol(y) < 2 || !identical(deparse(formula[[2]]), resp_name))
    stop("the left-hand side of 'formula' must be one numeric matrix column of 'data' holding ",
         "all item responses, e.g. data$resp <- as.matrix(items); resp ~ z1 + z2", call. = FALSE)
  y <- .as_items(y)
  .check_coding(y)
  nas <- vapply(data[vars[-1]], function(v) sum(is.na(v)), 0)
  if (any(nas > 0))
    stop("the partitioning covariate(s) ", paste(names(nas)[nas > 0], collapse = ", "),
         " have missing values; drop those rows first", call. = FALSE)
  all_na <- rowSums(!is.na(y)) == 0
  if (any(all_na)) {
    message(sum(all_na), " row(s) without any item response dropped")
    data <- data[!all_na, , drop = FALSE]; y <- y[!all_na, , drop = FALSE]
  }
  K <- apply(as.matrix(y), 2, max, na.rm = TRUE) + 1                  # categories per item
  items <- colnames(y)
  if (is.null(items)) items <- paste0("I", seq_len(ncol(y)))
  I <- ncol(y); nsteps <- sum(K - 1)
  n_slope <- switch(model, tppcm = nsteps, gpcm = I, pcm = 1)
  block <- switch(parm, all = "all", disc = , slope = "slope", diff = , int = "location")
  parm_idx <- switch(block, all = NULL, slope = seq_len(n_slope), location = n_slope + seq_len(nsteps))
  if (!is.null(impact)) {
    g <- if (is.character(impact) && length(impact) == 1) data[[impact]] else impact
    if (is.null(g) || length(g) != nrow(data) || anyNA(g))
      stop("'impact' must name a complete factor in 'data' (or be one of length nrow(data))", call. = FALSE)
    y <- cbind(unname(as.matrix(y)), .impact = as.integer(factor(g)))
  }
  y <- unname(as.matrix(y))
  y[is.na(y)] <- -1                      # mob() drops rows with NA in the response; TAM keeps them
  data[[resp_name]] <- y

  fitter <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ...,
                     estfun = FALSE, object = FALSE) {
    y <- as.matrix(y)
    y[y < 0] <- NA
    group <- NULL
    if (!is.null(impact)) { group <- y[, ncol(y)]; y <- y[, -ncol(y), drop = FALSE] }
    colnames(y) <- items
    if (!is.null(weights) && any(weights != 1)) {
      rep_i <- rep(seq_len(nrow(y)), weights)
      y <- y[rep_i, , drop = FALSE]; group <- group[rep_i]
    }
    miss <- which(apply(y, 2, function(v) length(unique(stats::na.omit(v)))) < K)   # K per item
    if (length(miss)) stop("not all categories observed for item(s) ",
                           paste(colnames(y)[miss], collapse = ", "),
                           " in a node; increase minsize (e.g. minsize = 100) or use model = \"pcm\"",
                           call. = FALSE)
    if (!is.null(group) && length(unique(group)) == 1) group <- NULL
    fit <- TAM::tam.mml.3pl(y, E = step_design(y, model, K = K), group = group,
                            est.variance = FALSE, verbose = FALSE)
    info <- get_parts(fit)
    S <- .step_scores(info, param)
    cf <- .step_coef(info, param)
    names(cf) <- colnames(S)
    if (!is.null(impact)) {                       # name group parameters by the global levels
      lev_node <- sort(unique(group))
      ren <- function(nm) vapply(nm, function(n) {
        j <- as.integer(sub("grp:(mean|var)", "", n)); sub("[0-9]+$", lev_node[j], n) }, "")
      is_g <- grepl("^grp:", names(cf))
      names(cf)[is_g] <- ren(names(cf)[is_g]); colnames(S)[is_g] <- names(cf)[is_g]
      full <- paste0("grp:", c("mean", "var"), rep(seq_len(n_impact)[-1], each = 2))
      cf_g <- stats::setNames(rep(NA_real_, length(full)), full)
      cf_g[intersect(names(cf), full)] <- cf[intersect(names(cf), full)]
      cf <- c(cf[!is_g], cf_g)
    }
    sc <- if (estfun) .orthogonalize(S, block, orthogonal) else NULL
    list(coefficients = cf, objfun = -as.numeric(stats::logLik(fit)),
         estfun = sc, object = if (object) fit else NULL)
  }
  n_impact <- if (is.null(impact)) 0 else length(unique(y[, ncol(y)]))
  ctrl <- partykit::mob_control(parm = parm_idx, ytype = "matrix", ...)
  tr <- partykit::mob(formula, data = data, fit = fitter, control = ctrl)
  tr$info$tppcm <- list(model = model, parm = parm, param = param,
                        impact = if (is.character(impact)) impact else if (!is.null(impact)) "(vector)")
  class(tr) <- c("tppcmtree", class(tr))
  tr
}

#' @export
print.tppcmtree <- function(x, ...) {
  tp <- x$info$tppcm
  cat(sprintf("TPPCM tree | node model: %s | splits driven by: %s (%s parameterization)%s\n",
              tp$model, tp$parm, tp$param,
              if (is.null(tp$impact)) "" else sprintf(" | impact: %s", tp$impact)))
  cat(sprintf("Formula: %s\n", deparse(x$info$formula)))
  NextMethod(FUN = function(info) sprintf(": n = %d, -logLik = %.1f", info$nobs, info$objfun),
             header = FALSE, footer = FALSE)
  cat(sprintf("%d terminal node(s), %d parameters per node; coef(x) gives the node parameters\n",
              partykit::width(x), length(x$node$info$coefficients)))
  invisible(x)
}

# Efficient scores: residualize the tested block on the other block (OPG information)
.orthogonalize <- function(S, block, orthogonal) {
  if (!orthogonal || block == "all") return(S)
  test <- attr(S, "block") == block
  Iob <- crossprod(S[, !test, drop = FALSE])
  Iot <- crossprod(S[, !test, drop = FALSE], S[, test, drop = FALSE])
  S[, test] <- S[, test, drop = FALSE] - S[, !test, drop = FALSE] %*% solve(Iob, Iot)
  S
}
