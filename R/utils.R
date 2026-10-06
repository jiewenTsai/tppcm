# Category probabilities of one TPPCM item at the nodes in Theta (matrix, 1 column)
.tppcm_P <- function(a, d, Theta) {
  th <- Theta[, 1]; K <- length(a) + 1
  eta <- matrix(0, length(th), K)
  for (k in 2:K) eta[, k] <- eta[, k - 1] + a[k - 1] * th - d[k - 1]
  P <- exp(eta - apply(eta, 1, max))
  P / rowSums(P)
}

# One row per saturated parameter: item, itemnr, type ("a"/"d"), step, name.
# Order within an item: a_1..a_L, d_1..d_L (same as coef() of xxirt_tppcm()).
.tppcm_parindex <- function(items, ncat) {
  out <- do.call(rbind, lapply(seq_along(items), function(i) {
    L <- ncat[i] - 1
    data.frame(item = items[i], itemnr = i, type = rep(c("a", "d"), each = L),
               step = rep(seq_len(L), 2), stringsAsFactors = FALSE)
  }))
  out$name <- paste0(out$item, "_", out$type, out$step)
  out
}

.item_names <- function(dat) {
  nm <- colnames(dat)
  if (is.null(nm)) nm <- paste0("I", seq_len(ncol(dat)))
  nm
}

# Columns of J1 that are not in span(J0), after projecting out J0
.extra_dirs <- function(J0, J1, tol = 1e-7) {
  Q <- qr.Q(qr(J0))
  R <- J1 - Q %*% crossprod(Q, J1)
  keep <- sqrt(colSums(R^2)) > tol
  if (!any(keep)) stop("the alternative adds no direction to the restricted model", call. = FALSE)
  R <- R[, keep, drop = FALSE]
  qrR <- qr(R, tol = tol)
  R[, qrR$pivot[seq_len(qrR$rank)], drop = FALSE]
}

# Multiple-testing adjustment shared by all tables of tests: every table with
# one row per test has the columns p and p_adj (stats::p.adjust() over all rows)
# and names the method in its header.
.check_adjust <- function(adjust) {
  if (!is.character(adjust) || length(adjust) != 1 || !adjust %in% stats::p.adjust.methods)
    stop("'adjust' must be one of ", paste0('"', stats::p.adjust.methods, '"', collapse = ", "),
         " (see ?p.adjust)", call. = FALSE)
  adjust
}

.add_p_adj <- function(out, adjust) {
  if (!is.null(out) && nrow(out)) out$p_adj <- stats::p.adjust(out$p, adjust)
  out
}

.adj_label <- function(adjust) {
  lab <- c(holm = "Holm", hochberg = "Hochberg", hommel = "Hommel", bonferroni = "Bonferroni",
           BH = "Benjamini-Hochberg", BY = "Benjamini-Yekutieli", fdr = "Benjamini-Hochberg",
           none = "none (= p)")
  paste("p_adj:", lab[[adjust]])
}
