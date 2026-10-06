#' Limited-information goodness of fit: M2, RMSEA2, SRMSR, CFI, TLI
#'
#' Computes the limited-information statistic of Maydeu-Olivares and Joe
#' (2006) for any fit that [get_parts()] accepts, including step designs and
#' the product form, which other software cannot express. The derivatives of
#' the margins are taken with respect to the saturated parameters and mapped
#' to the fitted model with the Jacobian of [get_parts()].
#'
#' * `"M2"`: all univariate and bivariate margins (categories 1, ..., K-1).
#' * `"C2"`: univariate margins and bivariate cross-products
#'   \eqn{E(X_i X_j)} (Cai & Monroe, 2014).
#' * `"M2*"`: means \eqn{E(X_i)} and cross-products \eqn{E(X_i X_j)}.
#'
#' The reduced statistics have few degrees of freedom; for the saturated
#' TPPCM with few items they are not identified (`df <= 0`), so the full
#' `"M2"` is preferred when it is affordable. CFI and TLI compare against the
#' independence model (items independent, free category proportions). SRMSR
#' is the root mean squared difference between observed and model-implied
#' Pearson correlations of the item scores.
#'
#' **Missing responses.** The observed margins use the available cases and
#' their covariance is scaled by the numbers of persons observed on each pair
#' of margins. This is valid only when responses are missing completely at
#' random (MCAR), e.g. planned missingness by booklet. When non-response is
#' related to ability or to the responses themselves (omitted items, not
#' reached items), the statistic is invalid: in a simulation with
#' missingness depending on the response to one item, M2 rejected the true
#' model in every replication. Score it or analyze complete cases in that
#' situation. N in RMSEA is the number of persons.
#'
#' Under misfit `"C2"` can differ from mirt's C2 by 10-20% with the same
#' log-likelihood: both use the same margins and derivatives, but this
#' function uses the exact model-implied covariance matrix of the margins
#' (Maydeu-Olivares & Joe, 2006), whereas mirt estimates parts of it from the
#' sample. Under the true model the two agree.
#'
#' @param x A fitted model (see [get_parts()]) or the object returned by `get_parts()`.
#' @param type `"auto"` (`"M2"` with up to 2000 margins, otherwise `"C2"`),
#'   `"M2"`, `"C2"` or `"M2*"`. The computations use all univariate and
#'   bivariate margins, so the number of margins is limited to 5000 (about 30
#'   items with 4 categories).
#' @param CI Coverage of the RMSEA confidence interval.
#' @return Object of class `tppcm_m2`: a one-row data frame with `M2`, `df`,
#'   `p`, `RMSEA`, `RMSEA_lower`, `RMSEA_upper`, `SRMSR`, `CFI`, `TLI`.
#' @references Maydeu-Olivares, A., & Joe, H. (2006). Limited information
#'   goodness-of-fit testing in multidimensional contingency tables.
#'   *Psychometrika, 71*, 713-732.
#'
#'   Cai, L., & Monroe, S. (2014). A new statistic for evaluating item
#'   response theory models for ordinal data. CRESST Report 839.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(1000, disc = outer(c(0.8, 1, 1.2, 1.5, 1), c(1, 1.5, 0.7)),
#'                  diff = matrix(c(-1, 0, 1), 5, 3, byrow = TRUE))
#' m2(tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE))
#' m2(tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE))
#' }
#' @export
m2 <- function(x, type = c("auto", "M2", "C2", "M2*"), CI = 0.90) {
  type <- match.arg(type)
  x <- get_parts(x)
  dat <- as.matrix(x$data)
  if (!is.null(x$groups) && nrow(x$groups) > 1) stop("m2() supports single-group fits only", call. = FALSE)
  if (isTRUE(x$weighted)) stop("m2() does not support person weights", call. = FALSE)
  if (isTRUE(x$latreg)) stop("m2() does not support latent regression", call. = FALSE)
  N <- nrow(dat); I <- ncol(dat); ncat <- x$ncat
  w <- x$weights
  mg <- .margins(ncat)
  if (nrow(mg) > 5000) stop("too many margins (", nrow(mg), ") for m2(); at most 5000 (about 30 items with 4 categories)", call. = FALSE)
  obs <- .observed_margins(dat, mg)               # observed proportions and their scaling
  if (any(obs$n == 0)) {                          # pairs never observed together
    keep <- obs$n > 0
    mg <- mg[keep, ]; obs <- list(p = obs$p[keep], n = obs$n[keep],
                                  scale = if (is.matrix(obs$scale)) obs$scale[keep, keep] else obs$scale)
  }
  if (type == "auto") type <- if (nrow(mg) <= 2000) "M2" else "C2"
  T <- .reduce(mg, type, ncat)

  # model: margins, their covariance and derivatives
  P <- lapply(seq_len(I), function(i) t(x$probs[i, seq_len(ncat[i]), , drop = TRUE]))   # Q x K_i
  mod <- .margin_moments(P, mg, w)
  dP <- .dprob(P, x$theta, ncat)
  D <- .margin_jacobian(P, dP, mg, w, x$parindex) %*% x$J
  st <- .m2_stat(T, obs$p, mod$pi, mod$Xi * obs$scale, D)
  # independence model: category proportions free
  P0 <- lapply(seq_len(I), function(i) {
    v <- dat[!is.na(dat[, i]), i]; matrix(tabulate(v + 1, ncat[i]) / length(v), 1)
  })
  m0 <- .margin_moments(P0, mg, 1)
  D0 <- .null_jacobian(P0, mg)
  st0 <- .m2_stat(T, obs$p, m0$pi, m0$Xi * obs$scale, D0)

  rmsea <- function(X2, df) sqrt(max(X2 - df, 0) / (N * df))
  ci <- .rmsea_ci(st$stat, st$df, N, CI)
  # SRMSR from model-implied moments of the item scores
  sr <- .srmsr(dat, P, w, ncat)
  d_m <- max(st$stat - st$df, 0); d_0 <- max(st0$stat - st0$df, d_m)
  out <- data.frame(type = type, M2 = st$stat, df = st$df,
                    p = stats::pchisq(st$stat, st$df, lower.tail = FALSE),
                    RMSEA = rmsea(st$stat, st$df), RMSEA_lower = ci[1], RMSEA_upper = ci[2],
                    SRMSR = sr,
                    CFI = if (d_0 > 0) 1 - d_m / d_0 else 1,
                    TLI = (st0$stat / st0$df - st$stat / st$df) / (st0$stat / st0$df - 1))
  structure(out, class = c("tppcm_m2", "data.frame"), CI = CI, null = st0,
            missing = anyNA(dat))
}

#' @export
print.tppcm_m2 <- function(x, digits = 3, ...) {
  y <- as.data.frame(unclass(x))
  cat(sprintf("%s = %.2f, df = %d, p = %s\n", y$type, y$M2, y$df, format.pval(y$p, digits = digits)))
  cat(sprintf("RMSEA = %.*f [%.*f, %.*f], SRMSR = %.*f, CFI = %.*f, TLI = %.*f\n",
              digits, y$RMSEA, digits, y$RMSEA_lower, digits, y$RMSEA_upper,
              digits, y$SRMSR, digits, y$CFI, digits, y$TLI))
  if (isTRUE(attr(x, "missing")))
    cat("(missing responses: margins from available cases, valid under MCAR only)\n")
  invisible(x)
}

# --- internals ---------------------------------------------------------------

# Univariate (i, k) and bivariate (i, k, j, l) margins, categories >= 1
.margins <- function(ncat) {
  I <- length(ncat)
  uni <- do.call(rbind, lapply(seq_len(I), function(i)
    data.frame(i = i, k = seq_len(ncat[i] - 1), j = NA_integer_, l = NA_integer_)))
  bi <- do.call(rbind, lapply(seq_len(I - 1), function(i) do.call(rbind, lapply((i + 1):I, function(j)
    expand.grid(k = seq_len(ncat[i] - 1), l = seq_len(ncat[j] - 1), i = i, j = j)[, c("i", "k", "j", "l")]))))
  rbind(uni, bi)
}

# Linear map from the full margins to the margins used by the statistic
.reduce <- function(mg, type, ncat) {
  n <- nrow(mg); uni <- is.na(mg$j)
  if (type == "M2") return(diag(n))
  rows <- list()
  if (type == "C2") {
    for (r in which(uni)) { v <- numeric(n); v[r] <- 1; rows[[length(rows) + 1]] <- v }
  } else {
    for (i in unique(mg$i[uni])) {
      v <- numeric(n); sel <- uni & mg$i == i; v[sel] <- mg$k[sel]; rows[[length(rows) + 1]] <- v
    }
  }
  pairs <- unique(mg[!uni, c("i", "j")])
  for (r in seq_len(nrow(pairs))) {
    v <- numeric(n); sel <- !uni & mg$i == pairs$i[r] & mg$j == pairs$j[r]
    v[sel] <- mg$k[sel] * mg$l[sel]; rows[[length(rows) + 1]] <- v
  }
  do.call(rbind, rows)
}

# Margin probabilities at the nodes (margins x Q)
.margin_nodes <- function(P, mg) {
  uni <- is.na(mg$j)
  F <- matrix(0, nrow(mg), nrow(P[[1]]))
  for (r in seq_len(nrow(mg))) {
    f <- P[[mg$i[r]]][, mg$k[r] + 1]
    if (!uni[r]) f <- f * P[[mg$j[r]]][, mg$l[r] + 1]
    F[r, ] <- f
  }
  F
}

# Model-implied margins pi and covariance Xi = E(uu') - pi pi'
.margin_moments <- function(P, mg, w) {
  F <- .margin_nodes(P, mg)
  pi <- as.numeric(F %*% w)
  E <- tcrossprod(F * rep(sqrt(w), each = nrow(F)))          # assumes disjoint items
  cells <- rbind(cbind(mg$i, mg$k), cbind(mg$j, mg$l))
  idx <- rep(seq_len(nrow(mg)), 2)
  ok <- !is.na(cells[, 1])
  cells <- cells[ok, , drop = FALSE]; idx <- idx[ok]
  # margins sharing item i with the same category k: count P_ik once
  for (key in unique(paste(cells[, 1], cells[, 2]))) {
    sel <- paste(cells[, 1], cells[, 2]) == key
    m <- unique(idx[sel]); ik <- cells[which(sel)[1], ]
    Pik <- P[[ik[1]]][, ik[2] + 1]
    E[m, m] <- tcrossprod(F[m, , drop = FALSE] * rep(sqrt(w / Pik), each = length(m)))
  }
  # same item, different categories: impossible
  for (i in unique(cells[, 1])) {
    sel <- cells[, 1] == i
    m <- idx[sel]; k <- cells[sel, 2]
    E[m, m][outer(k, k, "!=")] <- 0
  }
  diag(E) <- pi
  list(pi = pi, Xi = E - tcrossprod(pi))
}

# dP_ik / d(a_i1..L, d_i1..L) at the nodes: list of arrays Q x K_i x 2L
.dprob <- function(P, theta, ncat) {
  lapply(seq_along(P), function(i) {
    Pi <- P[[i]]; L <- ncat[i] - 1; Q <- nrow(Pi)
    S <- t(apply(Pi, 1, function(p) rev(cumsum(rev(p)))))[, -1, drop = FALSE]   # Q x L, P(X >= l)
    out <- array(0, c(Q, L + 1, 2 * L))
    for (k in 0:L) for (l in seq_len(L)) {
      g <- Pi[, k + 1] * ((k >= l) - S[, l])
      out[, k + 1, l] <- theta * g
      out[, k + 1, L + l] <- -g
    }
    out
  })
}

# d pi / d saturated parameters (margins x p)
.margin_jacobian <- function(P, dP, mg, w, pidx) {
  D <- matrix(0, nrow(mg), nrow(pidx))
  cols <- split(seq_len(nrow(pidx)), pidx$itemnr)
  for (r in seq_len(nrow(mg))) {
    i <- mg$i[r]; k <- mg$k[r]
    if (is.na(mg$j[r])) {
      D[r, cols[[i]]] <- colSums(dP[[i]][, k + 1, ] * w)
    } else {
      j <- mg$j[r]; l <- mg$l[r]
      D[r, cols[[i]]] <- colSums(dP[[i]][, k + 1, ] * (w * P[[j]][, l + 1]))
      D[r, cols[[j]]] <- colSums(dP[[j]][, l + 1, ] * (w * P[[i]][, k + 1]))
    }
  }
  D
}

# Independence model: derivatives with respect to the category proportions 1..K-1
.null_jacobian <- function(P0, mg) {
  off <- cumsum(c(0, sapply(P0, ncol) - 1))
  D <- matrix(0, nrow(mg), off[length(off)])
  for (r in seq_len(nrow(mg))) {
    i <- mg$i[r]; k <- mg$k[r]
    if (is.na(mg$j[r])) D[r, off[i] + k] <- 1 else {
      j <- mg$j[r]; l <- mg$l[r]
      D[r, off[i] + k] <- P0[[j]][1, l + 1]; D[r, off[j] + l] <- P0[[i]][1, k + 1]
    }
  }
  D
}

# Observed margins from the available cases. With missing data (assumed MCAR),
# Cov(p_a, p_b) = Xi_ab * n_ab / (n_a n_b), where n_ab counts persons observed on
# all items of both margins; with complete data the scale is 1 / N.
.observed_margins <- function(dat, mg) {
  N <- nrow(dat)
  if (!anyNA(dat)) {
    p <- sapply(seq_len(nrow(mg)), function(r) {
      u <- dat[, mg$i[r]] == mg$k[r]
      if (!is.na(mg$j[r])) u <- u & dat[, mg$j[r]] == mg$l[r]
      mean(u)
    })
    return(list(p = p, n = rep(N, nrow(mg)), scale = 1 / N))
  }
  R <- !is.na(dat)
  O <- sapply(seq_len(nrow(mg)), function(r) {
    o <- R[, mg$i[r]]
    if (!is.na(mg$j[r])) o <- o & R[, mg$j[r]]
    o
  })
  U <- sapply(seq_len(nrow(mg)), function(r) {
    u <- dat[, mg$i[r]] %in% mg$k[r]
    if (!is.na(mg$j[r])) u <- u & dat[, mg$j[r]] %in% mg$l[r]
    u
  })
  n <- colSums(O)
  p <- ifelse(n > 0, colSums(U & O) / pmax(n, 1), 0)
  Nab <- crossprod(O * 1)
  list(p = p, n = n, scale = Nab / outer(pmax(n, 1), pmax(n, 1)))
}

# Sigma: covariance of the observed margins (Xi / N with complete data)
.m2_stat <- function(T, p, pi, Sigma, D) {
  e <- T %*% (p - pi); Xr <- T %*% Sigma %*% t(T); Dr <- T %*% D
  # C2 = Dc (Dc' Xi Dc)^-1 Dc' with Dc the orthogonal complement of the derivatives
  q <- qr(Dr)
  Q <- qr.Q(q, complete = TRUE)
  if (q$rank >= ncol(Q)) stop("M2 cannot be computed: too few degrees of freedom", call. = FALSE)
  Dc <- Q[, (q$rank + 1):ncol(Q), drop = FALSE]
  C2 <- Dc %*% solve(crossprod(Dc, Xr %*% Dc), t(Dc))
  df <- ncol(Dc)
  list(stat = as.numeric(crossprod(e, C2 %*% e)), df = df)
}

.rmsea_ci <- function(X2, df, N, CI) {
  f <- function(lambda, q) stats::pchisq(X2, df, ncp = lambda) - q
  bound <- function(q) {
    if (f(0, q) < 0) return(0)
    hi <- max(X2, 1); while (f(hi, q) > 0) hi <- hi * 2
    stats::uniroot(f, c(0, hi), q = q)$root
  }
  lam <- c(bound((1 + CI) / 2), bound((1 - CI) / 2))
  sqrt(lam / (N * df))
}

.srmsr <- function(dat, P, w, ncat) {
  I <- ncol(dat)
  Ex  <- sapply(seq_len(I), function(i) P[[i]] %*% (0:(ncat[i] - 1)))      # Q x I
  Ex2 <- sapply(seq_len(I), function(i) P[[i]] %*% (0:(ncat[i] - 1))^2)
  m <- colSums(Ex * w); v <- colSums(Ex2 * w) - m^2
  C <- crossprod(Ex * sqrt(w)) - tcrossprod(m)
  diag(C) <- v
  R_mod <- stats::cov2cor(C); R_obs <- stats::cor(dat, use = "pairwise.complete.obs")
  sqrt(mean((R_obs - R_mod)[lower.tri(R_obs)]^2))
}
