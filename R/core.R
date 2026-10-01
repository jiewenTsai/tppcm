# Internal mathematics. Saturated parameter order within an item:
# a_1..a_L, d_1..d_L (L = K - 1 steps).

# Person scores, gradient and information (Louis, OPG) of the saturated TPPCM,
# evaluated at the item response probabilities of any fitted TPPCM-family model.
# grp (optional): list(g = group index per person, MEAN, VAR per group,
# free = groups whose mean and variance are free parameters); their scores and
# information are appended after the item parameters (names "grp:mean<g>",
# "grp:var<g>"). The prior of group g is the normal density at the nodes,
# normalized to sum to one, as in TAM.
.score_core <- function(dat, post, probs, theta, ncat, weights = NULL, grp = NULL) {
  dat <- as.matrix(dat)
  N <- nrow(dat); I <- ncol(dat); TP <- length(theta)
  if (is.null(weights)) weights <- rep(1, N)
  post <- as.matrix(post)
  pidx <- .tppcm_parindex(.item_names(dat), ncat)
  p <- nrow(pidx)
  cols <- split(seq_len(p), pidx$itemnr)
  resp <- !is.na(dat)
  ng <- if (is.null(grp)) 0 else 2 * length(grp$free)
  gnames <- if (ng) paste0("grp:", c("mean", "var"), rep(grp$free, each = 2)) else character(0)
  ginfo <- if (ng) lapply(grp$free, function(g) {
    mu <- grp$MEAN[g]; v <- grp$VAR[g]
    phi <- stats::dnorm(theta, mu, sqrt(v)); w <- phi / sum(phi)
    s <- cbind((theta - mu) / v, -1 / (2 * v) + (theta - mu)^2 / (2 * v^2))   # d log phi_q
    H <- cbind(-1 / v, -(theta - mu) / v^2, 1 / (2 * v^2) - (theta - mu)^2 / v^3)  # d2 log phi_q: mm, mv, vv
    sbar <- colSums(w * s)
    EwH <- matrix(colSums(w * H)[c(1, 2, 2, 3)], 2)
    Varws <- crossprod(s * sqrt(w)) - tcrossprod(sbar)
    list(sel = which(grp$g == g), dlogw = sweep(s, 2, sbar), negH = -H, const = EwH + Varws)
  })

  scores  <- matrix(0, N, p + ng)
  VarPart <- matrix(0, p + ng, p + ng)
  EPart   <- matrix(0, p + ng, p + ng)
  for (q in seq_len(TP)) {
    th <- theta[q]
    wq <- post[, q] * weights
    Uq <- matrix(0, N, p + ng)
    for (j in seq_along(ginfo)) {
      gi <- ginfo[[j]]; gc <- p + 2 * j - 1:0
      Uq[gi$sel, gc] <- rep(gi$dlogw[q, ], each = length(gi$sel))
      EPart[gc, gc] <- EPart[gc, gc] + sum(wq[gi$sel]) * (matrix(gi$negH[q, c(1, 2, 2, 3)], 2) + gi$const)
    }
    for (i in seq_len(I)) {
      L  <- ncat[i] - 1
      Pi <- probs[i, seq_len(L + 1), q]
      S  <- rev(cumsum(rev(Pi)))[-1]                 # S_l = P(X >= l)
      r  <- resp[, i]
      C  <- outer(ifelse(r, dat[, i], 0), seq_len(L), ">=") * 1
      Dm <- sweep(C, 2, S)
      Dm[!r, ] <- 0
      ci <- cols[[i]]
      Uq[, ci] <- cbind(th * Dm, -Dm)                # complete-data score
      Cm <- outer(seq_len(L), seq_len(L), function(l, m) S[pmax(l, m)]) - tcrossprod(S)
      Hi <- rbind(cbind(th^2 * Cm, -th * Cm), cbind(-th * Cm, Cm))
      EPart[ci, ci] <- EPart[ci, ci] + sum(wq[r]) * Hi
    }
    scores  <- scores + post[, q] * Uq               # Fisher's identity
    VarPart <- VarPart + crossprod(Uq * sqrt(wq))
  }
  info_opg   <- crossprod(scores * sqrt(weights))
  info_louis <- EPart - VarPart + info_opg           # Louis' formula
  nms <- c(pidx$name, gnames)
  dimnames(info_louis) <- dimnames(info_opg) <- list(nms, nms)
  grad <- colSums(scores * weights)
  names(grad) <- colnames(scores) <- nms
  list(scores = scores, grad = grad, info = info_louis, info_opg = info_opg,
       parindex = pidx, nobs = sum(weights), gnames = gnames)
}

# Pad a matrix of directions in the saturated item-parameter space with zero
# rows for the group parameters appended by .score_core()
.pad_rows <- function(M, n) {
  if (nrow(M) >= n) return(M)
  rbind(M, matrix(0, n - nrow(M), ncol(M), dimnames = list(NULL, colnames(M))))
}

# log(P_k / P_{k-1}) = a_k theta - d_k, fitted over the nodes
# Stops when the logits are not linear in theta, i.e. the model is not a
# divide-by-total (TPPCM-family) model (e.g. guessing parameters).
.par_from_probs <- function(probs, theta, ncat, items, tol = 1e-3) {
  worst <- 0
  par <- unlist(lapply(seq_len(dim(probs)[1]), function(i) {
    L <- ncat[i] - 1
    ad <- sapply(seq_len(L), function(l) {
      y  <- log(probs[i, l + 1, ]) - log(probs[i, l, ])
      ok <- is.finite(y) & abs(theta) <= 4
      m  <- stats::lm(y[ok] ~ theta[ok])
      worst <<- max(worst, abs(stats::resid(m)))
      cf <- stats::coef(m)
      c(cf[2], -cf[1])
    })
    c(ad[1, ], ad[2, ])
  }))
  if (worst > tol)
    stop("the fitted item response functions are not of the TPPCM (divide-by-total) form ",
         "(e.g. guessing parameters); this model is not supported", call. = FALSE)
  names(par) <- .tppcm_parindex(items, ncat)$name
  par
}

# Jacobian of the saturated parameters with respect to the restricted ones
.jacobian <- function(model, parindex, par, E = NULL, fixed = NULL, loc = NULL) {
  p <- nrow(parindex)
  ia <- which(parindex$type == "a"); id <- which(parindex$type == "d")
  items <- unique(parindex$item)
  unit_cols <- function(sel) {
    M <- diag(p)[, sel, drop = FALSE]; colnames(M) <- parindex$name[sel]; M
  }
  col_by <- function(f, labels) matrix(sapply(labels, f), nrow = p)

  if (model == "E") {
    if (is.null(E)) stop("model 'E' needs the design array E")
    Ja <- matrix(0, p, dim(E)[4])
    for (r in ia) {
      i <- parindex$itemnr[r]; l <- parindex$step[r]
      Ja[r, ] <- E[i, l + 1, 1, ] - E[i, l, 1, ]
    }
    colnames(Ja) <- if (!is.null(dimnames(E)[[4]])) dimnames(E)[[4]] else paste0("a", seq_len(ncol(Ja)))
    if (length(fixed)) Ja <- Ja[, -fixed, drop = FALSE]
  } else if (model == "tppcm") {
    Ja <- unit_cols(ia)
  } else if (model == "gpcm") {
    Ja <- col_by(function(it) as.numeric(parindex$item == it & parindex$type == "a"), items)
    colnames(Ja) <- paste0(items, "_alpha")
  } else if (model == "step") {                      # a_il = gamma_l, shared by all items
    Ja <- col_by(function(l) as.numeric(parindex$type == "a" & parindex$step == l),
                 sort(unique(parindex$step)))
    colnames(Ja) <- paste0("step", sort(unique(parindex$step)))
  } else if (model == "pcm") {                       # common slope = free latent variance
    Ja <- matrix(as.numeric(parindex$type == "a"), p, 1, dimnames = list(NULL, "slope"))
  } else if (model == "pcm_fixed") {
    Ja <- NULL
  } else if (model == "rank1") {
    L <- max(parindex$step)
    if (length(unique(table(parindex$itemnr))) > 1) stop("rank1 needs equal category numbers")
    A <- matrix(par[ia], ncol = L, byrow = TRUE)
    alpha <- A[, 1]; gamma <- colMeans(A / alpha)
    Jal <- col_by(function(it) {
      v <- numeric(p); sel <- parindex$item == it & parindex$type == "a"
      v[sel] <- gamma[parindex$step[sel]]; v }, items)
    Jga <- col_by(function(l) {
      v <- numeric(p); sel <- parindex$type == "a" & parindex$step == l
      v[sel] <- alpha[parindex$itemnr[sel]]; v }, 2:L)
    colnames(Jal) <- paste0(items, "_alpha"); colnames(Jga) <- paste0("gamma", 2:L)
    Ja <- cbind(Jal, Jga)
  } else stop("unknown model: ", model)
  J <- cbind(Ja, if (is.null(loc)) unit_cols(id) else loc)
  rownames(J) <- parindex$name
  J
}

# Curvature term sum_k g_k d2h_k of the product form (zero for linear models)
.curvature <- function(model, J, parindex, grad) {
  curv <- matrix(0, ncol(J), ncol(J), dimnames = list(colnames(J), colnames(J)))
  if (model == "rank1") {
    for (r in which(parindex$type == "a" & parindex$step >= 2)) {
      ca <- paste0(parindex$item[r], "_alpha"); cg <- paste0("gamma", parindex$step[r])
      curv[ca, cg] <- curv[ca, cg] + grad[r]
      curv[cg, ca] <- curv[cg, ca] + grad[r]
    }
  }
  curv
}

.detect_model <- function(fit) {
  if (inherits(fit, "tam.mml.mfr") || (is.list(fit) && !is.null(fit[["facets"]])))
    stop("many-facet models (tam.mml.mfr) are not supported", call. = FALSE)
  if (inherits(fit, "tam.jml"))
    stop("JML fits (tam.jml) are not supported; use tam.mml(), tam.mml.2pl() or tam.mml.3pl()", call. = FALSE)
  if (inherits(fit, "tam.mml.3pl")) return("E")
  if (inherits(fit, "tam.mml")) {
    if (identical(fit$irtmodel, "2PL")) return("tppcm")
    if (identical(fit$irtmodel, "GPCM")) return("gpcm")
    if (!is.null(fit$variance.fixed)) return("pcm_fixed")
    return("pcm")
  }
  if (inherits(fit, "xxirt")) {
    pt <- fit$partable
    if (any(grepl("^RANK1TPPCM", pt$type))) return("rank1")
    pa <- pt[grepl("^a[0-9]+$", pt$parname), ]
    if (!any(pa$est)) return("pcm_fixed")
    pa <- pa[pa$est, ]
    if (all(tapply(pa$parindex, pa$item, function(x) length(unique(x))) == 1)) return("gpcm")
    return("tppcm")
  }
  stop("expected a fitted model from TAM (tam.mml(), tam.mml.2pl(), tam.mml.3pl()) or ",
       "xxirt_tppcm(), not an object of class '", class(fit)[1], "'", call. = FALSE)
}
