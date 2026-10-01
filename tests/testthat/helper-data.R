make_data <- function(N = 600, seed = 11) {
  set.seed(seed)
  a <- outer(c(0.8, 1.1, 1.4, 1.0), c(1, 1.5, 0.7))
  d <- matrix(c(-1, 0, 1, -0.5, 0.3, 1.2, -1.2, -0.2, 0.8, -0.8, 0.5, 1.5), 4, byrow = TRUE)
  sim_tppcm(N, a, d / a)        # d given as a * diff
}

# marginal log-likelihood of the saturated TPPCM on TAM's default nodes
ll_sat <- function(par, dat, nodes = seq(-6, 6, length.out = 21)) {
  X <- as.matrix(dat); L <- max(X) ; w <- dnorm(nodes); w <- w / sum(w)
  lik <- matrix(1, nrow(X), length(nodes))
  for (i in seq_len(ncol(X))) {
    pp <- par[(i - 1) * 2 * L + 1:(2 * L)]
    P <- tppcm:::.tppcm_P(pp[1:L], pp[L + 1:L], matrix(nodes, ncol = 1))
    lik <- lik * t(P[, X[, i] + 1])
  }
  sum(log(lik %*% w))
}
