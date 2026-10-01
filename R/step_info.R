#' Step, item and test information
#'
#' Fisher information of a TPPCM-family fit, split into the contributions of
#' the steps. With the item score on the step scale
#' \eqn{A_i(X) = \sum_{l \le X} a_{il}}, the item information is
#' \eqn{I_i(\theta) = \mathrm{Var}(A_i(X) \mid \theta)} and step \eqn{l}
#' contributes
#' \deqn{c_{il}(\theta) = a_{il} \sum_m a_{im} \{S_{\max(l,m)} - S_l S_m\},
#'   \quad S_l = P(X_i \ge l \mid \theta),}
#' so that \eqn{\sum_l c_{il} = I_i}; the contributions are non-negative when
#' all \eqn{a_{il} \ge 0}. For the PCM (all \eqn{a = 1}) this is the usual
#' item information, the variance of the item score.
#'
#' @param x A fitted model (see [get_parts()]) or a `get_parts` object.
#' @param theta Trait values, on the scale of the fit.
#' @return Object of class `step_info` with `theta`, `step` (array
#'   theta x items x steps), `item` (theta x items) and `test`.
#'   `plot()` draws the test, item or step information.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(1000, disc = outer(c(0.8, 1, 1.2, 1.5, 1), c(1, 1.5, 0.7)),
#'                  diff = matrix(c(-1, 0, 1), 5, 3, byrow = TRUE))
#' fit <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
#' inf <- step_info(fit)
#' inf
#' plot(inf)                 # test and item information
#' plot(inf, item = "I4")    # step contributions of one item
#' }
#' @export
step_info <- function(x, theta = seq(-4, 4, by = 0.05)) {
  x <- get_parts(x)
  if (!is.numeric(theta) || !length(theta) || anyNA(theta))
    stop("'theta' must be a numeric vector of trait values", call. = FALSE)
  a <- x$a; b <- x$diff
  I <- nrow(a); L <- ncol(a)
  step <- array(NA_real_, c(length(theta), I, L),
                dimnames = list(NULL, rownames(a), colnames(a)))
  for (i in seq_len(I)) {
    ok <- which(!is.na(a[i, ])); ai <- a[i, ok]; bi <- b[i, ok]
    P <- .tppcm_P(ai, ai * bi, matrix(theta, ncol = 1))
    S <- t(apply(P, 1, function(p) rev(cumsum(rev(p)))))[, -1, drop = FALSE]   # P(X >= l)
    for (l in seq_along(ok)) {
      cov_lm <- sapply(seq_along(ok), function(m) S[, max(l, m)] - S[, l] * S[, m])
      step[, i, ok[l]] <- ai[l] * as.numeric(matrix(cov_lm, ncol = length(ok)) %*% ai)
    }
  }
  item <- apply(step, c(1, 2), sum, na.rm = TRUE)
  structure(list(theta = theta, step = step, item = item, test = rowSums(item)),
            class = "step_info")
}

#' @export
print.step_info <- function(x, digits = 3, ...) {
  cat(sprintf("Fisher information on %d theta values in [%g, %g]\n",
              length(x$theta), min(x$theta), max(x$theta)))
  cat("max: maximum of the item (test) information and where it is reached;\n",
      "max_step*: maximum of each step's contribution (each at its own theta)\n", sep = "")
  pk <- function(v) c(max = max(v), at = x$theta[which.max(v)])
  tab <- rbind(t(apply(x$item, 2, pk)), test = pk(x$test))
  steps <- apply(x$step, c(2, 3), function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE))
  out <- cbind(round(tab, digits), round(rbind(steps, NA), digits))
  colnames(out) <- c("max", "theta_at_max", paste0("max_", dimnames(x$step)[[3]]))
  out[nrow(out), -(1:2)] <- NA
  print(out, na.print = "")
  invisible(x)
}

#' @export
plot.step_info <- function(x, item = NULL, ...) {
  if (is.null(item)) {
    nI <- ncol(x$item)
    graphics::matplot(x$theta, cbind(x$test, x$item), type = "l", lty = c(1, rep(2, nI)),
                      lwd = c(2, rep(1, nI)), col = c(1, seq_len(nI) + 1),
                      xlab = expression(theta), ylab = "Information", ...)
    if (nI <= 10)
      graphics::legend("topright", c("test", colnames(x$item)), lty = c(1, rep(2, nI)),
                       col = c(1, seq_len(nI) + 1), bty = "n", cex = 0.8)
    else graphics::legend("topright", c("test", "items"), lty = 1:2, lwd = 2:1, bty = "n", cex = 0.8)
  } else {
    if (length(item) != 1 || !(item %in% colnames(x$item) || (is.numeric(item) && item %in% seq_len(ncol(x$item)))))
      stop("'item' must be one item name (", paste(colnames(x$item), collapse = ", "),
           ") or its number", call. = FALSE)
    if (is.numeric(item)) item <- colnames(x$item)[item]
    s <- x$step[, item, , drop = TRUE]
    graphics::matplot(x$theta, cbind(x$item[, item], s), type = "l", lty = c(1, rep(2, ncol(s))),
                      lwd = c(2, rep(1, ncol(s))), col = c(1, seq_len(ncol(s)) + 1),
                      xlab = expression(theta), ylab = "Information", main = item, ...)
    graphics::legend("topright", c("item", colnames(s)), lty = c(1, rep(2, ncol(s))),
                     col = c(1, seq_len(ncol(s)) + 1), bty = "n", cex = 0.8)
  }
  invisible(x)
}
