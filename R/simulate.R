#' Simulate item responses from a TPPCM
#'
#' @param N Number of persons.
#' @param disc Items x steps matrix of step discriminations \eqn{a_{il}}.
#'   Trailing `NA`s give items with fewer steps; a vector gives binary items.
#' @param diff Items x steps matrix of step difficulties \eqn{b_{il}}; the
#'   kernel of step \eqn{l} is \eqn{a_{il}(\theta - b_{il})}. Same shape as
#'   `disc`.
#' @param theta Person parameters (default standard normal).
#' @return Data frame of responses coded `0, ..., K-1`, columns `I1, I2, ...`.
#' @examples
#' dat <- sim_tppcm(200, disc = matrix(1, 3, 2), diff = matrix(c(-0.5, 0.5), 3, 2, byrow = TRUE))
#' @export
sim_tppcm <- function(N, disc, diff, theta = stats::rnorm(N)) {
  a <- as.matrix(disc); b <- as.matrix(diff)
  if (!identical(dim(a), dim(b))) stop("disc and diff must have the same dimensions (items x steps)", call. = FALSE)
  if (any(is.na(a) != is.na(b))) stop("disc and diff must have NA in the same places", call. = FALSE)
  if (length(theta) != N) stop("theta must have length N (", N, ")", call. = FALSE)
  I <- nrow(a)
  dat <- matrix(NA_integer_, N, I)
  for (i in 1:I) {
    ok <- unname(which(!is.na(a[i, ])))
    if (!identical(ok, seq_along(ok))) stop("NA steps must come last (item ", i, ")", call. = FALSE)
    K <- length(ok) + 1
    P <- .tppcm_P(a[i, ok], a[i, ok] * b[i, ok], matrix(theta, ncol = 1))
    cum <- t(apply(P, 1, cumsum))[, -K, drop = FALSE]
    dat[, i] <- rowSums(stats::runif(N) > cum)
  }
  colnames(dat) <- paste0("I", 1:I)
  as.data.frame(dat)
}
