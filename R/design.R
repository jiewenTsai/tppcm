#' Step design for estimating the TPPCM with TAM::tam.mml.3pl()
#'
#' `tppcm()` and `step_design()` are the same function.
#'
#' Builds the slope design array `E` for [TAM::tam.mml.3pl()] so that the
#' slope parameters (`gammaslope`) are the step discriminations \eqn{a_{il}}.
#' The category slope of category \eqn{k} is \eqn{\sum_{l \le k} a_{il}}.
#'
#' Linear restrictions on the step discriminations are written into `index`:
#' steps with the same index share one parameter. To fix a step, use the
#' `gammaslope.fixed` argument of [TAM::tam.mml.3pl()] with its index.
#' Do not use `gammaslope.constr.V` for equality constraints (it does not
#' reach the maximum likelihood estimate) and do not use
#' `userfct.gammaslope` to impose \eqn{a_{il} \ge 0} (truncation is not the
#' constrained estimate); use [xxirt_tppcm()] for bounds.
#'
#' @param dat Data frame or matrix of item responses coded `0, 1, ..., K_i - 1`.
#'   Items may have different numbers of categories; every category from 0
#'   to the item's maximum must be observed.
#' @param model Shortcut for common designs; ignored when `index` is given.
#'   * `"tppcm"`: every step its own discrimination (saturated TPPCM).
#'   * `"gpcm"`: one discrimination per item, \eqn{a_{il} = \alpha_i}.
#'   * `"step"`: one discrimination per step, shared by all items,
#'     \eqn{a_{il} = \gamma_l} (columns of `index` equal).
#'   * `"pcm"`: one common discrimination, \eqn{a_{il} = a}. With the latent
#'     variance fixed at 1 this is the PCM with free variance
#'     (same log-likelihood as [TAM::tam.mml()]).
#' @param index Optional integer matrix (items x steps, steps up to the
#'   largest number of steps) of parameter indices `1, ..., H`. Entries for
#'   steps an item does not have are ignored (use `NA`). Indices are renumbered
#'   to 1, ..., H in order.
#' @param K Number of categories per item (one number for all items, or one
#'   per item). Defaults to the observed maximum + 1 of each item; give it
#'   when a subsample may lack the highest categories (e.g. in tree nodes).
#' @return An array of dimension items x categories x 1 x H with attribute
#'   `"index"`.
#' @seealso [irt_pars()] for the estimates with correct standard errors;
#'   [tppcm-package] for an overview and `vignette("tppcm-tutorial")`.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(500, disc = matrix(c(1, 1.5, 0.7), 4, 3, byrow = TRUE),
#'                       diff = matrix(c(-1, 0, 1), 4, 3, byrow = TRUE))
#' # saturated TPPCM
#' mod <- TAM::tam.mml.3pl(dat, E = step_design(dat), est.variance = FALSE,
#'                         verbose = FALSE)
#' # same thing, other name
#' mod <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
#' # items 1-2 saturated, items 3-4 GPCM
#' idx <- rbind(1:3, 4:6, rep(7, 3), rep(8, 3))
#' mod2 <- TAM::tam.mml.3pl(dat, E = step_design(dat, index = idx),
#'                          est.variance = FALSE, verbose = FALSE)
#' }
#' @export
step_design <- function(dat, model = c("tppcm", "gpcm", "step", "pcm"), index = NULL,
                        K = NULL) {
  model <- match.arg(model)
  custom <- !is.null(index)
  dat <- .as_items(dat)
  I <- ncol(dat)
  if (is.null(K)) {
    .check_coding(dat)
    K <- apply(dat, 2, max, na.rm = TRUE) + 1
  }
  K <- unname(rep(K, length.out = I))
  Li <- K - 1; L <- max(Li)
  has <- outer(seq_len(I), seq_len(L), function(i, l) l <= Li[i])     # existing steps
  if (is.null(index)) {
    index <- switch(model,
                    tppcm = { tm <- matrix(NA_integer_, L, I); tm[t(has)] <- seq_len(sum(has)); t(tm) },
                    gpcm  = matrix(rep(seq_len(I), L), I, L),
                    step  = matrix(seq_len(L), I, L, byrow = TRUE),
                    pcm   = matrix(1L, I, L))
  }
  index <- as.matrix(index)
  if (!all(dim(index) == c(I, L))) stop("index must be a ", I, " x ", L, " matrix", call. = FALSE)
  index[!has] <- NA
  if (anyNA(index[has])) stop("index has missing entries for existing steps", call. = FALSE)
  used <- sort(unique(index[has]))
  index[has] <- match(index[has], used)                               # renumber 1..H
  H <- length(used)
  E <- array(0, c(I, L + 1, 1, H))
  for (i in seq_len(I)) for (l in seq_len(Li[i])) {
    cats <- (l + 1):(Li[i] + 1)
    E[i, cats, 1, index[i, l]] <- E[i, cats, 1, index[i, l]] + 1
  }
  items <- .item_names(dat)
  pnames <- if (custom) paste0("a", seq_len(H)) else switch(model,
    tppcm = paste0(items[row(index)[has]], "_disc", col(index)[has])[order(index[has])],
    gpcm  = paste0(items, "_disc"),
    step  = paste0("step", seq_len(H)),
    pcm   = "disc")
  dimnames(E) <- list(items, paste0("Cat", 0:L), "Dim1", pnames)
  attr(E, "index") <- index
  E
}

#' @rdname step_design
#' @export
tppcm <- step_design

# Categories must be coded 0, 1, ..., K_i - 1 with every category observed
.check_coding <- function(dat) {
  if (!is.numeric(dat)) stop("responses must be numeric codes 0, 1, ..., K-1 ",
                             "(convert factors or character columns first)", call. = FALSE)
  nm <- .item_names(dat)
  empty <- colSums(!is.na(dat)) == 0
  if (any(empty)) stop("no responses for item(s) ", paste(nm[empty], collapse = ", "), call. = FALSE)
  mins <- apply(dat, 2, min, na.rm = TRUE)
  maxs <- apply(dat, 2, max, na.rm = TRUE)
  if (any(maxs > 20)) stop("item(s) ", paste(nm[maxs > 20], collapse = ", "), " have more than 20 ",
                           "categories; is this an ID or another non-item column?", call. = FALSE)
  if (all(mins > 0)) {
    stop("categories must be coded 0, 1, ..., K-1, but the lowest observed value is ",
         min(mins), "; recode, e.g. dat - ", min(mins), call. = FALSE)
  }
  if (any(mins > 0)) {
    stop("category 0 is not observed for item(s) ", paste(.item_names(dat)[mins > 0], collapse = ", "),
         "; collapse or recode the categories", call. = FALSE)
  }
  if (any(maxs == 0)) stop("item(s) ", paste(nm[maxs == 0], collapse = ", "),
                           " have a single category; drop them", call. = FALSE)
  gaps <- apply(dat, 2, function(v) { v <- v[!is.na(v)]; !all(seq.int(0, max(v)) %in% v) })
  if (any(gaps)) {
    stop("unobserved middle categories for item(s) ", paste(.item_names(dat)[gaps], collapse = ", "),
         "; collapse or recode the categories", call. = FALSE)
  }
  invisible(TRUE)
}

# Data frame or matrix of item responses -> numeric matrix (factors refused)
.as_items <- function(dat) {
  if (is.data.frame(dat)) {
    bad <- !vapply(dat, is.numeric, logical(1))
    if (any(bad)) stop("column(s) ", paste(names(dat)[bad], collapse = ", "), " are not numeric; ",
                       "code responses as numbers 0, 1, ..., K-1", call. = FALSE)
  }
  as.matrix(dat)
}
