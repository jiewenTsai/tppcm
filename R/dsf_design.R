#' Multigroup TPPCM design with group.equal and group.partial
#'
#' Builds virtual items and the arguments of [TAM::tam.mml.3pl()] for a
#' multigroup TPPCM in which the step parameters are constrained across groups
#' in the way of the `group.equal` and `group.partial` arguments of lavaan.
#' Every item is split into one virtual item per group; a person answers only
#' the virtual items of their own group.
#'
#' The step parameters are named as in the CBD literature (see
#' [tppcm-package]):
#' `"cbd"` (category boundary discrimination, the step discrimination or
#' slope, cf. loadings) and `"ibd"` (intercept boundary difference, the step
#' intercept, cf. intercepts or thresholds). With the CBDs of a step equal,
#' equal IBDs and equal CBLs (category boundary locations, CBL = -IBD / CBD)
#' are the same restriction, so `"cbl"` is accepted for `"ibd"`.
#'
#' * `group.equal = NULL`: configural model (each group its own parameters;
#'   the latent mean and variance of every group are fixed at 0 and 1).
#' * `group.equal = "cbd"`: boundary-level metric model (CBDs equal, IBDs free;
#'   latent variances of groups 2, ..., G free, means fixed).
#' * `group.equal = c("cbd", "ibd")` (default): boundary-level scalar model
#'   (latent means and variances of groups 2, ..., G free).
#'
#' `group.partial` releases single steps from these restrictions, as
#' `list(cbd = "SO3:2", ibd = c("SO3:2", "SO5"))` (`"item:step"` or an item
#' name for all its steps); a character vector releases the steps from all
#' restrictions in `group.equal`. Releasing the CBD of a step while its IBD is
#' kept equal makes the boundary curves of the groups cross at `theta = 0`
#' (the mean of the reference group), which depends on the origin of the
#' scale; release both, or use `"cbl"`, in which case the IBD of a step with a
#' released CBD is released as well (with unequal CBDs, equal CBLs are not a
#' linear restriction).
#'
#' The reference group is `sort(unique(group))[1]`, as in TAM. The
#' identification of the latent means and variances is returned as
#' `variance.fixed` and `beta.fixed`; pass both to [TAM::tam.mml.3pl()]
#' together with `group`.
#'
#' @param dat Data frame or matrix of item responses coded `0, 1, ..., K_i - 1`.
#' @param group Vector of group memberships, one per row of `dat`, without
#'   missing values.
#' @param group.equal `NULL`, `"cbd"`, or `c("cbd", "ibd")` (default);
#'   `"cbl"` may be used for `"ibd"`.
#' @param group.partial Steps released from `group.equal`: a list with
#'   elements `cbd` and/or `ibd` (or `cbl`), or a character vector for all
#'   restrictions.
#' @param design Structure of the step discriminations within each group, a
#'   linear design formula of [tppcm()]: `~ item * step` (default), `~ item`,
#'   `~ step` or `~ 1`.
#' @param index Optional items x steps matrix of discrimination indices within
#'   a group, as in [tppcm()]; overrides `design`.
#' @param K Number of categories per item, as in [tppcm()].
#' @return A list of class `tppcm_dsf_design` with `resp` (virtual items,
#'   named `item_group`), `A`, `E`, `variance.fixed`, `beta.fixed` (arguments
#'   of [TAM::tam.mml.3pl()]), `index` and `index_int` (virtual items x steps
#'   parameter numbers of the CBDs and IBDs), `groups`, `equal`, `partial`.
#' @seealso [dsf()] for likelihood ratio tests step by step.
#' @examples
#' \donttest{
#' set.seed(1)
#' a <- matrix(c(1, 1.5, 1.2), 5, 3, byrow = TRUE)
#' b <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
#' g <- rep(1:2, each = 500)
#' dat <- rbind(sim_tppcm(500, disc = a, diff = b),
#'              sim_tppcm(500, disc = a, diff = b, theta = rnorm(500, 0.3, 0.8)))
#' fit_mg <- function(d) TAM::tam.mml.3pl(d$resp, A = d$A, E = d$E, group = g,
#'   variance.fixed = d$variance.fixed, beta.fixed = d$beta.fixed, verbose = FALSE)
#' configural <- fit_mg(dsf_design(dat, g, group.equal = NULL))
#' metric     <- fit_mg(dsf_design(dat, g, group.equal = "cbd"))
#' scalar     <- fit_mg(dsf_design(dat, g, group.equal = c("cbd", "ibd")))
#' anova(configural, metric)
#' anova(metric, scalar)
#' partial <- fit_mg(dsf_design(dat, g, group.partial = list(cbd = "I2:2", ibd = "I2:2")))
#' anova(scalar, partial)
#' }
#' @export
dsf_design <- function(dat, group, group.equal = c("cbd", "ibd"), group.partial = NULL,
                       design = ~ item * step, index = NULL, K = NULL) {
  model <- .parse_design(design)
  if (model == "rank1" && is.null(index))
    stop("design = ~ item + step (Rank-1) is multiplicative and has no linear multigroup design; ",
         "use ~ item * step, ~ item, ~ step or ~ 1", call. = FALSE)
  dat <- .as_items(dat)
  N <- nrow(dat); I <- ncol(dat); items <- .item_names(dat)
  if (length(group) != N) stop("group must have one value per row of dat", call. = FALSE)
  if (anyNA(group)) stop("group has missing values", call. = FALSE)
  groups <- sort(unique(group)); G <- length(groups)
  if (G < 2) stop("group must have at least two values", call. = FALSE)
  gi <- match(group, groups)
  ge <- .dsf_types(group.equal, "group.equal")
  use_cbl <- "cbl" %in% group.equal
  if (ge[["ibd"]] && !ge[["cbd"]] && use_cbl)
    stop("equal CBLs with unequal CBDs are not a linear restriction; use group.equal = c(\"cbd\", \"cbl\")",
         call. = FALSE)
  E0 <- .design_array(dat, model, index = index, K = K)
  base <- attr(E0, "index"); L <- ncol(base); has <- !is.na(base)
  Kv <- rowSums(has) + 1
  part <- .dsf_partial(group.partial, items, has, ge)
  if (use_cbl) part$ibd <- part$ibd | (part$cbd & ge[["ibd"]])
  eq_cbd <- if (ge[["cbd"]]) has & !part$cbd else has & FALSE
  eq_ibd <- if (ge[["ibd"]]) has & !part$ibd else has & FALSE
  gl <- as.character(groups)
  vitem <- rep(seq_len(I), each = G); vgrp <- rep(seq_len(G), I)
  vnames <- paste0(items[vitem], "_", gl[vgrp])
  resp <- matrix(NA_real_, N, I * G, dimnames = list(NULL, vnames))
  for (j in seq_along(vitem)) { rows <- gi == vgrp[j]; resp[rows, j] <- dat[rows, vitem[j]] }
  # CBD indices: the reference group uses the base design; a step of group g > 1
  # shares it when equal; otherwise it gets a new index (one per base index and
  # group in the configural/IBD-only case, so that within-group restrictions of
  # the base design are kept; one per released step otherwise)
  H0 <- max(base, na.rm = TRUE); anames0 <- dimnames(E0)[[4]]
  idx <- base[vitem, , drop = FALSE]; anames <- anames0; key <- list()
  for (j in seq_along(vitem)) if (vgrp[j] > 1) for (l in which(has[vitem[j], ])) {
    if (eq_cbd[vitem[j], l]) next
    k <- paste(base[vitem[j], l], vgrp[j], if (ge[["cbd"]]) paste(vitem[j], l) else "")
    if (is.null(key[[k]])) {
      key[[k]] <- H0 + length(key) + 1L
      anames <- c(anames, if (ge[["cbd"]]) paste0(items[vitem[j]], "_cbd", l, "_", gl[vgrp[j]])
                          else paste0(anames0[base[vitem[j], l]], "_", gl[vgrp[j]]))
    }
    idx[j, l] <- key[[k]]
  }
  E <- .design_array(resp, index = idx, K = Kv[vitem])
  dimnames(E)[[4]] <- anames
  # IBD indices: one per existing step, new ones for unequal steps of groups g > 1
  bi <- matrix(NA_integer_, L, I); bi[t(has)] <- seq_len(sum(has)); bi <- t(bi)
  dnames <- paste0(items[row(bi)[has]], "_ibd", col(bi)[has])[order(bi[has])]
  idx_int <- bi[vitem, , drop = FALSE]
  nxt <- max(bi, na.rm = TRUE)
  for (j in seq_along(vitem)) if (vgrp[j] > 1) for (l in which(has[vitem[j], ])) {
    if (eq_ibd[vitem[j], l]) next
    nxt <- nxt + 1; idx_int[j, l] <- nxt
    dnames <- c(dnames, paste0(items[vitem[j]], "_ibd", l, "_", gl[vgrp[j]]))
  }
  A <- array(0, c(I * G, L + 1, nxt), list(vnames, paste0("Category", 0:L), dnames))
  for (j in seq_along(vitem)) {
    Lj <- Kv[vitem[j]] - 1
    for (k in seq_len(Lj)) for (l in seq_len(k)) A[j, k + 1, idx_int[j, l]] <- -1
    if (Lj < L) A[j, (Lj + 2):(L + 1), ] <- NA
  }
  # identification: variances of groups > 1 free if some CBD is shared, means if some IBD is shared
  vg <- if (any(eq_cbd)) 1L else seq_len(G)
  mg <- if (any(eq_ibd)) 1L else seq_len(G)
  sn <- paste0("step", seq_len(L))
  idx_E <- attr(E, "index"); dimnames(idx_E) <- dimnames(idx_int) <- list(vnames, sn)
  dimnames(part$cbd) <- dimnames(part$ibd) <- list(items, sn)
  out <- list(resp = as.data.frame(resp), A = A, E = E,
              variance.fixed = unname(cbind(vg, 1, 1, 1)), beta.fixed = unname(cbind(mg, 1, 0)),
              index = idx_E, index_int = idx_int, groups = groups,
              equal = c(cbd = ge[["cbd"]], ibd = ge[["ibd"]]), partial = part, items = items,
              design = if (is.null(index)) .design_formula(model) else NULL)
  class(out) <- "tppcm_dsf_design"
  out
}

# group.equal -> c(cbd = TRUE/FALSE, ibd = TRUE/FALSE)
.dsf_types <- function(x, arg) {
  bad <- setdiff(x, c("cbd", "ibd", "cbl"))
  if (length(bad)) stop(arg, ": unknown parameter type(s) ", paste0("'", bad, "'", collapse = ", "),
                        "; use \"cbd\", \"ibd\" or \"cbl\"", call. = FALSE)
  c(cbd = "cbd" %in% x, ibd = any(c("ibd", "cbl") %in% x))
}

# group.partial -> list(cbd = logical items x steps, ibd = ...)
.dsf_partial <- function(x, items, has, ge) {
  out <- list(cbd = has & FALSE, ibd = has & FALSE)
  if (is.null(x)) return(out)
  if (!is.list(x)) x <- stats::setNames(rep(list(x), sum(ge)), c("cbd", "ibd")[ge])
  if (is.null(names(x))) stop("group.partial must be a list with elements cbd and/or ibd (or cbl)", call. = FALSE)
  names(x)[names(x) == "cbl"] <- "ibd"
  if (length(setdiff(names(x), c("cbd", "ibd"))))
    stop("group.partial must be a list with elements cbd and/or ibd (or cbl)", call. = FALSE)
  for (nm in names(x)) out[[nm]] <- out[[nm]] | .dsf_steps(x[[nm]], items, has, "group.partial")
  out
}

# "item:step" / item names / logical matrix -> logical items x steps matrix
.dsf_steps <- function(x, items, has, arg) {
  out <- matrix(FALSE, nrow(has), ncol(has))
  if (is.null(x)) return(out)
  if (is.matrix(x) && is.logical(x)) {
    if (!all(dim(x) == dim(has))) stop(arg, " must be a ", nrow(has), " x ", ncol(has), " matrix", call. = FALSE)
    out[] <- x & has
    return(out)
  }
  for (s in as.character(x)) {
    p <- strsplit(s, ":", fixed = TRUE)[[1]]
    i <- match(p[1], items)
    if (is.na(i)) stop(arg, ": unknown item '", p[1], "'", call. = FALSE)
    if (length(p) == 1) { out[i, ] <- has[i, ]; next }
    l <- suppressWarnings(as.integer(p[2]))
    if (is.na(l) || l < 1 || l > ncol(has) || !has[i, l])
      stop(arg, ": item ", p[1], " has no step ", p[2], call. = FALSE)
    out[i, l] <- TRUE
  }
  out
}

#' @export
print.tppcm_dsf_design <- function(x, ...) {
  G <- length(x$groups)
  eq <- c("cbd", "ibd")[x$equal]
  lvl <- if (!length(eq)) "configural" else if (identical(eq, "cbd")) "metric (CBD equal)"
         else if (identical(eq, "ibd")) "IBD equal" else "scalar (CBD and IBD equal)"
  cat(sprintf("Multigroup TPPCM design: %d items x %d groups (reference group: %s)\n",
              length(x$items), G, format(x$groups[1])))
  cat(sprintf("  group.equal: %s -> %s; within groups: %s\n",
              if (length(eq)) paste(eq, collapse = ", ") else "none", lvl,
              if (is.null(x$design)) "custom index" else paste("design", .deparse_design(x$design))))
  lab <- function(m) { w <- which(m, arr.ind = TRUE)
    if (!nrow(w)) "none" else paste0(x$items[w[, 1]], ":", w[, 2], collapse = ", ") }
  cat("  group.partial cbd:", lab(x$partial$cbd), "\n")
  cat("  group.partial ibd:", lab(x$partial$ibd), "\n")
  cat(sprintf("  %d CBD and %d IBD parameters; free latent means: groups %s; variances: groups %s\n",
              dim(x$E)[4], dim(x$A)[3],
              .free_lab(G, x$beta.fixed[, 1]), .free_lab(G, x$variance.fixed[, 1])))
  cat("Fit with TAM::tam.mml.3pl(d$resp, A = d$A, E = d$E, group = <group>,\n",
      "        variance.fixed = d$variance.fixed, beta.fixed = d$beta.fixed)\n", sep = "")
  invisible(x)
}

.free_lab <- function(G, fixed) {
  f <- setdiff(seq_len(G), fixed)
  if (length(f)) paste(f, collapse = ", ") else "none"
}
