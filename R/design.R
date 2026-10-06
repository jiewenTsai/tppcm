#' Step design for the TPPCM and its submodels
#'
#' Builds the slope design array `E` for [TAM::tam.mml.3pl()] so that the
#' slope parameters (`gammaslope`) are the step discriminations (CBDs)
#' \eqn{a_{il}}, restricted as described by `design`. The category slope of
#' category \eqn{k} is \eqn{\sum_{l \le k} a_{il}}.
#'
#' @section The design formula (an item x step layout on the log scale):
#' `design` describes the step discriminations of the items as a two-way
#' item x step layout, read on the log scale like a two-way ANOVA:
#' \deqn{\log a_{il} = \mu + \alpha_i + \gamma_l + (\alpha \gamma)_{il}.}
#'
#' | `design` | Model | \eqn{a_{il}} | Fitted with |
#' |---|---|---|---|
#' | `~ 1` | PCM (common slope) | \eqn{a} | `tam.mml.3pl(E = tppcm(dat, design = ~ 1))` |
#' | `~ item` | GPCM | \eqn{\alpha_i} | `tppcm(dat, design = ~ item)` |
#' | `~ step` | step model | \eqn{\gamma_l}, common to all items | `tppcm(dat, design = ~ step)` |
#' | `~ item + step` | Rank-1 (product form) | \eqn{\alpha_i \gamma_l} | [xxirt_tppcm()] (not linear) |
#' | `~ item * step` | saturated TPPCM (default) | \eqn{a_{il}} | `tppcm(dat)` |
#'
#' `~ item:step` and `~ item + step + item:step` are the same model as
#' `~ item * step`. The models form a hierarchy: `~ 1` is contained in
#' `~ item` and `~ step`, both are contained in `~ item + step`, which is
#' contained in `~ item * step`. Read as main effects and an interaction:
#'
#' * the item x step interaction means that items differ in the *shape* of
#'   their step discriminations: the ratios \eqn{a_{il} / a_{il'}} differ
#'   between items (saturated TPPCM);
#' * Rank-1 (main effects only): every item has the same shape
#'   \eqn{\gamma_l}, scaled to its own level \eqn{\alpha_i};
#' * GPCM (`~ item`): the shape is flat, every step of an item has the
#'   item's discrimination;
#' * `~ step`: the same shape and the same level for all items.
#'
#' Accordingly the score test of the GPCM against Rank-1
#' (`score_test(gpcm_fit, against = ~ item + step)`) asks whether there is a
#' step effect, and the test of Rank-1 against the saturated model
#' (`score_test(rank1_fit)`) whether there is an item x step interaction.
#' Rank-1 is multiplicative, so it cannot be written as a linear design for
#' TAM; `tppcm()` stops and points to [xxirt_tppcm()]. For the four linear
#' models the log scale does not matter (the same equality restrictions on
#' \eqn{a_{il}}), and `tppcm()` returns the design array. The log scale
#' presumes \eqn{a_{il} > 0}: Yu (1991) restricted the step discriminations
#' to be positive, TAM does not impose it ([xxirt_tppcm()] does, with
#' `lower0 = TRUE`).
#'
#' The same item x step language describes the locations in TAM's and
#' ConQuest's facet formulas: `formulaA = ~ item + step` is the rating scale
#' model (step parameters common to all items) and `~ item + item:step` the
#' partial credit model. `design = ~ step` is the rating-scale idea applied to
#' the discriminations.
#'
#' @section Custom restrictions:
#' Any other linear restriction on the step discriminations is written into
#' `index`: steps with the same index share one parameter. `index` takes
#' precedence over `design`. To fix a step, use the `gammaslope.fixed`
#' argument of [TAM::tam.mml.3pl()] with its index. Do not use
#' `gammaslope.constr.V` for equality constraints (it does not reach the
#' maximum likelihood estimate) and do not use `userfct.gammaslope` to impose
#' \eqn{a_{il} \ge 0} (truncation is not the constrained estimate); use
#' [xxirt_tppcm()] for bounds.
#'
#' @param dat Data frame or matrix of item responses coded `0, 1, ..., K_i - 1`.
#'   Items may have different numbers of categories; every category from 0
#'   to the item's maximum must be observed.
#' @param design One-sided formula in the terms `1`, `item`, `step` and
#'   `item:step` describing the step discriminations (see the section
#'   above): `~ 1`, `~ item`, `~ step` or `~ item * step` (default). Ignored
#'   when `index` is given.
#' @param index Optional integer matrix (items x steps, steps up to the
#'   largest number of steps) of parameter indices `1, ..., H`. Entries for
#'   steps an item does not have are ignored (use `NA`). Indices are renumbered
#'   to 1, ..., H in order. Overrides `design`.
#' @param K Number of categories per item (one number for all items, or one
#'   per item). Defaults to the observed maximum + 1 of each item; give it
#'   when a subsample may lack the highest categories (e.g. in tree nodes).
#' @return An array of dimension items x categories x 1 x H with attributes
#'   `"index"` (items x steps parameter numbers) and `"design"` (the design
#'   formula, `NULL` for a custom `index`).
#' @seealso [irt_pars()] for the estimates with correct standard errors;
#'   [xxirt_tppcm()] for Rank-1; [tppcm-package] for an overview and
#'   `vignette("tppcm-tutorial")`.
#' @references Yu, M.-N. (1991). *A two-parameter partial credit model*
#'   (Doctoral dissertation). University of Illinois at Urbana-Champaign.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- sim_tppcm(500, disc = matrix(c(1, 1.5, 0.7), 4, 3, byrow = TRUE),
#'                       diff = matrix(c(-1, 0, 1), 4, 3, byrow = TRUE))
#' fit <- function(E) TAM::tam.mml.3pl(dat, E = E, est.variance = FALSE, verbose = FALSE)
#' m_sat  <- fit(tppcm(dat))                       # saturated: ~ item * step
#' m_gpcm <- fit(tppcm(dat, design = ~ item))      # GPCM
#' m_step <- fit(tppcm(dat, design = ~ step))      # step model
#' m_pcm  <- fit(tppcm(dat, design = ~ 1))         # PCM (common slope)
#' score_test(m_gpcm, against = ~ item + step)     # step effect?
#' # items 1-2 saturated, items 3-4 GPCM
#' idx <- rbind(1:3, 4:6, rep(7, 3), rep(8, 3))
#' m_mix <- fit(tppcm(dat, index = idx))
#' }
#' @export
tppcm <- function(dat, design = ~ item * step, index = NULL, K = NULL) {
  code <- .parse_design(design)
  if (code == "rank1" && is.null(index))
    stop("design = ~ item + step (Rank-1, a_il = alpha_i * gamma_l) is multiplicative and cannot be ",
         "written as a linear design for TAM; fit it with xxirt_tppcm(dat, design = ~ item + step)",
         call. = FALSE)
  .design_array(dat, code, index, K)
}

# The design formulas of the five models (internal codes as names)
.design_codes <- c(pcm = "~ 1", gpcm = "~ item", step = "~ step", rank1 = "~ item + step",
                   tppcm = "~ item * step")
.design_names <- c(pcm = "PCM", gpcm = "GPCM", step = "step model", rank1 = "Rank-1 product form",
                   tppcm = "saturated TPPCM")

.design_formula <- function(code) {
  if (!code %in% names(.design_codes)) return(NULL)
  stats::as.formula(.design_codes[[code]], env = baseenv())
}

# design formula -> internal code "pcm", "gpcm", "step", "rank1", "tppcm"
.parse_design <- function(design, arg = "design") {
  if (is.character(design))
    stop("'", arg, "' must be a formula: ~ 1, ~ item, ~ step, ~ item + step or ~ item * step ",
         "(model names such as \"", design[1], "\" are no longer used)", call. = FALSE)
  if (!inherits(design, "formula") || length(design) != 2)
    stop("'", arg, "' must be a one-sided formula in item and step: ~ 1, ~ item, ~ step, ",
         "~ item + step or ~ item * step", call. = FALSE)
  tt <- stats::terms(design)
  lab <- attr(tt, "term.labels")
  lab[lab == "step:item"] <- "item:step"
  bad <- setdiff(lab, c("item", "step", "item:step"))
  if (length(bad) || !is.null(attr(tt, "offset")))
    stop("'", arg, "' may only contain the terms 1, item, step and item:step; not ",
         paste(c(bad, if (!is.null(attr(tt, "offset"))) "offset()"), collapse = ", "), call. = FALSE)
  if (attr(tt, "intercept") == 0)
    stop("'", arg, "' must keep the intercept (the overall level mu of log a_il); remove '- 1' or '+ 0'",
         call. = FALSE)
  if ("item:step" %in% lab) "tppcm"
  else if (all(c("item", "step") %in% lab)) "rank1"
  else if ("item" %in% lab) "gpcm"
  else if ("step" %in% lab) "step"
  else "pcm"
}

# Step design array for an internal design code ("pcm", "gpcm", "step",
# "tppcm") or a custom index
.design_array <- function(dat, code = "tppcm", index = NULL, K = NULL) {
  model <- code
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
  if (!all(dim(index) == c(I, L))) stop("index must be a ", I, " x ", L, " matrix (items x steps)", call. = FALSE)
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
  attr(E, "design") <- if (custom) NULL else .design_formula(model)
  E
}

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

.deparse_design <- function(f) .design_codes[[.parse_design(f)]]
