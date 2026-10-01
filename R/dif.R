#' Score-based DIF tests for steps, items or parameter blocks
#'
#' Tests whether parameters of a single-group fit vary with a person
#' covariate, without fitting a multiple-group model (Merkle & Zeileis, 2013;
#' Strobl, Kopf & Zeileis, 2015). The case-wise scores (`get_parts(fit, "estfun")`) are
#' ordered by the covariate and aggregated with a functional from
#' \pkg{strucchange}. Fit the saturated TPPCM (`tam.mml.3pl(E = tppcm(dat))`)
#' to test single step discriminations; a GPCM fit gives item-level tests.
#'
#' Each tested unit is replaced by its efficient scores (residuals from the
#' regression on all other parameters, outer-product information), so DIF in
#' other parameters, e.g. in the difficulties, does not leak into the test.
#'
#' **Impact is not DIF.** A single-group fit assumes one latent distribution
#' for everybody. When the covariate is related to the latent trait (groups
#' that differ in ability, a numeric covariate correlated with it), the scores
#' of all parameters move with the covariate and most items are flagged. Fit
#' the multigroup model (`tam.mml.3pl(..., group = )`) or the latent regression
#' (`Y = `) with that covariate first and test on that fit; the group means and
#' variances (or regression coefficients) are then part of the model, and only
#' item-level differences remain in the scores.
#'
#' @param x A fitted model (see [get_parts()]) or a `get_parts` object.
#' @param covariate Person covariate in the order of the fitted responses
#'   (no missing values). A factor gives the LM test (`"LM"`), an ordered factor
#'   the weighted maximum LM test (`"maxLMo"`), a numeric covariate the
#'   Cramer-von Mises type test (`"CvM"`).
#' @param parm Parameters to test: `"disc"` / `"slope"` (slope block),
#'   `"diff"` / `"int"` (location block), `"all"`, or column names of
#'   `get_parts(fit, "estfun")`.
#' @param by Testing unit: every parameter (`"param"`), every item (`"item"`)
#'   or all selected parameters together (`"block"`).
#' @param functional `NULL` (chosen from the covariate type), one of `"LM"`,
#'   `"maxLMo"`, `"DM"`, `"CvM"`, `"maxLM"`, or a \pkg{strucchange}
#'   functional.
#' @param param `"irt"` (`disc`, `diff`) or `"si"` (`slope`, `int`).
#' @return Data frame (class `tppcm_table`) with `unit`, `npar`, `statistic`,
#'   `p` and Holm-adjusted `p_holm`, in the order of the items.
#' @references Merkle, E. C., & Zeileis, A. (2013). Tests of measurement
#'   invariance without subgroups: A generalization of classical methods.
#'   *Psychometrika, 78*, 59-82.
#' @examples
#' \donttest{
#' set.seed(1)
#' a <- matrix(c(1, 1.5, 0.7), 5, 3, byrow = TRUE)
#' b <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
#' a2 <- a; a2[3, 2] <- 0.4                       # step DIF: item 3, step 2
#' dat <- rbind(sim_tppcm(800, a, b), sim_tppcm(800, a2, b))
#' g <- factor(rep(1:2, each = 800))
#' fit <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
#' dif_test(fit, g, parm = "disc")                # every step discrimination
#' dif_test(fit, g, parm = "disc", by = "item")
#' dif_test(fit, g, parm = "diff", by = "item")
#' }
#' @export
dif_test <- function(x, covariate, parm = "disc", by = c("param", "item", "block"),
                     functional = NULL, param = c("irt", "si")) {
  if (!requireNamespace("strucchange", quietly = TRUE)) stop("dif_test() needs the strucchange package", call. = FALSE)
  by <- match.arg(by); param <- match.arg(param)
  x <- get_parts(x)
  if (isTRUE(x$n_loc > 0)) stop("dif_test() needs free step locations; RSM-type fits are not supported", call. = FALSE)
  if (isTRUE(x$weighted)) stop("dif_test() does not support person weights", call. = FALSE)
  S <- .step_scores(x, param)
  if (length(covariate) != nrow(S)) stop("'covariate' must have one value per person of the fit", call. = FALSE)
  if (anyNA(covariate)) stop("'covariate' has missing values", call. = FALSE)
  if (is.character(covariate) || is.logical(covariate)) covariate <- factor(covariate)
  if (is.factor(covariate) && (is.null(x$groups) || nrow(x$groups) == 1) && !isTRUE(x$latreg))
    message("single-group fit: a difference in ability between the groups of 'covariate' (impact) ",
            "is reported as DIF in the difficulties; refit with group = covariate if that is possible (see ?dif_test)")
  blk <- attr(S, "block")
  S <- sweep(S, 2, colMeans(S))                       # scores sum to zero at the MLE
  tested <- if (length(parm) == 1 && parm %in% c("disc", "slope", "diff", "int", "all")) {
    switch(parm, disc = , slope = which(blk == "slope"), diff = , int = which(blk == "location"),
           all = seq_len(ncol(S)))
  } else match(parm, colnames(S))
  if (anyNA(tested)) stop("unknown parameters in 'parm'", call. = FALSE)
  nm <- colnames(S)[tested]
  it <- sub("_.*$", "", nm)
  units <- switch(by,
    param = stats::setNames(as.list(tested), nm),
    item  = split(tested, factor(it, levels = unique(it))),   # keep the item order
    block = list(all = tested))
  fun_name <- if (is.character(functional) || is.null(functional)) .dif_functional_name(covariate, functional) else "user"
  out <- do.call(rbind, lapply(names(units), function(u) {
    t <- units[[u]]
    St <- .efficient(S, t)
    gp <- strucchange::gefp(St, fit = NULL, scores = function(m) m, order.by = covariate)
    fn <- if (fun_name == "user") functional else .dif_functional(fun_name, gp)
    st <- strucchange::sctest(gp, functional = fn)
    data.frame(unit = u, npar = length(t), statistic = unname(st$statistic),
               p = unname(st$p.value))
  }))
  out$p_holm <- stats::p.adjust(out$p, "holm")
  rownames(out) <- NULL
  .as_table(out, sprintf("Score-based DIF test (%s), unit: %s", fun_name, by))
}

# Efficient scores of the columns t given all other columns (OPG information)
.efficient <- function(S, t) {
  St <- S[, t, drop = FALSE]
  if (length(t) == ncol(S)) return(St)
  So <- S[, -t, drop = FALSE]
  St - So %*% solve(crossprod(So), crossprod(So, St))
}

.dif_functional_name <- function(z, functional) {
  if (!is.null(functional)) return(match.arg(functional, c("LM", "maxLMo", "DM", "CvM", "maxLM")))
  if (is.ordered(z)) "maxLMo" else if (is.factor(z) || is.character(z) || is.logical(z)) "LM" else "CvM"
}

.dif_functional <- function(name, gp) {
  switch(name,
         LM     = strucchange::catL2BB(gp),
         maxLMo = strucchange::ordwmax(gp),
         DM     = strucchange::maxBB,
         CvM    = strucchange::meanL2BB,
         maxLM  = strucchange::supLM(0.1))
}

