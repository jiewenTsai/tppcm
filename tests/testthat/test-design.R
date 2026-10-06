test_that("tppcm() designs reproduce the NRM and the GPCM in TAM", {
  dat <- make_data()
  E <- tppcm(dat)
  expect_equal(dim(E), c(4, 4, 1, 12))
  m_sat <- TAM::tam.mml.3pl(dat, E = E, est.variance = FALSE, verbose = FALSE)
  m_nrm <- TAM::tam.mml.2pl(dat, irtmodel = "2PL", verbose = FALSE)
  expect_equal(as.numeric(logLik(m_sat)), as.numeric(logLik(m_nrm)), tolerance = 1e-3)
  m_g3 <- TAM::tam.mml.3pl(dat, E = tppcm(dat, design = ~ item), est.variance = FALSE, verbose = FALSE)
  m_g2 <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  expect_equal(as.numeric(logLik(m_g3)), as.numeric(logLik(m_g2)), tolerance = 1e-3)
})

test_that("tppcm() checks its index", {
  dat <- make_data(100)
  expect_error(tppcm(dat, index = matrix(1, 2, 3)))
  expect_equal(dim(tppcm(dat, index = matrix(c(1, 3), 4, 3)))[4], 2)   # renumbered
})

test_that("step and pcm designs", {
  dat <- make_data()
  expect_equal(attr(tppcm(dat, design = ~ step), "index")[, 2], rep(2L, 4))
  m_pcmE <- TAM::tam.mml.3pl(dat, E = tppcm(dat, design = ~ 1), est.variance = FALSE, verbose = FALSE)
  m_pcm  <- TAM::tam.mml(dat, verbose = FALSE)
  expect_equal(as.numeric(logLik(m_pcmE)), as.numeric(logLik(m_pcm)), tolerance = 0.01)
  expect_equal(unname(m_pcmE$gammaslope), sqrt(m_pcm$variance[1]), tolerance = 0.01)
  st <- score_test(m_pcm, against = ~ step)
  expect_equal(st$df, 2)
})

test_that("coding checks and misuse messages", {
  dat <- as.matrix(make_data())
  expect_error(tppcm(dat + 1), "coded 0, 1")
  d2 <- dat; d2[d2[, 2] == 0, 2] <- 1
  expect_error(tppcm(d2), "not observed for item")
  expect_error(irt_pars(dat), "expected a fitted model")
  f <- TAM::tam.mml.3pl(dat, E = tppcm(dat), verbose = FALSE)        # variance not fixed
  expect_warning(get_parts(f), "est.variance = FALSE")
})

test_that("items with different numbers of categories", {
  set.seed(5)
  A <- matrix(c(1, 1.4, .8), 5, 3, byrow = TRUE); B <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
  d <- as.matrix(sim_tppcm(1200, A, B)); d[, 4] <- pmin(d[, 4], 2); d[, 5] <- pmin(d[, 5], 1)
  E <- tppcm(d)
  expect_equal(dim(E)[4], 3 * 3 + 2 + 1)
  expect_equal(dimnames(E)[[4]][10:12], c("I4_disc1", "I4_disc2", "I5_disc1"))
  f <- TAM::tam.mml.3pl(d, E = E, est.variance = FALSE, verbose = FALSE)
  n <- TAM::tam.mml.2pl(d, irtmodel = "2PL", verbose = FALSE)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(n)), tolerance = 1e-3)
  p <- irt_pars(f, fit = FALSE)
  expect_true(all(is.na(p$items["I5", c("disc2", "disc3")])))
  expect_equal(nrow(wald_test(f)), 4)                       # binary item skipped
  expect_equal(dimnames(tppcm(d, design = ~ step))[[4]], paste0("step", 1:3))
  d2 <- d; d2[d2[, 1] == 1, 1] <- 2
  expect_error(tppcm(d2), "middle categories")
})

test_that("sim_tppcm checks shapes and accepts named matrices and NA steps", {
  A <- matrix(c(1.0, 1.1, 1.2, 0.9, 1.0, 1.3, 1.1, 1.0, 1.2), 3, 3, dimnames = list(paste0("Q", 1:3), NULL))
  B <- matrix(c(-1, 0, 1, -0.5, 0.2, 0.9, -1.2, -0.1, 0.8), 3, 3, byrow = TRUE)
  expect_equal(dim(sim_tppcm(20, A, B)), c(20, 3))
  expect_error(sim_tppcm(10, matrix(1, 3, 2), matrix(0, 3, 3)), "same dimensions")
  d <- sim_tppcm(400, rbind(c(1, 1, 1), c(1, 1, NA)), rbind(c(-1, 0, 1), c(-1, 0, NA)))
  expect_lte(max(d[, 2]), 2)
  expect_error(sim_tppcm(10, rbind(c(NA, 1)), rbind(c(NA, 0))), "come last")
  expect_error(tppcm(data.frame(a = factor(c(0, 1)), b = c(0, 1))), "not numeric")
  expect_error(tppcm(cbind(id = 1:50, b = rep(0:1, 25))), "non-item column")
  expect_error(tppcm(cbind(a = rep(0, 10), b = rep(0:1, 5))), "single category")
})

test_that("design formulas are parsed as an item x step layout", {
  dat <- make_data(200)
  code <- function(f) tppcm:::.parse_design(f)
  expect_equal(code(~ 1), "pcm")
  expect_equal(code(~ item), "gpcm")
  expect_equal(code(~ step), "step")
  expect_equal(code(~ item + step), "rank1")
  expect_equal(code(~ step + item), "rank1")
  for (f in list(~ item * step, ~ item:step, ~ step:item, ~ item + step + item:step, ~ (item + step)^2))
    expect_equal(code(f), "tppcm")
  expect_identical(tppcm(dat, design = ~ item:step), tppcm(dat))
  expect_identical(tppcm(dat, design = ~ item + step + item:step), tppcm(dat, design = ~ item * step))
  expect_error(tppcm(dat, design = ~ item + step), "xxirt_tppcm\\(dat, design = ~ item \\+ step\\)")
  expect_error(tppcm(dat, design = ~ item + rater), "only contain the terms")
  expect_error(tppcm(dat, design = ~ log(item)), "only contain the terms")
  expect_error(tppcm(dat, design = ~ item - 1), "intercept")
  expect_error(tppcm(dat, design = y ~ item), "one-sided formula")
  expect_error(tppcm(dat, design = "gpcm"), "must be a formula")
  expect_error(tppcm(dat, model = "gpcm"), "unused argument")
  expect_false(exists("step_design", envir = asNamespace("tppcm")))
  expect_error(score_test(TAM::tam.mml(dat, verbose = FALSE), against = "gpcm"), "must be a formula")
  expect_equal(attr(tppcm(dat, design = ~ item), "design"), ~ item, ignore_attr = TRUE)
  expect_null(attr(tppcm(dat, index = matrix(1, 4, 3)), "design"))
})

test_that("each design formula equals the corresponding index design and its logLik", {
  dat <- make_data()
  idx <- list("~ 1" = matrix(1, 4, 3), "~ item" = matrix(1:4, 4, 3),
              "~ step" = matrix(1:3, 4, 3, byrow = TRUE), "~ item * step" = matrix(1:12, 4, 3, byrow = TRUE))
  fit <- function(E) TAM::tam.mml.3pl(dat, E = E, est.variance = FALSE, verbose = FALSE)
  for (f in names(idx)) {
    E1 <- tppcm(dat, design = stats::as.formula(f)); E2 <- tppcm(dat, index = idx[[f]])
    expect_equal(unname(unclass(E1)[seq_along(E1)]), unname(unclass(E2)[seq_along(E2)]))
    m <- fit(E1)
    expect_equal(as.numeric(logLik(m)), as.numeric(logLik(fit(E2))), tolerance = 1e-6)
    x <- get_parts(m)
    expect_equal(get_parts(x, "formula"), stats::as.formula(f), ignore_attr = TRUE)
    expect_match(tppcm:::.label(x), f, fixed = TRUE)
  }
  expect_null(get_parts(fit(tppcm(dat, index = rbind(1:3, 4:6, 7, 8))), "formula"))
  # TAM's own fits report their design too
  expect_equal(get_parts(TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE), "formula"), ~ item,
               ignore_attr = TRUE)
})

test_that("score_test, xxirt_tppcm, dsf_design and tppcmtree take design formulas", {
  dat <- make_data()
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  st <- score_test(g, against = ~ item + step)
  expect_equal(st$df, 2)
  expect_output(print(st), "against the Rank-1 product form \\(~ item \\+ step\\)")
  expect_equal(score_test(g)$df, score_test(g, against = ~ item:step)$df)
  r1 <- xxirt_tppcm(dat, design = ~ item + step)
  expect_equal(get_parts(r1, "model"), "rank1")
  expect_output(print(r1), "design ~ item \\+ step")
  expect_error(xxirt_tppcm(dat, design = ~ step), "tam.mml.3pl")
  expect_error(xxirt_tppcm(dat, model = "rank1"), "no 'model' argument")
  d <- dsf_design(dat, rep(1:2, 300), design = ~ item)
  expect_output(print(d), "design ~ item")
  expect_error(dsf_design(dat, rep(1:2, 300), design = ~ item + step), "multiplicative")
  expect_error(dsf_design(dat, rep(1:2, 300), model = "gpcm"), "unused argument")
})
