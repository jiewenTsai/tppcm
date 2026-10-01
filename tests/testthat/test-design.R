test_that("step_design reproduces the NRM and the GPCM in TAM", {
  dat <- make_data()
  E <- step_design(dat)
  expect_equal(dim(E), c(4, 4, 1, 12))
  m_sat <- TAM::tam.mml.3pl(dat, E = E, est.variance = FALSE, verbose = FALSE)
  m_nrm <- TAM::tam.mml.2pl(dat, irtmodel = "2PL", verbose = FALSE)
  expect_equal(as.numeric(logLik(m_sat)), as.numeric(logLik(m_nrm)), tolerance = 1e-3)
  m_g3 <- TAM::tam.mml.3pl(dat, E = step_design(dat, "gpcm"), est.variance = FALSE, verbose = FALSE)
  m_g2 <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  expect_equal(as.numeric(logLik(m_g3)), as.numeric(logLik(m_g2)), tolerance = 1e-3)
})

test_that("step_design checks its index", {
  dat <- make_data(100)
  expect_error(step_design(dat, index = matrix(1, 2, 3)))
  expect_equal(dim(step_design(dat, index = matrix(c(1, 3), 4, 3)))[4], 2)   # renumbered
})

test_that("step and pcm designs; tppcm() is an alias", {
  dat <- make_data()
  expect_identical(tppcm(dat, "step"), step_design(dat, "step"))
  expect_equal(attr(tppcm(dat, "step"), "index")[, 2], rep(2L, 4))
  m_pcmE <- TAM::tam.mml.3pl(dat, E = tppcm(dat, "pcm"), est.variance = FALSE, verbose = FALSE)
  m_pcm  <- TAM::tam.mml(dat, verbose = FALSE)
  expect_equal(as.numeric(logLik(m_pcmE)), as.numeric(logLik(m_pcm)), tolerance = 0.01)
  expect_equal(unname(m_pcmE$gammaslope), sqrt(m_pcm$variance[1]), tolerance = 0.01)
  st <- score_test(m_pcm, against = "step")
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
  expect_equal(dimnames(tppcm(d, "step"))[[4]], paste0("step", 1:3))
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
