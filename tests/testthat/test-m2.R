test_that("m2 matches mirt for dichotomous 2PL", {
  skip_if_not_installed("mirt")
  set.seed(5)
  a <- c(.8, 1, 1.2, 1.5, 1, .9); b <- seq(-1, 1, length.out = 6)
  th <- rnorm(1200)
  X <- sapply(1:6, function(i) rbinom(1200, 1, plogis(a[i] * (th - b[i]))))
  colnames(X) <- paste0("I", 1:6)
  nodes <- seq(-6, 6, length.out = 41)
  f <- TAM::tam.mml.2pl(X, irtmodel = "2PL", verbose = FALSE,
                        control = list(nodes = nodes, conv = 1e-7, deviance.conv = 1e-9, maxiter = 3000))
  mm <- mirt::mirt(X, 1, "2PL", verbose = FALSE, quadpts = 41, theta_lim = c(-6, 6), TOL = 1e-7)
  ours <- m2(f); ref <- mirt::M2(mm, type = "M2", quadpts = 41, theta_lim = c(-6, 6))
  expect_equal(ours$df, ref$df)
  expect_equal(ours$M2, ref$M2[1], tolerance = 1e-3)
  expect_equal(ours$SRMSR, ref$SRMSR, tolerance = 1e-3)
})

test_that("m2 degrees of freedom, types and fit header", {
  dat <- as.matrix(make_data(800))
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  r <- m2(g)
  expect_equal(r$type, "M2")
  expect_equal(r$df, 4 * 3 + choose(4, 2) * 9 - (4 + 12))     # margins - parameters
  expect_equal(m2(g, "C2")$df, 4 * 3 + choose(4, 2) - 16)
  expect_error(m2(g, "M2*"), "degrees of freedom")
  expect_true(r$CFI <= 1 && r$RMSEA >= 0)
  f <- irt_pars(g)$fit
  expect_equal(f$npar, 16)
  expect_equal(f$AIC, CDM::IRT.IC(g)[["AIC"]], tolerance = 1e-6)
  expect_output(print(irt_pars(g)), "M2\\(")
  d2 <- dat; d2[1, 1] <- NA
  r2 <- m2(TAM::tam.mml.2pl(d2, irtmodel = "GPCM", verbose = FALSE))    # one missing response
  expect_equal(r2$df, r$df)
  expect_equal(r2$M2, r$M2, tolerance = 0.05)
})
