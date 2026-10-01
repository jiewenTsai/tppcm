test_that("stored score and information match numerical derivatives", {
  skip_if_not_installed("numDeriv")
  dat <- make_data()
  fit <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  tp <- get_parts(fit)
  expect_s3_class(tp, "tppcm_parts")
  expect_identical(get_parts(tp), tp)
  expect_equal(ll_sat(tp$par, dat), as.numeric(logLik(fit)), tolerance = 1e-6)
  expect_equal(unname(tp$grad), numDeriv::grad(ll_sat, tp$par, dat = dat), tolerance = 1e-4)
  H <- numDeriv::hessian(ll_sat, tp$par, dat = dat)
  expect_lt(max(abs(tp$info + H)) / max(abs(H)), 1e-4)
  expect_equal(dim(tp$a), c(4, 3))
  expect_equal(ncol(tp$J), 4 + 12)
})

test_that("model detection", {
  dat <- make_data(300)
  expect_equal(get_parts(TAM::tam.mml(dat, verbose = FALSE))$model, "pcm")
  expect_equal(get_parts(TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE))$model, "gpcm")
  expect_equal(get_parts(xxirt_tppcm(dat, "rank1"))$model, "rank1")
  expect_equal(get_parts(xxirt_tppcm(dat, "gpcm"))$model, "gpcm")
})

test_that("standard errors agree with the numerical Hessian of xxirt", {
  dat <- make_data()
  x_sat <- xxirt_tppcm(dat, "tppcm")
  v <- sqrt(diag(xxirt_vcov(x_sat)))
  expect_equal(unname(irt_pars(x_sat, long = TRUE)$se[irt_pars(x_sat, long = TRUE)$par == "disc"]), unname(v[grep("_a", names(v))]), tolerance = 1e-3)

  x_r1 <- xxirt_tppcm(dat, "rank1")
  tp <- get_parts(x_r1)
  cx <- coef(x_r1); vx <- sqrt(diag(xxirt_vcov(x_r1)))
  ix <- c(grep("alpha", names(cx)), match(c("I1_g2", "I1_g3"), names(cx)))
  expect_equal(irt_pars(x_r1, restricted = TRUE)$se, unname(vx[ix]), tolerance = 1e-3)   # needs the curvature term
  expect_equal(tp$restricted$est, unname(cx[ix]), tolerance = 1e-3)
})

test_that("fixed steps in tam.mml.3pl are detected", {
  dat <- make_data()
  fit <- TAM::tam.mml.3pl(dat, E = step_design(dat), est.variance = FALSE, verbose = FALSE,
                          gammaslope.fixed = cbind(2, 1.5))
  tp <- get_parts(fit)
  expect_equal(tp$fixed, 2)
  expect_true(is.na(irt_pars(fit, long = TRUE)$se[2]))
  expect_equal(dim(tp$vcov), c(11 + 12, 11 + 12))
  expect_output(print(tp), "get_parts")
})

test_that("get_parts works as an accessor", {
  dat <- make_data()
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  x <- get_parts(g)
  expect_identical(get_parts(x), x)
  expect_equal(get_parts(g, "disc"), x$a)
  expect_equal(unname(get_parts(x, "se")), x$restricted$se)
  expect_equal(get_parts(x, "npar"), 16)
  expect_equal(as.numeric(logLik(x)), as.numeric(logLik(g)))
  expect_equal(attr(logLik(x), "df"), 16)
  expect_named(get_parts(x, c("loglik", "nobs")), c("loglik", "nobs"))
  expect_error(get_parts(x, "nonsense"))
})
