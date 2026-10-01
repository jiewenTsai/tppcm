test_that("step_info gives Fisher information split by steps", {
  dat <- as.matrix(make_data())
  fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
  inf <- step_info(fit, theta = seq(-3, 3, by = 0.5))
  expect_s3_class(inf, "step_info")
  expect_equal(inf$item, apply(inf$step, c(1, 2), sum))
  expect_true(all(inf$step >= -1e-12))
  sp <- get_parts(fit); th <- inf$theta; h <- 1e-5
  num <- sapply(seq_len(ncol(dat)), function(i) {
    P <- function(t) tppcm:::.tppcm_P(sp$a[i, ], sp$a[i, ] * sp$diff[i, ], matrix(t))
    dP <- (P(th + h) - P(th - h)) / (2 * h)
    rowSums(dP^2 / P(th))
  })
  expect_equal(unname(inf$item), unname(num), tolerance = 1e-6)
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  ic <- TAM::IRT.informationCurves(g)
  expect_equal(step_info(g, theta = ic$theta[, 1])$test, as.numeric(ic$test_info_curve), tolerance = 1e-3)
  expect_output(print(inf), "Fisher information")
})

