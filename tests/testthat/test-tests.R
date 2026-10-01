test_that("score tests are close to likelihood ratio tests", {
  dat <- make_data(1000)
  m_sat <- TAM::tam.mml.2pl(dat, irtmodel = "2PL", verbose = FALSE)
  x_r1  <- xxirt_tppcm(dat, "rank1")
  lr <- 2 * (as.numeric(logLik(m_sat)) - as.numeric(logLik(x_r1)))
  s1 <- score_test(x_r1)
  expect_equal(s1$df, 6)
  expect_lt(abs(s1$statistic - lr), 1)

  tp <- get_parts(TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE))
  expect_equal(score_test(tp)$df, 8)
  expect_equal(score_test(tp, "rank1")$df, 2)
  expect_equal(nrow(mi(tp)), 12)
  expect_true(all(item_test(tp)$df == 2))
})

test_that("Wald tests skip items with equal steps", {
  dat <- make_data()
  w <- wald_test(TAM::tam.mml.3pl(dat, E = step_design(dat), est.variance = FALSE, verbose = FALSE))
  expect_equal(nrow(w), 4); expect_true(all(w$df == 2))
  idx <- rbind(1:3, 4:6, rep(7, 3), rep(8, 3))
  w2 <- wald_test(TAM::tam.mml.3pl(dat, E = step_design(dat, index = idx), est.variance = FALSE,
                              verbose = FALSE))
  expect_equal(w2$item, c("I1", "I2"))
})
