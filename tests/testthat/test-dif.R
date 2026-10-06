test_that("dif_test LM equals the closed form for a binary covariate", {
  skip_if_not_installed("strucchange")
  dat <- make_data()
  fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
  info <- get_parts(fit)
  set.seed(3); g <- factor(sample(1:2, nrow(dat), TRUE))
  res <- suppressMessages(dif_test(info, g))
  S <- get_parts(info, "estfun"); S <- sweep(S, 2, colMeans(S))
  s <- tppcm:::.efficient(S, 5)[, 1]
  n1 <- sum(g == 1); n2 <- sum(g == 2)
  stat <- sum(s[g == 1])^2 * (1 / n1 + 1 / n2) / mean(s^2)
  expect_equal(res$statistic[5], stat, tolerance = 1e-8)
  expect_equal(res$unit[5], "I2_disc2")
  expect_equal(nrow(suppressMessages(dif_test(info, g, by = "item"))), 4)
  expect_equal(suppressMessages(dif_test(info, g, parm = "diff", by = "block"))$npar, 12)
  expect_error(dif_test(info, g[-1]), "one value per person")
})

test_that("dif_test finds a single step DIF and not location DIF", {
  skip_if_not_installed("strucchange")
  set.seed(7)
  a <- matrix(c(1, 1.5, 0.7), 4, 3, byrow = TRUE); b <- matrix(c(-1, 0, 1), 4, 3, byrow = TRUE)
  a2 <- a; a2[3, 2] <- 0.3
  dat <- rbind(sim_tppcm(900, a, b), sim_tppcm(900, a2, b))
  g <- factor(rep(1:2, each = 900))
  fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
  res <- suppressMessages(dif_test(fit, g))
  expect_equal(res$unit[which.min(res$p)], "I3_disc2")
  expect_lt(min(res$p_adj), 0.05)
  expect_gt(min(suppressMessages(dif_test(fit, g, parm = "diff", by = "item"))$p_adj), 0.05)
})

test_that(".resid_span equals the normal-equations projection at full rank", {
  set.seed(1)
  X <- matrix(rnorm(200 * 4), 200); Y <- matrix(rnorm(200 * 2), 200)
  expect_equal(tppcm:::.resid_span(Y, X), Y - X %*% solve(crossprod(X), crossprod(X, Y)),
               tolerance = 1e-10)
  # a zero column and a duplicated column do not change the column space
  X2 <- cbind(X, 0, X[, 1] * 3)
  expect_equal(tppcm:::.resid_span(Y, X2), tppcm:::.resid_span(Y, X), tolerance = 1e-10)
})

test_that("dif_test works when a step discrimination is estimated at 0", {
  skip_if_not_installed("strucchange")
  skip_if_not_installed("sirt")
  set.seed(11)
  a <- matrix(c(1, 1.5, 1.2), 4, 3, byrow = TRUE); a[1, 1] <- -0.8
  b <- matrix(c(-1, 0, 1), 4, 3, byrow = TRUE)
  dat <- sim_tppcm(1000, a, b)
  fit <- suppressWarnings(xxirt_tppcm(dat, design = ~ item * step))        # bound a >= 0: I1_disc1 = 0
  expect_lt(abs(get_parts(fit, "disc")[1, 1]), 1e-6)
  g <- factor(rep(1:2, each = 500))
  expect_message(res <- dif_test(fit, g, parm = "all", by = "item"), "fewer testable directions in: I1")
  expect_equal(res$npar, c(5, 6, 6, 6))
  expect_true(all(is.finite(res$p)))
  res_p <- suppressMessages(dif_test(fit, g, parm = "all"))
  expect_true(is.na(res_p$p[res_p$unit == "I1_diff1"]))
  expect_true(all(is.finite(res_p$p[res_p$unit != "I1_diff1"])))
})

test_that("multigroup dif_test uses the split Louis score test", {
  skip_if_not_installed("strucchange")
  set.seed(5)
  a <- matrix(c(1, 1.8, 2), 4, 3, byrow = TRUE); b <- matrix(c(-1.5, 0, 1.5), 4, 3, byrow = TRUE)
  g <- factor(rep(1:2, each = 400))
  dat <- rbind(sim_tppcm(400, a, b), sim_tppcm(400, a, b, rnorm(400, 0.3)))
  fit <- suppressWarnings(TAM::tam.mml.3pl(dat, E = tppcm(dat), group = as.integer(g),
                                           est.variance = FALSE, verbose = FALSE))
  x <- get_parts(fit)
  # the information of the levels adds up to the total observed information
  L <- lapply(1:2, function(l) {
    gr <- x$grp; gr$g <- gr$g[g == l]
    tppcm:::.score_core(x$data[g == l, ], x$post[g == l, ], x$probs, x$theta, x$ncat, grp = gr)$info
  })
  expect_equal(L[[1]] + L[[2]], x$info, tolerance = 1e-8)
  res <- dif_test(x, g, parm = "all", by = "item")
  expect_match(capture.output(print(res))[1], "robust information")
  expect_match(capture.output(print(dif_test(x, g, parm = "I1_disc2", info = "louis")))[1], "louis information")
  expect_equal(res$npar[1:4], rep(6, 4))
  expect_false(any(startsWith(res$unit, "grp:")))           # group means/variances are not tested
  one <- dif_test(x, g, parm = "I1_disc2")
  expect_equal(one$npar, 1)
  expect_true(is.finite(one$p))
})

test_that("by = 'step' gives the joint 2-df test of a step, invariant to the parameterization", {
  skip_if_not_installed("strucchange")
  set.seed(11)
  a <- matrix(c(1, 1.5, 1.2), 4, 3, byrow = TRUE); b <- matrix(c(-1, 0, 1), 4, 3, byrow = TRUE)
  a2 <- a; a2[2, 2] <- 2.5
  dat <- rbind(sim_tppcm(400, a, b), sim_tppcm(400, a2, b, theta = stats::rnorm(400, .3, .8)))
  g <- factor(rep(c("a", "b"), each = 400))
  fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), group = g, est.variance = FALSE, verbose = FALSE)
  st <- dif_test(fit, g, by = "step")
  expect_equal(st$unit, paste0(rep(paste0("I", 1:4), each = 3), ":", 1:3))
  expect_equal(st$npar, rep(2, 12))
  expect_named(st, c("unit", "npar", "statistic", "df", "p", "p_adj"))
  expect_equal(st$df, rep(2L, 12))
  expect_equal(st$p_adj, p.adjust(st$p, "holm"))
  expect_equal(dif_test(fit, g, by = "step", adjust = "none")$p_adj, st$p)
  expect_match(attr(st, "title"), "p_adj: Holm")
  expect_error(dif_test(fit, g, adjust = "foo"), "'adjust' must be one of")
  for (inf in c("robust", "louis")) {
    irt <- dif_test(fit, g, parm = c("I2_disc2", "I2_diff2"), by = "block", info = inf)
    si  <- dif_test(fit, g, parm = c("I2_slope2", "I2_int2"), by = "block", param = "si", info = inf)
    expect_equal(irt$statistic, si$statistic, tolerance = 1e-6)
    if (inf == "robust") expect_equal(st$statistic[5], si$statistic, tolerance = 1e-8)
  }
  expect_identical(dif_test(fit, g), dif_test(fit, g, parm = "all"))
  expect_identical(dif_test(fit, g), dif_test(fit, g, parm = c("diff", "disc")))
  expect_false(any(startsWith(dif_test(fit, g)$unit, "grp:")))
  # all item parameters split: the group mean and variance drop out, without a message
  expect_silent(bl <- dif_test(fit, g, by = "block"))
  expect_equal(bl$npar, 22)
  expect_error(dif_test(fit, g, parm = "disc", by = "step"), "jointly")
  fg <- TAM::tam.mml.3pl(dat, E = tppcm(dat, design = ~ item), group = g, est.variance = FALSE, verbose = FALSE)
  expect_error(dif_test(fg, g, by = "step"), "one discrimination per step")
})
