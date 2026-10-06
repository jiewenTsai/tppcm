dsf_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  a <- matrix(c(1, 1.5, 1.2), 4, 3, byrow = TRUE)
  b <- matrix(c(-1, 0, 1), 4, 3, byrow = TRUE)
  a2 <- a; a2[2, 2] <- 2.5
  list(dat = rbind(sim_tppcm(n, disc = a, diff = b),
                   sim_tppcm(n, disc = a2, diff = b, theta = stats::rnorm(n, 0.3, 0.8))),
       g = rep(c("a", "b"), each = n))
}

fit_mg <- function(d, g) TAM::tam.mml.3pl(d$resp, A = d$A, E = d$E, group = g,
                                          variance.fixed = d$variance.fixed,
                                          beta.fixed = d$beta.fixed, verbose = FALSE)
ll <- function(m) as.numeric(logLik(m))

test_that("scalar design equals the multigroup TPPCM with est.variance = FALSE", {
  s <- dsf_data()
  d0 <- dsf_design(s$dat, group = s$g)
  expect_s3_class(d0, "tppcm_dsf_design")
  expect_equal(dim(d0$resp), c(800, 8))
  expect_equal(dim(d0$E)[4], 12)
  expect_equal(dim(d0$A)[3], 12)
  expect_true(all(is.na(d0$resp[s$g == "a", grepl("_b$", names(d0$resp))])))
  expect_equal(d0$variance.fixed[, 1], 1); expect_equal(d0$beta.fixed[, 1], 1)
  m0 <- fit_mg(d0, s$g)
  mg <- TAM::tam.mml.3pl(s$dat, E = tppcm(s$dat), group = s$g, est.variance = FALSE, verbose = FALSE)
  expect_equal(ll(m0), ll(mg), tolerance = 1e-3)
  expect_true(abs(m0$variance[2] - 1) > 1e-3)
  expect_output(print(d0), "scalar")
})

test_that("configural equals separate calibrations; metric and scalar are nested", {
  s <- dsf_data(300)
  dc <- dsf_design(s$dat, s$g, group.equal = NULL)
  expect_equal(dc$variance.fixed[, 1], 1:2); expect_equal(dc$beta.fixed[, 1], 1:2)
  expect_equal(dim(dc$E)[4], 24); expect_equal(dim(dc$A)[3], 24)
  mc <- fit_mg(dc, s$g)
  sep <- sum(sapply(c("a", "b"), function(k) {
    xx <- s$dat[s$g == k, ]
    ll(TAM::tam.mml.3pl(xx, E = tppcm(xx, K = 4), est.variance = FALSE, verbose = FALSE))
  }))
  expect_equal(ll(mc), sep, tolerance = 1e-2)
  dm <- dsf_design(s$dat, s$g, group.equal = "cbd")
  expect_equal(dm$variance.fixed[, 1], 1); expect_equal(dm$beta.fixed[, 1], 1:2)
  mm <- fit_mg(dm, s$g); ms <- fit_mg(dsf_design(s$dat, s$g), s$g)
  expect_gte(ll(mc), ll(mm) - 1e-3); expect_gte(ll(mm), ll(ms) - 1e-3)
  dg <- dsf_design(s$dat, s$g, group.equal = NULL, design = ~ item)   # configural GPCM
  expect_equal(dim(dg$E)[4], 8)
})

test_that("group.partial releases single steps; cbl releases the IBD with the CBD", {
  s <- dsf_data()
  d1 <- dsf_design(s$dat, s$g, group.partial = list(cbd = "I2:2", ibd = "I2:2"))
  expect_equal(dim(d1$E)[4], 13); expect_equal(dim(d1$A)[3], 13)
  expect_equal(dimnames(d1$E)[[4]][13], "I2_cbd2_b")
  expect_equal(unname(d1$index["I2_b", 2]), 13L)
  expect_identical(dsf_design(s$dat, s$g, group.partial = "I2:2")$index, d1$index)
  dl <- dsf_design(s$dat, s$g, group.equal = c("cbd", "cbl"), group.partial = list(cbd = "I2:2"))
  expect_true(dl$partial$ibd["I2", 2])
  dm <- fit_mg(dsf_design(s$dat, s$g), s$g)
  expect_gt(2 * (ll(fit_mg(d1, s$g)) - ll(dm)), 3.84)          # simulated DSF in I2 step 2
  expect_error(dsf_design(s$dat, s$g, group.equal = "cbl"), "not a linear restriction")
  expect_error(dsf_design(s$dat, s$g, group.equal = "loadings"), "unknown parameter")
  expect_error(dsf_design(s$dat, s$g, group.partial = list(foo = "I1")), "cbd and/or ibd")
  expect_error(dsf_design(s$dat, s$g, group.partial = "I9:1"), "unknown item")
  expect_error(dsf_design(s$dat, s$g, group.partial = "I1:4"), "has no step")
  expect_error(dsf_design(s$dat, group = rep(1, 800)), "at least two")
})

test_that("items with fewer categories get NA rows in A", {
  s <- dsf_data(150)
  dat <- s$dat; dat[, 4] <- pmin(dat[, 4], 2)
  d <- dsf_design(dat, group = s$g)
  expect_true(all(is.na(d$A["I4_a", "Category3", ])))
  mg <- TAM::tam.mml.3pl(dat, E = tppcm(dat), group = s$g, est.variance = FALSE, verbose = FALSE)
  expect_equal(ll(fit_mg(d, s$g)), ll(mg), tolerance = 1e-3)
})

test_that("dsf() with the scalar baseline: location and location-free CBD tests", {
  s <- dsf_data()
  fit <- TAM::tam.mml.3pl(s$dat, E = tppcm(s$dat), group = s$g, est.variance = FALSE, verbose = FALSE)
  r <- dsf(fit, items = "I2", verbose = FALSE)
  expect_s3_class(r, "tppcm_table")
  expect_named(r, c("item", "step", "par", "a", "b", "LR", "df", "p", "p_adj"))
  expect_equal(nrow(r), 6); expect_equal(r$par, rep(c("cbd", "ibd"), 3))
  cb <- r[r$par == "cbd", ]
  expect_lt(cb$p[2], .05); expect_gt(cb$b[2], cb$a[2])          # simulated DSF in I2 step 2
  # CBD row = (CBD and IBD released) vs (IBD released); joint = sum of the two rows
  ctl <- fit$control; ctl$conv <- 1e-6; ctl$convD <- 1e-5; ctl$maxiter <- 3000; ctl$progress <- FALSE
  f <- function(d) TAM::tam.mml.3pl(d$resp, A = d$A, E = d$E, group = s$g, variance.fixed = d$variance.fixed,
                                    beta.fixed = d$beta.fixed, control = ctl, verbose = FALSE)
  joint <- 2 * (ll(f(dsf_design(s$dat, s$g, group.partial = "I2:2"))) - ll(f(dsf_design(s$dat, s$g))))
  expect_equal(sum(r$LR[r$step == 2]), joint, tolerance = 1e-3)
  rl <- dsf(fit, group.equal = c("cbd", "cbl"), items = "I2", steps = 1, verbose = FALSE)
  expect_equal(rl$par, c("cbd", "cbl")); expect_equal(rl$LR, r$LR[1:2], tolerance = 1e-6)
  imp <- attr(r, "impact")
  expect_equal(imp$VAR[1], 1); expect_equal(imp$MEAN[1], 0)
})

test_that("dsf() with the metric baseline, partial baselines and bad input", {
  s <- dsf_data(200)
  fit <- TAM::tam.mml.3pl(s$dat, E = tppcm(s$dat), group = s$g, est.variance = FALSE, verbose = FALSE)
  rm <- dsf(fit, group.equal = "cbd", items = "I1", verbose = FALSE)
  expect_equal(rm$par, rep("cbd", 3))
  expect_match(attr(rm, "title"), "metric")
  rp <- dsf(fit, group.partial = "I1:1", items = "I1", verbose = FALSE)
  expect_false(any(rp$step == 1))                               # released steps are not tested
  fg <- TAM::tam.mml.3pl(s$dat, E = tppcm(s$dat, design = ~ item), group = s$g, est.variance = FALSE, verbose = FALSE)
  expect_equal(nrow(dsf(fg, items = "I1", steps = 1, verbose = FALSE)), 2)
  f1 <- TAM::tam.mml.3pl(s$dat, E = tppcm(s$dat), est.variance = FALSE, verbose = FALSE)
  expect_error(dsf(f1), "single group")
  expect_error(dsf(TAM::tam.mml(s$dat, verbose = FALSE)), "tam.mml.3pl")
  expect_error(dsf(fg, items = "I9"), "unknown item")
  r2 <- dsf(fg, items = c("I2:1", "I1:3"), verbose = FALSE)   # single steps, as from dif_test(by = "step")
  expect_equal(paste0(r2$item, ":", r2$step), c("I1:3", "I1:3", "I2:1", "I2:1"))
  expect_equal(r2$p_adj, p.adjust(r2$p, "holm"))
  expect_error(dsf(fg, group.equal = "ibd"), "must contain")
})
