# Tests added after the code review (October 2026)

test_that("multigroup fits carry the free group parameters", {
  skip_on_cran()
  data(data.timssAusTwn.scored, package = "TAM")
  dt <- data.timssAusTwn.scored[, 1:11]; grp <- data.timssAusTwn.scored$IDCNTRY
  f1 <- TAM::tam.mml.2pl(dt, irtmodel = "GPCM", group = grp, verbose = FALSE)
  x <- get_parts(f1)
  expect_equal(x$gnames, c("grp:mean2", "grp:var2"))
  expect_equal(tppcm:::.npar(x), f1$ic$Npars)
  g <- get_parts(x, "groups")
  expect_true(all(is.finite(g$SE_MEAN[2]), is.finite(g$SE_VAR[2])))
  # group block of the Louis information equals the numerical second derivative
  # of the marginal log-likelihood in the group parameters (item parameters fixed)
  lev <- sort(unique(grp)); p <- nrow(x$parindex)
  ll <- function(m2, v2) {
    tot <- 0
    for (gg in 1:2) {
      sel <- grp == lev[gg]
      w <- stats::dnorm(x$theta, c(0, m2)[gg], sqrt(c(1, v2)[gg])); w <- w / sum(w)
      lik <- rep(1, sum(sel))
      lik <- matrix(1, sum(sel), length(x$theta))
      for (i in seq_len(ncol(dt))) {
        L <- x$ncat[i] - 1; a <- x$a[i, 1:L]; d <- a * x$diff[i, 1:L]
        P <- tppcm:::.tppcm_P(a, d, matrix(x$theta))
        r <- !is.na(dt[sel, i]); lik[r, ] <- lik[r, ] * t(P[, dt[sel, i][r] + 1, drop = FALSE])
      }
      tot <- tot + sum(log(lik %*% w))
    }
    tot
  }
  h <- 1e-3; m <- g$MEAN[2]; v <- g$VAR[2]
  num_mm <- -(ll(m + h, v) - 2 * ll(m, v) + ll(m - h, v)) / h^2
  num_vv <- -(ll(m, v + h) - 2 * ll(m, v) + ll(m, v - h)) / h^2
  expect_equal(unname(x$info[p + 1, p + 1]), num_mm, tolerance = 1e-3)
  expect_equal(unname(x$info[p + 2, p + 2]), num_vv, tolerance = 1e-3)
  # tests and DIF run on the multigroup fit
  expect_s3_class(score_test(x), "tppcm_test")
  expect_true(nrow(mi(x)) > 0)
  expect_equal(as.character(dif_test(x, factor(grp), parm = "diff", by = "item")$unit), colnames(dt))
})

test_that("score_test checks containment and category numbers", {
  dat <- make_data()
  g <- TAM::tam.mml.3pl(dat, E = tppcm(dat, design = ~ item), est.variance = FALSE, verbose = FALSE)
  expect_error(score_test(g, against = ~ step), "does not contain")
  d2 <- as.matrix(dat); d2[, 1] <- pmin(d2[, 1], 2)
  g2 <- TAM::tam.mml.3pl(d2, E = tppcm(d2, design = ~ item), est.variance = FALSE, verbose = FALSE)
  expect_error(score_test(g2, against = ~ item + step), "same number of categories")
  expect_output(print(score_test(g)), "Louis information")
  expect_match(tppcm:::.label(get_parts(g)), "GPCM, design ~ item (tam.mml.3pl)", fixed = TRUE)
})

test_that("irt_pars flags near-zero discriminations and validates level", {
  data(data.gpcm, package = "TAM")
  f <- TAM::tam.mml.3pl(data.gpcm, E = tppcm(data.gpcm), est.variance = FALSE, verbose = FALSE)
  p <- irt_pars(f, fit = FALSE)
  expect_equal(p$near_zero, "Work_disc3")
  expect_true(is.na(p$items["Work", "diff3"]))
  expect_false(is.na(irt_pars(f, fit = FALSE, IRTpars = FALSE)$items["Work", "int3"]))
  expect_error(irt_pars(f, level = 95), "between 0 and 1")
  expect_s3_class(as.data.frame(p), "data.frame")
  expect_equal(names(coef(get_parts(f)))[1], "Comfort_disc1")
})

test_that("xxirt_tppcm objects print a summary", {
  dat <- make_data()
  fit <- xxirt_tppcm(as.matrix(dat), design = ~ item)
  expect_s3_class(fit, "xxirt_tppcm")
  expect_output(print(fit), "xxirt fit of the GPCM (design ~ item)", fixed = TRUE)
  expect_equal(get_parts(fit, "model"), "gpcm")
})

test_that("tppcmtree keeps rows with missing responses and supports impact", {
  skip_on_cran()
  set.seed(5)
  dat <- as.matrix(make_data(N = 800))
  dat[sample(length(dat), 40)] <- NA
  df <- data.frame(g = factor(rep(1:2, 400)), z = rnorm(800))
  df$resp <- dat
  tr <- tppcmtree(resp ~ g + z, data = df, parm = "diff", minsize = 200)
  expect_equal(sum(sapply(partykit::nodeids(tr, terminal = TRUE),
                          function(i) tr[[i]]$node$info$nobs)), 800)
  expect_output(print(tr), "terminal node")
  tr2 <- tppcmtree(resp ~ g + z, data = df, parm = "diff", minsize = 200, impact = "g")
  cf <- coef(tr2); nm <- if (is.matrix(cf)) colnames(cf) else names(cf)
  expect_true(any(grepl("^grp:mean", nm)))
  tr3 <- tppcmtree(resp ~ g, data = df, design = ~ step, parm = "diff", minsize = 200)
  expect_output(print(tr3), "node design: ~ step", fixed = TRUE)
  expect_error(tppcmtree(resp ~ g, data = df, design = ~ item + step), "Rank-1")
  expect_error(tppcmtree(resp ~ g, data = df, model = "pcm"), "no 'model' argument")
  df2 <- df; df2$g[1] <- NA
  expect_error(tppcmtree(resp ~ g + z, data = df2), "missing values")
  expect_error(tppcmtree(z ~ g, data = df), "matrix column")
})

test_that("dif_test reminds about impact on single-group fits", {
  dat <- make_data()
  f <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
  expect_message(dif_test(f, factor(rep(1:2, 300)), parm = "diff", by = "item"), "impact")
  expect_silent(dif_test(f, rnorm(600), parm = "diff", by = "block"))
})
