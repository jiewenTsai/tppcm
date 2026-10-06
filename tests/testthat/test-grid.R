test_that("every function runs on every supported fit without warnings", {
  skip_on_cran()
  dat <- as.matrix(make_data(600))
  q <- function(x) suppressWarnings(suppressMessages(x))
  fits <- list(
    pcm = q(TAM::tam.mml(dat, verbose = FALSE)),
    rsm = q(TAM::tam.mml(dat, irtmodel = "RSM", verbose = FALSE)),
    gpcm = q(TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)),
    nrm = q(TAM::tam.mml.2pl(dat, irtmodel = "2PL", verbose = FALSE)),
    tppcm = q(TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)),
    step = q(TAM::tam.mml.3pl(dat, E = tppcm(dat, design = ~ step), est.variance = FALSE, verbose = FALSE)),
    xx_gpcm = q(xxirt_tppcm(dat, design = ~ item)))
  g <- factor(rep(1:2, 300))
  saturated <- c("nrm", "tppcm")
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_no_warning(x <- get_parts(f))
    expect_no_warning(capture.output(print(irt_pars(x))))
    expect_no_warning(irt_pars(x, long = TRUE))
    expect_no_warning(irt_pars(x, restricted = TRUE))
    expect_no_warning(step_info(x))
    expect_no_warning(m2(x))
    expect_no_warning(suppressMessages(wald_test(x)))
    expect_no_warning(suppressMessages(item_test(x)))
    expect_no_warning(suppressMessages(mi(x)))
    if (nm %in% saturated) expect_error(score_test(x), "nothing to test") else expect_no_warning(score_test(x))
    if (nm == "rsm") expect_error(dif_test(x, g), "free step locations") else expect_no_warning(suppressMessages(dif_test(x, g)))
  }
})

test_that("unsupported or special TAM fits are handled", {
  skip_on_cran()
  set.seed(9)
  d <- as.matrix(make_data(600))
  w <- sample(1:3, 600, TRUE)
  gw <- TAM::tam.mml.2pl(d, irtmodel = "GPCM", pweights = w, verbose = FALSE)
  ge <- TAM::tam.mml.2pl(d[rep(1:600, w), ], irtmodel = "GPCM", verbose = FALSE)
  expect_equal(irt_pars(gw, fit = FALSE)$se, irt_pars(ge, fit = FALSE)$se, tolerance = 1e-5)
  expect_error(m2(gw), "weights")
  th <- rnorm(1500)
  bin <- sapply(1:8, function(i) rbinom(1500, 1, 0.2 + 0.8 * plogis(1.3 * (th - seq(-1, 1, len = 8)[i]))))
  f3 <- TAM::tam.mml.3pl(bin, guess = rep(.2, 8), est.variance = FALSE, verbose = FALSE)
  expect_error(get_parts(f3), "divide-by-total")
  r <- TAM::tam.mml(d, Y = cbind(x = rnorm(600)), verbose = FALSE)
  expect_warning(x <- get_parts(r), "latent regression")
  expect_error(m2(x), "latent regression")
  mf <- TAM::tam.mml.mfr(d, facets = data.frame(r = rep(1:2, 300)), formulaA = ~ item + item:step + r,
                         verbose = FALSE)
  expect_error(get_parts(mf), "many-facet")
  g2 <- TAM::tam.mml.2pl(d, irtmodel = "GPCM", group = rep(1:2, 300), verbose = FALSE)
  expect_no_warning(get_parts(g2))
})
