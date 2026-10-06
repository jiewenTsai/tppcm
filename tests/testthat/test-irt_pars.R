test_that("irt_pars gives TAM's GPCM locations and handles the PCM", {
  dat <- make_data()
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  p <- irt_pars(g)
  expect_s3_class(p, "irt_pars")
  expect_equal(dim(p$items), c(4, 6))
  tam_b <- g$item_irt$beta + as.matrix(g$item_irt[, grep("tau", names(g$item_irt))])
  expect_equal(unname(p$items[, 4:6]), unname(tam_b), tolerance = 1e-3)
  expect_output(print(p), "\\$items")

  pc <- irt_pars(TAM::tam.mml(dat, verbose = FALSE))
  expect_true(all(pc$items[, 1:3] == 1))
  expect_true(all(is.na(pc$se[, 1:3])))
  expect_false(is.na(pc$group$SE_VAR))

  x <- irt_pars(xxirt_tppcm(dat, design = ~ item + step), IRTpars = FALSE)
  expect_equal(c(x$group$MEAN, x$group$VAR), c(0, 1))
})

test_that("irt_pars labels disc/diff and slope/int", {
  dat <- make_data()
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
  expect_equal(colnames(irt_pars(g)$items), c(paste0("disc", 1:3), paste0("diff", 1:3)))
  p2 <- irt_pars(g, IRTpars = FALSE)
  expect_equal(colnames(p2$items), c(paste0("slope", 1:3), paste0("int", 1:3)))
  # cumulative step intercepts = mirt / psychotools category intercepts = -TAM AXsi_
  expect_equal(unname(t(apply(p2$items[, 4:6], 1, cumsum))), unname(-g$AXsi_[, -1]), tolerance = 1e-3)
})

test_that("restricted TAM locations (RSM), groups and multidimensional fits", {
  dat <- as.matrix(make_data())
  rsm <- TAM::tam.mml(dat, irtmodel = "RSM", verbose = FALSE)
  x <- get_parts(rsm)
  expect_equal(x$n_loc, ncol(dat) + 2)                 # item locations + 2 free thresholds
  expect_equal(ncol(x$J), 1 + ncol(dat) + 2)
  pcm <- get_parts(TAM::tam.mml(dat, verbose = FALSE))
  expect_true(all(irt_pars(x)$se[, 4:6] < irt_pars(pcm)$se[, 4:6]))
  expect_error(get_parts(x, "estfun"), "free step locations")
  g <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", group = rep(1:2, length.out = nrow(dat)), verbose = FALSE)
  expect_equal(nrow(irt_pars(g)$group), 2)
  m2 <- TAM::tam.mml.2pl(dat, irtmodel = "GPCM", Q = cbind(c(1, 1, 0, 0), c(0, 0, 1, 1)), verbose = FALSE)
  expect_error(get_parts(m2), "unidimensional")
})
