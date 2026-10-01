test_that("estfun scores sum to the gradient and change coordinates correctly", {
  dat <- make_data()
  fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
  info <- get_parts(fit)
  ssi <- get_parts(info, "estfun.si"); sirt <- get_parts(info, "estfun")
  expect_equal(dim(ssi), c(nrow(dat), 24))
  g <- info$grad
  expect_equal(unname(colSums(ssi)), unname(c(g[info$parindex$type == "a"], -g[info$parindex$type == "d"])),
               tolerance = 1e-8)
  a <- info$par[info$parindex$type == "a"]
  expect_equal(unname(sirt[, 13:24]), unname(-ssi[, 13:24] * rep(a, each = nrow(dat))))
  expect_true(all(grepl("_diff", colnames(sirt)[13:24])))
  expect_true(all(grepl("_int", colnames(ssi)[13:24])))
  expect_equal(colnames(sirt)[1], "I1_disc1")
})

test_that("tppcmtree separates slope DIF from location DIF", {
  skip_if_not_installed("partykit")
  set.seed(42)
  a0 <- matrix(c(1, 1.4, 0.8), 5, 3, byrow = TRUE); b0 <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE) / a0
  N <- 2000; za <- rbinom(N, 1, .5); zb <- rbinom(N, 1, .5)
  resp <- matrix(NA_integer_, N, 5)
  for (ga in 0:1) for (gb in 0:1) {
    idx <- which(za == ga & zb == gb); a <- a0; b <- b0
    if (ga == 1) a[3, 2] <- 0.3
    if (gb == 1) b[5, ] <- b0[5, ] + 1
    resp[idx, ] <- as.matrix(sim_tppcm(length(idx), a, b))
  }
  colnames(resp) <- paste0("I", 1:5)
  df <- data.frame(za = factor(za), zb = factor(zb)); df$resp <- resp
  tr_a <- tppcmtree(resp ~ za + zb, data = df, parm = "disc", minsize = 300, maxdepth = 2)
  tr_b <- tppcmtree(resp ~ za + zb, data = df, parm = "diff", minsize = 300, maxdepth = 2)
  expect_equal(partykit::split_node(tr_a$node)$varid, 2L)   # za (variable 1 is resp)
  expect_equal(partykit::split_node(tr_b$node)$varid, 3L)   # zb
  expect_true("I3_disc2" %in% colnames(coef(tr_a)))
})

test_that("orthogonalized scores remove cross-block correlation", {
  dat <- make_data()
  fit <- TAM::tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
  S <- get_parts(fit, "estfun")
  So <- tppcm:::.orthogonalize(S, "slope", TRUE)
  bl <- attr(S, "block")
  expect_lt(max(abs(crossprod(So[, bl == "slope"], So[, bl == "location"]))), 1e-6)
  expect_identical(tppcm:::.orthogonalize(S, "all", TRUE), S)
})
