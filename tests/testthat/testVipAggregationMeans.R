# Focused contracts for all-time and cumulative-component-path VIP means.
.vip_mean_fixture <- function() {
  set.seed(123)
  X <- array(rnorm(12 * 4 * 5), c(12, 4, 5), dimnames = list(paste0('s',
    1:12), paste0('f', 1:4), paste0('t', 1:5)))
  Y <- matrix(rep(c(0, 1), 6), ncol = 1)
  X[, 1, ] <- X[, 1, ] + Y[, 1]
  list(X = X, Y = Y)
}
test_that('model1 averages every time, including t4 and t5, for H=2 and H=3', {
  V <- array(seq_len(4 * 5 * 3) / 10, c(4, 5, 3))
  features <- paste0('f', 1:4)
  times <- paste0('t', 1:5)
  results <- list()
  for (H in 2:3) {
    a <- V[, , seq_len(H), drop = FALSE]
    fit <- list(VIP = matrix(a, 20, H), ncomp = H)
    z <- .nested_vips(fit, features, times)
    results[[as.character(H)]] <- z
    expect_equal(dim(z$VIP3Dmodel1), c(4L, H))
    expect_identical(colnames(z$VIP3Dmodel1), paste0('Comp', 1:H, '_meanTime'))
    for (f in 1:4)
      for (h in 1:H)
      expect_equal(unname(z$VIP3Dmodel1[f,
      h]), mean(a[f, 1:5, h]))
    # Each late time contributes 1/k of its change, independently.
    for (ti in 4:5) {
      changed <- a
      changed[, ti, ] <- changed[, ti, ] + 5
      zz <- .nested_vips(list(VIP = matrix(changed, 20, H), ncomp = H), features,
        times)
      expect_equal(unname(zz$VIP3Dmodel1 - z$VIP3Dmodel1), matrix(1, 4,
        H))
    }
    expect_false(isTRUE(all.equal(unname(z$VIP3Dmodel1), apply(a[, 1:H, , drop = FALSE],
      c(1, 3), mean))))
  }
  expect_identical(results[['2']]$VIP3Dmodel1, results[['3']]$VIP3Dmodel1[,
    1:2, drop = FALSE])
})
test_that('M2_mean averages all fitted cumulative VIP columns and retains feature/time mapping', {
  V <- array(seq_len(60) / 10, c(4, 5, 3))
  features <- paste0('f', 1:4)
  times <- paste0('t', 1:5)
  fit <- list(VIP = matrix(V, 20, 3), ncomp = 3L)
  z <- .nested_vips(fit, features, times)
  expect_equal(dim(z$M2_mean), c(4L, 5L))
  for (f in 1:4)
    for (ti in 1:5)
    expect_equal(unname(z$M2_mean[f,
    ti]), mean(V[f, ti, 1:3]))
  changed <- V
  changed[, , 3] <- changed[, , 3] + 3
  zz <- .nested_vips(list(VIP = matrix(changed, 20, 3), ncomp = 3L),
    features, times)
  expect_equal(unname(zz$M2_mean - z$M2_mean), matrix(1, 4, 5))
  expect_identical(unname(z$VIP2D), fit$VIP)
  expect_identical(rownames(z$VIP2D), as.vector(outer(features, times, paste, sep = '_')))
  expect_identical(dimnames(z$M2_mean), list(features, times))
})
test_that('public global path and nested helper agree with unchanged base VIP and unfolding', {
  d <- .vip_mean_fixture()
  X <- d$X
  Y <- d$Y
  Xmat <- matrix(X, 12)
  for (f in 1:4)
    for (ti in 1:5)
    expect_identical(Xmat[, f +
    4 * (ti - 1)], as.vector(X[, f, ti]))
  for (H in 2:3) {
    fit <- .tensorpls_fit(Xmat, Y, H)
    public <- tensorpls_vips(fit_tensorpls_pls(X, Y, ncomp = H))
    nested <- .nested_vips(fit, dimnames(X)[[2]], dimnames(X)[[3]])
    for (name in c("VIP2D", "VIP3Dmodel1"))
      expect_identical(public[[name]], nested[[name]])
    expect_identical(public$M2_mean, nested$M2_mean)
    expect_identical(unname(public$VIP2D), fit$VIP)
    r <- colMeans(cor(Y, fit$x.scores)^2)
    formula <- sapply(1:H, function(h)
      sqrt(20 * as.vector((fit$raw.wgs[,
      1:h, drop = FALSE]^2) %*% r[1:h]) / sum(r[1:h])))
    expect_equal(unname(public$VIP2D), formula, tolerance = 1e-12)
    V <- array(public$VIP2D, c(4, 5, H))
    expect_equal(unname(public$VIP3Dmodel1), apply(V, c(1, 3), mean))
    expect_equal(unname(public$M2_mean), apply(V, c(1, 2), mean))
    for (f in 1:4)
      for (ti in 1:5)
      expect_identical(as.vector(V[f,
      ti, ]), as.vector(public$VIP2D[f + 4 * (ti - 1), ]))
  }
})
