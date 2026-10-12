tucker_fixture <- function(missing_train = TRUE) {
  .vip_with_seed(123, function() {
    A <- matrix(rnorm(48), 24, 2)
    B <- matrix(rnorm(24), 12, 2)
    C <- matrix(rnorm(8), 4, 2)
    G <- array(rnorm(8), c(2, 2, 2))
    X <- .tk_reconstruct(G, A, B, C)
    dimnames(X) <- list(paste0('S', 1:24), paste0('G', 1:12),
      paste0('T', 1:4))
    tr <- X[1:18, , , drop = FALSE]
    te <- X[19:24, , , drop = FALSE]
    if (missing_train)
      tr[sample.int(length(tr), round(.1 * length(tr)))] <- NA_real_
    te[sample.int(length(te), round(.2 * length(te)))] <- NA_real_
    list(tr = tr, te = te, truth = X, Y = matrix(rep(0:1, 12), 24,
      1, dimnames = list(dimnames(X)[[1]], 'Case')))
  })
}
test_that('3D train/apply Tucker recovers low rank and preserves every observation', {
  d <- tucker_fixture(FALSE)
  f <- fit_tucker_imputer(d$tr, c(2, 2, 2), tolerance = 1e-8)
  expect_true(f$convergence$converged)
  p <- predict_tucker_imputer(f, d$te, diagnostics = TRUE)
  expect_identical(dim(p$X), dim(d$te))
  expect_identical(dimnames(p$X), dimnames(d$te))
  expect_identical(p$X[!is.na(d$te)], d$te[!is.na(d$te)])
  expect_identical(f$X_train_imputed, d$tr)
  miss <- is.na(d$te)
  expect_lt(sqrt(mean((p$X[miss] - d$truth[19:24, , , drop = FALSE][miss])^2)),
    1e-6)
  # Verify the stated Kronecker formula independently against tensor contractions.
  rec <- .tk_reconstruct(f$G, f$A_train, f$B, f$C)
  expect_equal(unname(matrix(rec, nrow(d$tr))), unname(f$A_train %*% t(f$basis)), tolerance = 1e-12)
  d <- tucker_fixture(TRUE)
  f <- fit_tucker_imputer(d$tr, c(2, 2, 2), tolerance = 1e-7)
  p <- predict_tucker_imputer(f, d$te)
  expect_true(f$convergence$converged)
  expect_identical(f$X_train_imputed[!is.na(d$tr)], d$tr[!is.na(d$tr)])
  expect_identical(p[!is.na(d$te)], d$te[!is.na(d$te)])
  miss <- is.na(d$te)
  expect_lt(sqrt(mean((p[miss] - d$truth[19:24, , ,
    drop = FALSE][miss])^2)) / sd(as.vector(d$truth)),
    .15)
  expect_equal(f$preprocessing$center, apply(d$tr, c(2, 3), mean, na.rm = TRUE))
})
test_that('projection uses frozen factors and treats unseen subjects independently', {
  d <- tucker_fixture()
  f <- fit_tucker_imputer(d$tr, c(2, 2, 2))
  original <- f
  local_mocked_bindings(fit_tucker_imputer = function(...)
    stop('unexpected Tucker refit'), .package = 'TensorPLS')
  a <- predict_tucker_imputer(f, d$te[1, , , drop = FALSE])
  b <- predict_tucker_imputer(f, d$te)
  expect_identical(a, b[1, , , drop = FALSE])
  expect_identical(f, original)
  extreme <- d$te
  extreme[-1, , ] <- extreme[-1, , ] * 10000
  expect_identical(a, predict_tucker_imputer(f, extreme)[1, , , drop = FALSE])
  expect_identical(rownames(f$B), dimnames(d$tr)[[2]])
  expect_identical(f$training_ids, dimnames(d$tr)[[1]])
  expect_length(intersect(f$training_ids, dimnames(d$te)[[1]]), 0)
  empty <- d$te[1, , , drop = FALSE]
  empty[] <- NA_real_
  expect_error(predict_tucker_imputer(f, empty), 'no observed entries')
  one <- empty
  one[1, 1, 1] <- 2
  expect_warning(p <- predict_tucker_imputer(f, one, TRUE), 'rank deficient')
  expect_true(p$diagnostics$rank_deficient)
  expect_identical(p$X[1, 1, 1], 2)
})
test_that('initialization is deterministic and nonconvergence cannot be silently used', {
  d <- tucker_fixture()
  set.seed(34)
  rng <- .Random.seed
  a <- fit_tucker_imputer(d$tr, c(2, 2, 2), seed = 123)
  expect_identical(.Random.seed, rng)
  expect_identical(a, fit_tucker_imputer(d$tr, c(2, 2, 2), seed = 123))
  expect_warning(b <- fit_tucker_imputer(d$tr, c(2, 2, 2), max_iter = 1),
    'did not converge')
  expect_false(b$convergence$converged)
  expect_error(predict_tucker_imputer(b, d$te), 'non-converged')
  expect_error(fit_tucker_imputer(d$tr, c(20, 2, 2)), 'ranks')
  bad <- d$te[, , 4:1, drop = FALSE]
  expect_error(predict_tucker_imputer(a, bad), 'ordered feature/time')
})
test_that('Tucker nested folds isolate X and Y, keep 3D boundaries and choose ncomp internally', {
  d <- tucker_fixture()
  X <- d$truth
  set.seed(21)
  X[sample.int(length(X), round(.1 * length(X)))] <- NA_real_
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, dimnames(X)[[1]])
  run <- function(x = X, y = d$Y)
    nested_cv_tensorpls(x, y, folds = folds, ncomp = 1:2,
    vip_thr = 50, imputation = 'tucker', tucker_control = list(ranks = c(2,
      2, 2), tolerance = 1e-5), input_provenance = 'raw', verbose = FALSE)
  real <- fit_tucker_imputer
  seen <- list()
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ...) {
    expect_length(dim(X_train), 3L)
    seen[[length(seen) + 1]] <<- dimnames(X_train)[[1]]
    real(X_train, ...)
  }, .package = 'TensorPLS')
  a <- run()
  expect_true(all(a$outer_fold_metrics$status == 'ok'))
  # Three inner training fits plus one outer training fit per outer fold, cached across ncomp.
  expect_length(seen, 12L)
  for (f in 1:3) {
    expected <- names(folds$outer)[folds$outer != f]
    expect_identical(a$models[[f]]$imputer$training_ids, expected)
    expect_identical(a$models[[f]]$imputer$B, real(X[expected, , , drop = FALSE],
      ranks = c(2, 2, 2), tolerance = 1e-5)$B)
    for (z in a$inner_tuning_results[[f]]$imputation_diagnostics)
      expect_true(all(z$training_ids %in%
      expected))
    t <- a$inner_tuning_results[[f]]$tuning
    expect_equal(unname(a$selected_ncomp_per_outer_fold[f]), t$ncomp[which.max(t$selection_score)])
  }
  xt <- X
  xt[folds$outer == 1, , ] <- xt[folds$outer == 1, , ] * 100 + 20
  b <- run(xt)
  expect_identical(a$models[[1]], b$models[[1]])
  expect_identical(a$inner_tuning_results[[1]], b$inner_tuning_results[[1]])
  yt <- d$Y
  yt[folds$outer == 1, 1] <- 1 - yt[folds$outer == 1, 1]
  c <- run(y = yt)
  expect_identical(a$models[[1]], c$models[[1]])
  expect_identical(a$inner_tuning_results[[1]], c$inner_tuning_results[[1]])
  expect_identical(a$selected_ncomp_per_outer_fold[1], c$selected_ncomp_per_outer_fold[1])
  expect_identical(subset(a$outer_predictions, outer_fold == 1)$Y_pred, subset(c$outer_predictions,
    outer_fold == 1)$Y_pred)
  expect_identical(a, run())
  model <- a$models[[1]]
  test <- X[folds$outer == 1, , , drop = FALSE]
  expect_equal(.nested_predict_pipeline(model, test[1, , , drop = FALSE]),
    .nested_predict_pipeline(model,
    test)[1, , drop = FALSE], ignore_attr = TRUE)
})
test_that('QR-compressed core solve equals the explicit observed-entry design', {
  d <- tucker_fixture()
  f <- fit_tucker_imputer(d$tr, c(2, 2, 2))
  z <- sweep(d$tr, c(2, 3), f$preprocessing$center, '-')
  obs <- !is.na(z)
  z[!obs] <- 0
  index <- which(obs, arr.ind = TRUE)
  D <- t(vapply(seq_len(nrow(index)), function(i) {
    q <- index[i, ]
    as.vector(kronecker(f$C[q[3], ], kronecker(f$B[q[2], ], f$A_train[q[1],
      ])))
  }, numeric(8)))
  direct <- array(.tk_solve(D, matrix(z[obs], ncol = 1), 1e-12), c(2,
    2, 2))
  compressed <- .tk_weighted_core(z, obs, list(f$A_train, f$B, f$C), c(2, 2,
    2))
  expect_equal(compressed, direct, tolerance = 1e-10)
  expect_true(all(diff(f$history$objective) < 1e-9))
  r <- fit_tucker_imputer(d$tr, c(2, 2, 2), seed = 123, initialization = 'random')
  expect_identical(r, fit_tucker_imputer(d$tr, c(2, 2, 2), seed = 123,
    initialization = 'random'))
})
test_that('compressed time solve equals its dense observed design', {
  d <- tucker_fixture()
  f <- fit_tucker_imputer(d$tr, c(2, 2, 2))
  z <- sweep(d$tr, c(2, 3), f$preprocessing$center, '-')
  obs <- !is.na(z)
  z[!obs] <- 0
  Gu <- matrix(aperm(f$G, c(3, 1, 2)), 2)
  D <- kronecker(f$B, f$A_train) %*% t(Gu)
  patterns <- .tk_patterns(obs, 3)
  dense <- .tk_observed_coordinates(aperm(z, c(3, 1, 2)), D, patterns$patterns,
    patterns$groups, 1e-10, sqrt(f$ridge) * t(Gu))
  compressed <- .tk_time_coordinates(z, obs, list(f$A_train, f$B, f$C), f$G, f$ridge,
    1e-10)
  expect_equal(compressed, dense, tolerance = 1e-10)
})
