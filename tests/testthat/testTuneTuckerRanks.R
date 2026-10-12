.rank_fixture <- function() {
  X <- array(seq_len(8 * 6 * 3) / 20, c(8, 6, 3))
  X[c(10, 40, 99)] <- NA_real_
  X
}
.rank_mock_fit <- function(X_train, ranks, ...) {
  z <- X_train
  z[is.na(z)] <- 0
  list(X_train_imputed = z, convergence = list(converged = TRUE, iterations = 2L,
    final_delta = 0, observed_RMSE = 123))
}
test_that('shared reproducible masks protect mode members and preserve caller state', {
  X <- .rank_fixture()
  before <- serialize(X, NULL)
  seen <- list()
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ranks, ...) {
    seen[[length(seen) + 1L]] <<- X_train
    .rank_mock_fit(X_train, ranks, ...)
  }, .package = 'TensorPLS')
  set.seed(52)
  rng <- .Random.seed
  ranks <- rbind(c(2, 2, 2), c(1, 1, 1))
  a <- tune_tucker_ranks(X, ranks, repeats = 3, verbose = FALSE)
  b <- tune_tucker_ranks(X, ranks, repeats = 3, verbose = FALSE)
  expect_identical(a, b)
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(X, NULL), before)
  expect_length(seen, 12)
  expect_false(identical(a$masks[[1]]$indices, a$masks[[2]]$indices))
  for (rep in 1:3) {
    idx <- a$masks[[rep]]$indices
    expect_true(all(!is.na(X[idx])))
    expect_length(intersect(idx, which(is.na(X))), 0)
    expect_identical(seen[[2 * rep - 1]], seen[[2 * rep]])
    for (j in 1:3)
      expect_true(all(apply(!is.na(seen[[2 * rep]]),
      j, sum) > 0))
    z <- a$repetition_results[a$repetition_results$repetition == rep, ]
    expect_equal(z$RMSE_masked, rep(sqrt(mean(X[idx]^2)), 2))
    expect_equal(z$MAE_masked, rep(mean(abs(X[idx])), 2))
    expect_true(all(z$observed_unchanged))
    expect_true(all(z$observed_RMSE == 123))
  }
  expect_identical(unname(a$selected_ranks), c(2L, 2L, 2L))
  expect_identical(a$tuning$rank_subject, c(2L, 1L))
  expect_identical(a$tuning$mean_RMSE[1], min(a$tuning$mean_RMSE))
  expect_equal(a$tuning$sd_RMSE[1],
    sd(a$repetition_results$RMSE_masked[a$repetition_results$candidate ==
    1]))
  rm('.Random.seed', envir = .GlobalEnv)
  invisible(tune_tucker_ranks(X, ranks, repeats = 1, verbose = FALSE))
  expect_false(exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE))
})
test_that('failed or altered fits remain visible and cannot be selected', {
  X <- .rank_fixture()
  count <- 0
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ranks, ...) {
    count <<- count + 1
    if (ranks[1] == 1 && count > 2)
      stop('fixture failure')
    z <- .rank_mock_fit(X_train, ranks, ...)
    if (ranks[1] == 2) {
      z$convergence$converged <- FALSE
      warning('fixture nonconvergence')
    }
    z
  }, .package = 'TensorPLS')
  a <- tune_tucker_ranks(X, rbind(c(1, 1, 1), c(2, 2,
    2)), repeats = 2, verbose = FALSE)
  expect_true(all(is.na(a$selected_ranks)))
  expect_false(any(a$tuning$valid))
  expect_identical(a$tuning$valid_repetitions, c(1L, 0L))
  expect_match(a$repetition_results$diagnostic[2], 'fixture nonconvergence')
  expect_match(a$repetition_results$diagnostic[3], 'fixture failure')
  expect_equal(nrow(a$repetition_results), 4)
  expect_match(a$diagnostics$message, 'No candidate')
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ranks, ...) {
    z <- .rank_mock_fit(X_train, ranks)
    z$X_train_imputed[!is.na(X_train)] <- 999
    z
  }, .package = 'TensorPLS')
  a <- tune_tucker_ranks(X, matrix(c(1, 1, 1), 1), repeats = 1,
    verbose = FALSE)
  expect_false(a$tuning$valid)
  expect_match(a$repetition_results$diagnostic, 'Observed non-masked entries changed')
})
test_that('candidate and masking validation matches Tucker constraints', {
  X <- .rank_fixture()
  bad <- list(matrix(numeric(), 0, 3), list(c(1, 1, 1)),
    matrix(1, 2, 2), matrix('1', 1, 3), rbind(c(1,
      1, 1), c(1, 1, 1)), matrix(c(1.5, 2,
      2), 1), matrix(c(0, 1, 1), 1), matrix(c(9,
      2, 2), 1), matrix(c(3, 1, 1), 1),
    matrix(c(NA, 1, 1), 1))
  for (r in bad)
    expect_error(tune_tucker_ranks(X, r, verbose = FALSE), 'rank|Rank|Duplicate')
  r <- matrix(c(1, 1, 1), 1)
  expect_error(tune_tucker_ranks(X, r, mask_fraction = 0), 'mask_fraction')
  expect_error(tune_tucker_ranks(X, r, tucker_control = list(ranks = c(1, 1,
    1))), 'tucker_control')
  expect_error(tune_tucker_ranks(X, r, tucker_control = list(bogus = 1)), 'tucker_control')
  X[1, , ] <- NA
  expect_error(tune_tucker_ranks(X, r), 'observed entry')
  expect_error(tune_tucker_ranks(array(1, c(1, 1, 1)), r), 'eligible')
})
test_that('real public Tucker fits preserve data and select minimum masked RMSE', {
  set.seed(18)
  X <- array(rnorm(10 * 7 * 3), c(10, 7, 3))
  X[c(5, 20)] <- NA_real_
  original <- fit_tucker_imputer
  fits <- list()
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ranks, tolerance = 1e-6,
    max_iter = 1000, ...) {
    fit <- original(X_train, ranks, tolerance = tolerance, max_iter = max_iter, ...)
    fits[[length(fits) + 1L]] <<- fit
    fit
  }, .package = 'TensorPLS')
  z <- tune_tucker_ranks(X, rbind(c(1, 1, 1), c(2, 2,
    2)), repeats = 2, tucker_control = list(tolerance = 1e-4, max_iter = 1000),
    verbose = FALSE)
  expect_length(fits, 4)
  for (i in 1:4) {
    idx <- z$masks[[z$repetition_results$repetition[i]]]$indices
    kept <- setdiff(which(!is.na(X)), idx)
    expect_identical(fits[[i]]$X_train_imputed[kept], X[kept])
    if (z$repetition_results$valid[i]) {
      error <- fits[[i]]$X_train_imputed[idx] - X[idx]
      expect_identical(z$repetition_results$RMSE_masked[i], sqrt(mean(error^2)))
      expect_identical(z$repetition_results$MAE_masked[i], mean(abs(error)))
    }
  }
  good <- which(z$tuning$valid)
  expect_true(length(good) > 0)
  best <- good[which.min(z$tuning$mean_RMSE[good])]
  expect_identical(unname(z$selected_ranks), as.integer(z$tuning[best, 1:3]))
})
