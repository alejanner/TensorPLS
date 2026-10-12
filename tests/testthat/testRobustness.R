source(system.file('examples', 'robustness_helpers.R', package = 'TensorPLS'), local = TRUE)
robust_fixture <- function()
  readRDS(system.file('extdata', 'vip_example.rds', package = 'TensorPLS'))
robust_small <- function(d, seed = 123, permutation = NULL)
  robustness_case(d$X, d$Y,
  fold_seed = seed, ncomp = 1:2, permutation = permutation, outer_folds = 3,
  inner_folds = 3, vip_thr = 50)
test_that('distinct repeat seeds generate distinct reproducible fold sets', {
  d <- robust_fixture()
  ids <- dimnames(d$X)[[1]]
  folds <- lapply(1:20, function(s)
    make_tensorpls_folds(d$Y, 3, 3,
    s, ids))
  signatures <- vapply(folds, function(f)
    paste(f$outer, collapse = ','), character(1))
  expect_equal(length(unique(signatures)), 20L)
  for (s in 1:20)
    expect_identical(folds[[s]], make_tensorpls_folds(d$Y,
    3, 3, s, ids))
  for (f in folds)
    for (j in 1:3) {
    expect_identical(names(f$inner[[j]]), ids[f$outer != j])
    expect_length(intersect(names(f$inner[[j]]), ids[f$outer == j]), 0L)
  }
})
test_that('permutation generation preserves RNG and only relabels fixed subjects', {
  d <- robust_fixture()
  Y <- d$Y
  X <- d$X
  set.seed(4)
  old <- .Random.seed
  permutations <- robustness_permutations(Y, 100, 555)
  expect_identical(.Random.seed, old)
  expect_identical(permutations, robustness_permutations(Y, 100, 555))
  expect_equal(length(unique(vapply(permutations, paste, character(1), collapse = ','))),
    100L)
  for (index in permutations) {
    changed <- robustness_response(Y, index)
    expect_identical(dimnames(changed), dimnames(Y))
    expect_identical(sort(as.vector(changed)), sort(as.vector(Y)))
    expect_identical(unname(changed), unname(Y[index, , drop = FALSE]))
    expect_identical(d$Y, Y)
    expect_identical(d$X, X)
  }
  expect_error(robustness_response(Y, rep(1, nrow(Y))))
})
test_that('full repeated procedure reproduces numerical results and selected genes', {
  d <- robust_fixture()
  a <- robust_small(d)
  b <- robust_small(d)
  expect_identical(a$results, b$results)
  expect_identical(a$folds, b$folds)
  expect_true(a$labels_restored)
  expect_true(b$labels_restored)
  expect_identical(a$results[[1]]$folds, a$results[[2]]$folds)
  for (r in a$results)
    for (f in 1:3) {
    t <- r$inner_tuning_results[[f]]$tuning
    i <- order(-t$selection_score, t$ncomp, seq_len(nrow(t)), na.last = TRUE)[1]
    expect_equal(r$selected_ncomp_per_outer_fold[[f]], t$ncomp[i])
  }
})
test_that('permuted pipeline receives unchanged X and permuted Y, restored afterward', {
  d <- robust_fixture()
  original <- d
  perm <- robustness_permutations(d$Y, 1, 78)[[1]]
  expected <- robustness_response(d$Y, perm)
  real_fit <- nested_cv_tensorpls
  calls <- 0L
  local_mocked_bindings(nested_cv_tensorpls = function(X, Y, ...) {
    calls <<- calls + 1L
    expect_identical(X, d$X)
    expect_identical(Y, expected)
    real_fit(X, Y, ...)
  }, .package = 'TensorPLS')
  a <- robust_small(d, permutation = perm)
  b <- robust_small(d, permutation = perm)
  expect_equal(calls, 4L)
  expect_identical(d, original)
  expect_identical(a$results, b$results)
  expect_true(a$labels_restored)
  for (r in a$results) {
    expect_identical(as.numeric(r$outer_predictions$Y_true), as.numeric(expected))
    expect_true(all(r$outer_fold_metrics$status == 'ok'))
  }
})
test_that('plot data contain one genuine outer prediction per subject', {
  d <- robust_fixture()
  a <- robust_small(d)
  ids <- dimnames(d$X)[[1]]
  for (r in a$results) {
    p <- robustness_predictions(r, ids)
    expect_identical(p$subject_id, ids)
    expect_equal(nrow(p), length(ids))
    expect_identical(p$Y_pred, r$outer_predictions$Y_pred)
    expect_identical(as.integer(p$outer_fold), as.integer(r$outer_fold_assignments[ids]))
    expect_false(any(c('LV1', 'LV2', 'Comp1', 'Comp2') %in% names(p)))
    broken <- r
    broken$outer_predictions$outer_fold[1] <- 999L
    expect_error(robustness_predictions(broken, ids))
  }
})
test_that('held-out labels cannot select components or affect fitting on frozen folds', {
  d <- robust_fixture()
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, dimnames(d$X)[[1]])
  fit <- function(Y)
    nested_cv_tensorpls(d$X, Y, folds = folds, ncomp = 1:3,
    vip_thr = 50, imputation = 'training_mean', input_provenance = 'raw', verbose = FALSE)
  a <- fit(d$Y)
  changed <- d$Y
  changed[folds$outer == 1, 1] <- 1 - changed[folds$outer == 1,
    1]
  b <- fit(changed)
  expect_identical(a$models[[1]], b$models[[1]])
  expect_identical(a$inner_tuning_results[[1]], b$inner_tuning_results[[1]])
  expect_identical(a$selected_ncomp_per_outer_fold[1], b$selected_ncomp_per_outer_fold[1])
  expect_identical(a$selected_features_per_outer_fold[[1]], b$selected_features_per_outer_fold[[1]])
  expect_identical(subset(a$outer_predictions, outer_fold == 1)$Y_pred, subset(b$outer_predictions,
    outer_fold == 1)$Y_pred)
  expect_false(identical(a$outer_fold_metrics$PRESS[1], b$outer_fold_metrics$PRESS[1]))
})
