nested_test_data <- function() {
  .vip_with_seed(72, function() {
    X <- array(rnorm(36 * 14 * 3), c(36, 14, 3), dimnames = list(paste0("S",
      1:36), paste0("G", 1:14), paste0("T", 1:3)))
    y <- rep(0:1, 18)
    X[, 1, ] <- X[, 1, ] + 2 * y
    X[, 2, ] <- X[, 2, ] + y
    Y <- matrix(y, ncol = 1, dimnames = list(dimnames(X)[[1]], "Case"))
    list(X = X, Y = Y)
  })
}
nested_test_run <- function(d, folds = NULL, ...)
  nested_cv_tensorpls(d$X, d$Y, outer_folds = 3,
  inner_folds = 3, seed = 123, ncomp = 1:2, vip_thr = 65,
  input_provenance = "raw", verbose = FALSE, folds = folds, ...)
test_that("stratified folds preserve subject separation and RNG state", {
  d <- nested_test_data()
  set.seed(18)
  old <- .Random.seed
  ids <- dimnames(d$X)[[1]]
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, ids)
  expect_identical(.Random.seed, old)
  expect_identical(folds, make_tensorpls_folds(d$Y, 3, 3, 123, ids))
  expect_equal(length(folds$outer), 36L)
  expect_true(all(table(folds$outer) == 12))
  expect_true(all(table(folds$outer, d$Y[, 1]) == 6))
  for (f in 1:3) {
    tr <- ids[folds$outer != f]
    te <- ids[folds$outer == f]
    expect_length(intersect(tr, te), 0)
    expect_setequal(union(tr, te), ids)
    expect_identical(names(folds$inner[[f]]), tr)
    expect_length(intersect(names(folds$inner[[f]]), te), 0)
  }
  expect_error(make_tensorpls_folds(d$Y, 19, 3, 123, ids), "Each class")
})
test_that("changing OUTER TEST labels changes scoring only for that outer fold", {
  d <- nested_test_data()
  first <- nested_test_run(d)
  fold <- 1L
  idx <- which(first$folds$outer == fold)
  changed <- d
  changed$Y[idx, 1] <- 1 - changed$Y[idx, 1]
  second <- nested_test_run(changed, first$folds)
  expect_identical(first$inner_tuning_results[[fold]], second$inner_tuning_results[[fold]])
  expect_identical(first$models[[fold]], second$models[[fold]])
  expect_identical(first$selected_features_per_outer_fold[[fold]],
    second$selected_features_per_outer_fold[[fold]])
  expect_identical(first$selected_ncomp_per_outer_fold[fold],
    second$selected_ncomp_per_outer_fold[fold])
  a <- subset(first$outer_predictions, outer_fold == fold)
  b <- subset(second$outer_predictions, outer_fold == fold)
  expect_identical(a$Y_pred, b$Y_pred)
  expect_false(isTRUE(all.equal(first$outer_fold_metrics$PRESS[fold],
    second$outer_fold_metrics$PRESS[fold])))
})
test_that("extreme OUTER TEST X cannot affect development or training preprocessing", {
  d <- nested_test_data()
  first <- nested_test_run(d)
  idx <- which(first$folds$outer == 1)
  changed <- d
  changed$X[idx, , ] <- changed$X[idx, , ] * 100 + 1000
  second <- nested_test_run(changed, first$folds)
  expect_identical(first$inner_tuning_results[[1]], second$inner_tuning_results[[1]])
  expect_identical(first$models[[1]], second$models[[1]])
  expect_identical(first$models[[1]]$vip, second$models[[1]]$vip)
  expect_identical(first$preprocessing[[1]], second$preprocessing[[1]])
  a <- subset(first$outer_predictions, outer_fold == 1)
  b <- subset(second$outer_predictions, outer_fold == 1)
  expect_false(isTRUE(all.equal(a$Y_pred, b$Y_pred)))
})
test_that("fold-specific scaling and selected identities match direct training-only fits", {
  d <- nested_test_data()
  r <- nested_test_run(d)
  X <- matrix(d$X, dim(d$X)[1])
  features <- dimnames(d$X)[[2]]
  times <- dimnames(d$X)[[3]]
  for (f in 1:3) {
    tr <- which(r$folds$outer != f)
    model <- r$models[[f]]
    expect_equal(unname(model$imputer$means), unname(colMeans(X[tr, , drop = FALSE])))
    expect_equal(unname(model$fit$preprocessing$x_mean), unname(colMeans(X[tr, model$columns,
      drop = FALSE])))
    expect_equal(unname(model$fit$preprocessing$x_sd), unname(apply(X[tr, model$columns,
      drop = FALSE], 2, sd)))
    expect_false(isTRUE(all.equal(unname(model$fit$preprocessing$x_mean), unname(colMeans(X[,
      model$columns, drop = FALSE])))))
    base <- .tensorpls_fit(X[tr, , drop = FALSE], d$Y[tr, , drop = FALSE], model$config$ncomp)
    sets <- vip_feature_sets(.nested_vips(base, features, times), thr = 65)
    expect_setequal(model$selected_features, sets$sets$VIP2D_intersect_VIP3Dmodel2)
    expect_equal(model$vip, .nested_vips(base, features, times))
  }
})
test_that("complete nested analysis is reproducible and Q2 uses pooled held-out predictions", {
  d <- nested_test_data()
  a <- nested_test_run(d)
  b <- nested_test_run(d)
  expect_identical(a, b)
  expect_true(all(a$outer_fold_metrics$status == "ok"))
  p <- a$outer_predictions
  expect_equal(nrow(p), 36L)
  expect_equal(length(unique(p$subject_id)), 36L)
  press <- sum((p$Y_true - p$Y_pred)^2)
  denom <- sum((p$Y_true - p$training_mean)^2)
  expect_equal(a$Q2Y, 1 - press / denom)
  expect_equal(unname(a$PRESS), press)
  expect_equal(unname(a$Q2_denominator), denom)
  expect_equal(unname(a$classification_metrics["accuracy"]),
    mean(p$true_class == p$predicted_class))
  expect_equal(unname(a$classification_metrics["BER"]),
    1 - unname(a$classification_metrics["balanced_accuracy"]))
  expect_true(all(a$stability$frequency$frequency >= 0 & a$stability$frequency$frequency <=
    1))
  expect_equal(a$stability$jaccard, t(a$stability$jaccard))
  expect_equal(unname(diag(a$stability$jaccard)), rep(1, 3))
})
test_that("fixed and tuned strategies use inner selection only, with shared folds", {
  d <- nested_test_data()
  ids <- dimnames(d$X)[[1]]
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, ids)
  original <- nested_cv_tensorpls(d$X, d$Y, folds = folds, ncomp = 1:2,
    feature_strategy = "Original", input_provenance = "raw", verbose = FALSE)
  tuned <- nested_cv_tensorpls(d$X, d$Y, folds = folds, ncomp = 1:2, feature_strategy = c("VIP2D",
    "VIP3Dmodel2"), vip_thr = c(50, 70), strategy_mode = "tune", input_provenance = "raw",
    verbose = FALSE)
  expect_identical(original$folds, tuned$folds)
  expect_equal(unname(original$number_of_selected_features_per_fold), rep(14L, 3))
  for (f in 1:3) {
    table <- tuned$inner_tuning_results[[f]]$tuning
    chosen <- tuned$models[[f]]$config
    idx <- table$strategy == chosen$strategy & table$thr == chosen$thr & table$ncomp ==
      chosen$ncomp
    expect_equal(table$selection_score[idx], max(table$selection_score, na.rm = TRUE))
    expect_setequal(tuned$inner_tuning_results[[f]]$inner_results$inner_fold, 1:3)
  }
  expect_error(nested_cv_tensorpls(d$X, d$Y, feature_strategy = c("VIP2D", "VIP3Dmodel2")),
    "Fixed mode")
})
test_that("training-only mean imputation handles missing and unusable columns", {
  d <- nested_test_data()
  d$X[1:4, 3, 1] <- NA
  d$X[, 14, 3] <- NA
  r <- nested_test_run(d, imputation = "training_mean")
  expect_true(all(r$outer_fold_metrics$status == "ok"))
  f <- 1
  tr <- which(r$folds$outer != f)
  X <- matrix(d$X, 36)
  expected <- colMeans(X[tr, , drop = FALSE], na.rm = TRUE)
  expected[is.nan(expected)] <- 0
  expect_equal(unname(r$models[[f]]$imputer$means), unname(expected))
  expect_true(tail(r$models[[f]]$imputer$all_missing, 1))
  expect_true(r$convergence[[f]]$imputation$converged)
  expect_equal(r$convergence[[f]]$imputation$iterations, 0L)
  expect_error(nested_test_run(d), "Raw missing X")
  changed <- d
  changed$X[which(r$folds$outer == f), 14, 3] <- 1e9
  r2 <- nested_test_run(changed, r$folds, imputation = "training_mean")
  expect_identical(r$models[[f]], r2$models[[f]])
  expect_equal(subset(r$outer_predictions, outer_fold == f)$Y_pred, subset(r2$outer_predictions,
    outer_fold == f)$Y_pred)
})
test_that("multiclass and repeated dummy responses are explicit and aligned", {
  d <- nested_test_data()
  Y <- matrix(0, 36, 3, dimnames = list(dimnames(d$X)[[1]], c("A",
    "B", "C")))
  Y[cbind(1:36, rep(1:3, 12))] <- 1
  r <- nested_cv_tensorpls(d$X, Y, outer_folds = 3, inner_folds = 3, ncomp = 1:2,
    feature_strategy = "Original", input_provenance = "raw", verbose = FALSE)
  expect_equal(nrow(r$outer_predictions), 108L)
  expect_equal(length(r$Q2Y_per_response), 3L)
  expect_equal(r$Q2Y, mean(r$Q2Y_per_response))
  expect_true(is.na(r$classification_metrics["AUC"]))
  expect_true(all(r$outer_predictions$predicted_class %in% c("A", "B", "C")))
  Y3 <- array(rep(d$Y, 3), c(36, 1, 3), dimnames = list(dimnames(d$X)[[1]],
    "Case", dimnames(d$X)[[3]]))
  repeated <- d
  repeated$Y <- Y3
  expect_equal(nested_test_run(d)$outer_predictions, nested_test_run(repeated)$outer_predictions)
  Y3[1, 1, 2] <- 1 - Y3[1, 1, 2]
  expect_error(nested_cv_tensorpls(d$X, Y3), "time-invariant")
  expect_error(nested_cv_tensorpls(d$X, d$Y[36:1, , drop = FALSE]), "identical order")
})
test_that("invalid and failed candidates are visible, not silently dropped from outer metrics", {
  d <- nested_test_data()
  bad <- nested_cv_tensorpls(d$X, d$Y, outer_folds = 3, inner_folds = 3,
    ncomp = 10, vip_thr = 100, input_provenance = "raw", verbose = FALSE)
  expect_true(all(bad$outer_fold_metrics$status == "failed"))
  expect_true(is.na(bad$Q2Y))
  expect_true(all(is.na(bad$outer_predictions$Y_pred)))
  expect_length(bad$errors, 3L)
  expect_true(all(is.na(bad$stability$frequency$frequency)))
  expect_true(all(is.na(bad$stability$jaccard)))
  expect_true(any(bad$inner_tuning_results[[1]]$inner_results$status == "failed"))
  expect_warning(unknown <- nested_cv_tensorpls(d$X, d$Y, outer_folds = 3, inner_folds = 3,
    ncomp = 1, feature_strategy = "Original", verbose = FALSE), "provenance")
  expect_false(unknown$settings$inductive_conditional_on_raw_input)
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, dimnames(d$X)[[1]])
  names(folds$inner[[1]])[1] <- names(folds$outer)[folds$outer == 1][1]
  expect_error(nested_test_run(d, folds), "OUTER TRAIN")
})
test_that("binary classification metrics and tied AUC match hand calculations", {
  enc <- list(binary = TRUE, levels = c("0", "1"))
  y <- matrix(c(0, 0, 1, 1), ncol = 1)
  pred <- matrix(c(.1, .7, .7, .9), ncol = 1)
  m <- .nested_metrics(y, pred, matrix(.5, 4, 1), as.character(y[,
    1]), enc)
  expect_equal(unname(m$classification), c(.75, .75, 1, .5, .25,
    .875))
  expect_equal(m$Q2Y, 1 - sum((y - pred)^2))
})
test_that("synthetic fixture integrates nested selection with shared folds", {
  d <- readRDS(system.file("extdata", "vip_example.rds", package = "TensorPLS"))
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, dimnames(d$X)[[1]])
  original <- nested_cv_tensorpls(d$X, d$Y, folds = folds, ncomp = 2:3,
    feature_strategy = "Original", input_provenance = "raw", verbose = FALSE)
  vip <- nested_cv_tensorpls(d$X, d$Y, folds = folds, ncomp = 2:3, vip_thr = 50,
    input_provenance = "raw", verbose = FALSE)
  expect_identical(original$folds, vip$folds)
  expect_true(all(original$outer_fold_metrics$status == "ok"))
  expect_true(all(vip$outer_fold_metrics$status == "ok"))
  expect_true(is.finite(original$Q2Y))
  expect_true(is.finite(vip$Q2Y))
  expect_equal(nrow(vip$outer_predictions), 24L)
  expect_true(all(vip$number_of_selected_features_per_fold <= 16L))
  expect_s3_class(vip, "nested_tensorpls")
})
test_that("fit boundary receives only INNER TRAIN or OUTER TRAIN rows and no held-out encoding", {
  d <- nested_test_data()
  folds <- make_tensorpls_folds(d$Y, 3, 3, 123, dimnames(d$X)[[1]])
  original_fit <- .nested_fit_pipeline
  original_develop <- .nested_develop
  calls <- list()
  local_mocked_bindings(.nested_fit_pipeline = function(X, Y, features, times, config,
    imputation) {
    calls[[length(calls) + 1L]] <<- match(X[, 1], d$X[, 1, 1])
    original_fit(X, Y, features, times, config, imputation)
  }, .nested_develop = function(X, Y, classes, encoding, features, times, inner, configurations,
    imputation, selection_metric) {
    expect_named(encoding, c("levels", "binary"))
    expect_equal(length(classes), nrow(X))
    expect_equal(nrow(Y), nrow(X))
    original_develop(X, Y, classes, encoding, features, times, inner, configurations,
      imputation, selection_metric)
  })
  r <- nested_test_run(d, folds)
  expected <- list()
  for (f in 1:3) {
    tr <- which(folds$outer != f)
    for (c in 1:2)
      for (j in 1:3)
      expected[[length(expected) +
      1L]] <- tr[folds$inner[[f]] != j]
    expected[[length(expected) + 1L]] <- tr
  }
  expect_equal(lapply(calls, unname), lapply(expected, unname))
  expect_true(all(r$outer_fold_metrics$status == "ok"))
})
