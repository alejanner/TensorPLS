.impute_reference <- function(name) {
  readRDS(test_path("fixtures", "imputation-api", paste0(name, ".rds")))
}
test_that("mean completion is identical to the sealed pre-change outputs", {
  X <- .impute_reference("X_original")
  completed <- ImputeTensor(X, method = "mean")
  expect_identical(as.vector(completed[is.na(X)]), .impute_reference("mean_imputed_values_before"))
  expect_identical(as.vector(completed), as.vector(.impute_reference("X_mean_before")))
  expect_identical(completed, .impute_reference("X_mean_before"))
  expect_identical(dim(completed), dim(X))
  expect_identical(dimnames(completed), dimnames(X))
  expect_identical(completed[!is.na(X)], X[!is.na(X)])
  expect_false(anyNA(completed))
  predictors <- matrix(X, dim(X)[1])
  old <- .nested_impute(.nested_imputer(predictors, "training_mean"), predictors)
  expect_identical(as.vector(completed), as.vector(old))
})
test_that("Tucker completion is identical to the sealed pre-change outputs", {
  X <- .impute_reference("X_original")
  configuration <- .impute_reference("controls")
  completed <- ImputeTensor(X, "tucker", ranks = configuration$ranks,
    tucker_control = configuration$controls)
  expect_identical(as.vector(completed[is.na(X)]),
    .impute_reference("tucker_imputed_values_before"))
  expect_identical(as.vector(completed), as.vector(.impute_reference("X_tucker_before")))
  expect_identical(completed, .impute_reference("X_tucker_before"))
  expect_identical(dim(completed), dim(X))
  expect_identical(dimnames(completed), dimnames(X))
  expect_identical(completed[!is.na(X)], X[!is.na(X)])
  expect_false(anyNA(completed))
  expect_identical(ImputeTensor(X, "tucker", configuration$ranks), completed)
  old <- do.call(fit_tucker_imputer, c(list(X_train = X, ranks = configuration$ranks),
    configuration$controls))
  expect_identical(old, .impute_reference("tucker_fit_before"))
  expect_identical(completed, old$X_train_imputed)
})
test_that("protected scientific formals and executable bodies are unchanged", {
  reference <- .impute_reference("functions")
  protected <- setdiff(names(reference), c("tune_tucker_ranks", "nested_cv_tensorpls"))
  for (name in protected) {
    current <- get(name, asNamespace("TensorPLS"))
    expect_identical(formals(current), reference[[name]]$formals, info = name)
    expect_identical(deparse(body(current), width.cutoff = 500L), reference[[name]]$body,
      info = name)
  }
})
test_that("mean completion preserves all-missing semantics and array axes", {
  X <- .impute_reference("X_original")
  X[, 2, 3] <- NA_real_
  X[1, 3, 2] <- NaN
  completed <- ImputeTensor(X, "mean", ranks = "irrelevant")
  old_matrix <- matrix(X, dim(X)[1])
  old <- .nested_impute(.nested_imputer(old_matrix, "training_mean"), old_matrix)
  expect_identical(as.vector(completed), as.vector(old))
  expect_identical(completed[, 2, 3], setNames(rep(0, dim(X)[1]),
    dimnames(X)[[1]]))
  expect_identical(completed[!is.na(X)], X[!is.na(X)])
  expect_identical(dimnames(completed), dimnames(X))
  expect_identical(names(dimnames(completed)), names(dimnames(X)))
  expect_false(anyNA(completed))
  expect_identical(ImputeTensor(completed, "mean"), completed)
  all_missing <- array(NA_real_, c(2, 3, 4))
  expect_identical(ImputeTensor(all_missing, "mean"), array(0, c(2, 3,
    4)))
  integers <- array(1:24, c(2, 3, 4))
  predictors <- matrix(integers, 2)
  old_integer_result <- .nested_impute(.nested_imputer(predictors, "training_mean"), predictors)
  expect_identical(as.vector(ImputeTensor(integers, "mean")), as.vector(old_integer_result))
})
test_that("Tucker defaults and complete inputs use the existing fitter", {
  X <- array(2, c(4, 4, 3), dimnames = list(paste0("s",
    1:4), paste0("g", 1:4), 1:3))
  expect_identical(ImputeTensor(X, "tucker"), fit_tucker_imputer(X)$X_train_imputed)
  expect_identical(ImputeTensor(X, "tucker"), X)
  X[3] <- NA_real_
  expect_identical(ImputeTensor(X, "tucker"), fit_tucker_imputer(X)$X_train_imputed)
  expect_identical(ImputeTensor(X, "tucker", c(1L, 1L, 1L)), fit_tucker_imputer(X,
    c(1L, 1L, 1L))$X_train_imputed)
})
test_that("wrapper validation does not silently change ranks or controls", {
  X <- .impute_reference("X_original")
  expect_error(ImputeTensor(X, "other"), "arg")
  expect_error(ImputeTensor(matrix(1, 3, 4)), "tensor")
  invalid <- X
  invalid[2] <- Inf
  expect_error(ImputeTensor(invalid, "mean"), "infinite")
  expect_error(ImputeTensor(invalid, "tucker"), "infinite")
  for (ranks in list(c(1, 1), c(1.5, 2, 2), c(13, 2,
    2), c(3, 1, 1), c(NA, 2, 2)))
      expect_error(ImputeTensor(X,
    "tucker", ranks), "rank|Rank")
  expect_error(ImputeTensor(array(1, c(2, 2, 2)), "tucker"),
    "supply valid ranks or run tune_tucker_ranks")
  expect_error(ImputeTensor(X, "tucker", c(2, 2, 2), list(ranks = c(1,
    1, 1))), "tucker_control")
  expect_error(ImputeTensor(X, "tucker", c(2, 2, 2), list(unknown = 1)),
    "tucker_control")
  expect_error(ImputeTensor(X, "mean", tucker_control = list(ridge = 1)), "only")
})
test_that("the wrapper delegates without fitting PLS or tuning ranks", {
  X <- .impute_reference("X_original")
  blocked <- function(...)
    stop("Unexpected scientific computation")
  local_mocked_bindings(.tensorpls_fit = blocked, tune_tucker_ranks = blocked,
    nested_cv_tensorpls = blocked,
    .package = "TensorPLS")
  seen <- NULL
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ranks) {
    seen <<- list(X = X_train, ranks = ranks)
    X_train[is.na(X_train)] <- 10
    list(X_train_imputed = X_train)
  }, .package = "TensorPLS")
  completed <- ImputeTensor(X, "tucker")
  expect_identical(seen, list(X = X, ranks = c(4L, 4L, 3L)))
  expect_true(all(completed[is.na(X)] == 10))
  local_mocked_bindings(fit_tucker_imputer = blocked, .package = "TensorPLS")
  expect_identical(ImputeTensor(X, "mean"), .impute_reference("X_mean_before"))
})
