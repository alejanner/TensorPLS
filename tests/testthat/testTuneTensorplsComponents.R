.tune_fixture <- function() {
  set.seed(72)
  X <- array(rnorm(36 * 14 * 3), c(36, 14, 3), dimnames = list(paste0('S',
    1:36), paste0('G', 1:14), paste0('T', 1:3)))
  Y <- matrix(rep(0:1, 18), ncol = 1, dimnames = list(dimnames(X)[[1]],
    'Case'))
  X[, 2, ] <- X[, 2, ] + 2 * Y[, 1]
  folds <- setNames(rep(rep(1:3, each = 2), 6), dimnames(X)[[1]])
  list(X = X, Y = Y, folds = folds)
}
test_that('ordinary CV reuses pooled development scores and preserves RNG', {
  d <- .tune_fixture()
  set.seed(18)
  rng <- .Random.seed
  z <- tune_tensorpls_components(d$X, d$Y, ncomp = c(3, 1, 2), cv_folds = 3,
    folds = d$folds, verbose = FALSE)
  expect_s3_class(z, 'tensorpls_component_tuning')
  expect_identical(.Random.seed, rng)
  expect_identical(z$tuning$ncomp, 1:3)
  expect_true(all(is.na(z$tuning$thr)))
  expect_equal(nrow(z$fold_results), 9L)
  expect_equal(z$selected_ncomp, z$tuning$ncomp[which.max(z$tuning$selection_score)])
  enc <- .nested_response(d$Y)
  Xm <- matrix(d$X, 36)
  colnames(Xm) <- as.vector(outer(dimnames(d$X)[[2]], dimnames(d$X)[[3]],
    paste, sep = '_'))
  direct <- .nested_develop(Xm, enc$Y, enc$class, enc[c('levels', 'binary')], dimnames(d$X)[[2]],
    dimnames(d$X)[[3]], d$folds, z$tuning[, c('strategy', 'thr', 'ncomp')], 'none',
    'Q2Y')
  expect_identical(z$tuning, direct$tuning)
  expect_identical(z$fold_results, direct$inner_results)
  expect_equal(z$selected_ncomp, direct$model$config$ncomp)
  a <- tune_tensorpls_components(d$X, d$Y, ncomp = 1:2, cv_folds = 3,
    verbose = FALSE)
  b <- tune_tensorpls_components(d$X, d$Y, ncomp = 1:2, cv_folds = 3,
    verbose = FALSE)
  expect_identical(a, b)
  expect_identical(.Random.seed, rng)
})
test_that('ties and failed development keep shared-engine behavior', {
  d <- .tune_fixture()
  local_mocked_bindings(.nested_metrics = function(...)
    list(Q2Y = .25, classification = c(BER = .5)),
    .package = 'TensorPLS')
  z <- tune_tensorpls_components(d$X, d$Y, ncomp = c(3, 2, 1), cv_folds = 3,
    folds = d$folds, verbose = FALSE)
  expect_identical(z$selected_ncomp, 1L)
  local_mocked_bindings(.nested_fit_pipeline = function(...)
    stop('fixture failure'), .package = 'TensorPLS')
  z <- tune_tensorpls_components(d$X, d$Y, ncomp = 1, cv_folds = 3, folds = d$folds,
    verbose = FALSE)
  expect_true(is.na(z$selected_ncomp))
  expect_null(z$refit)
  expect_match(z$diagnostics$error, 'No candidate')
  expect_true(all(!is.finite(z$tuning$selection_score)))
})
test_that('selection and mean imputation occur within training splits before final refit', {
  d <- .tune_fixture()
  d$X[1, 3, 2] <- NA_real_
  seen <- selections <- list()
  current <- NULL
  original_imputer <- .nested_imputer
  original_selection <- vip_feature_sets
  local_mocked_bindings(.nested_imputer = function(X, method) {
    current <<- match(X[, 1], d$X[, 1, 1])
    seen[[length(seen) + 1]] <<- current
    original_imputer(X, method)
  }, vip_feature_sets = function(...) {
    selections[[length(selections) + 1]] <<- current
    original_selection(...)
  }, .package = 'TensorPLS')
  z <- tune_tensorpls_components(d$X, d$Y, ncomp = 1:2,
    feature_strategy = 'VIP2D_intersect_VIP3Dmodel2',
    vip_thr = 50, cv_folds = 3, folds = d$folds, imputation = 'training_mean',
    verbose = FALSE)
  expect_null(z$diagnostics$error)
  expected <- c(rep(lapply(1:3, function(f)
    unname(which(d$folds != f))),
    2), list(1:36))
  expect_identical(seen, expected)
  expect_identical(selections, expected)
  expect_true(all(z$tuning$strategy == 'VIP2D_intersect_VIP3Dmodel2'))
  expect_true(all(z$tuning$thr == 50))
  expect_equal(unname(z$refit$imputer$means), colMeans(matrix(d$X, 36), na.rm = TRUE))
})
test_that('Tucker uses the existing train/apply route and development cache', {
  d <- .tune_fixture()
  trained <- applied <- list()
  local_mocked_bindings(fit_tucker_imputer = function(X_train, ...) {
    trained[[length(trained) + 1]] <<- dimnames(X_train)[[1]]
    structure(list(X_train_imputed = X_train, training_ids = dimnames(X_train)[[1]],
      convergence = list(converged = TRUE)), class = 'tucker_imputer')
  }, predict_tucker_imputer = function(object, X, ...) {
    expect_length(intersect(object$training_ids, dimnames(X)[[1]]), 0)
    applied[[length(applied) + 1]] <<- dimnames(X)[[1]]
    X
  }, .package = 'TensorPLS')
  z <- tune_tensorpls_components(d$X, d$Y, ncomp = 1:2, cv_folds = 3,
    folds = d$folds, imputation = 'tucker', tucker_control = list(ranks = c(2L,
      2L, 2L)), verbose = FALSE)
  ids <- dimnames(d$X)[[1]]
  expect_identical(trained, c(lapply(1:3, function(f)
    ids[d$folds != f]),
    list(ids)))
  expect_identical(applied, rep(lapply(1:3, function(f)
    ids[d$folds == f]),
    2))
  expect_null(z$diagnostics$error)
  expect_identical(z$settings$tucker_control$seed, 123L)
})
test_that('invalid inputs and misaligned folds fail before development', {
  d <- .tune_fixture()
  call <- function(...)
    tune_tensorpls_components(d$X, d$Y, verbose = FALSE, ...)
  for (h in list(numeric(), 0, 11, 1.5, NA_real_, Inf, '2', 1i))
    expect_error(call(ncomp = h),
    'ncomp')
  expect_error(call(feature_strategy = c('Original', 'VIP2D')), 'one supported')
  expect_error(call(vip_thr = c(95, 99)), 'one percentile')
  expect_error(call(cv_folds = 1), 'cv_folds')
  for (f in list(unname(d$folds), rev(d$folds), d$folds + 1, replace(d$folds,
    1, NA), d$folds[-1]))
      expect_error(call(ncomp = 1, cv_folds = 3,
    folds = f), 'folds')
  bad <- setNames(rep(1:2, 18), names(d$folds))
  expect_error(call(ncomp = 1, cv_folds = 2, folds = bad), 'all classes')
  y <- d$Y
  rownames(y) <- rev(rownames(y))
  expect_error(tune_tensorpls_components(d$X, y), 'identical order')
  x <- d$X
  x[1, 2, 1] <- NA
  expect_error(tune_tensorpls_components(x, d$Y), 'Missing X')
  expect_warning(call(ncomp = 1, cv_folds = 3, input_provenance = 'preimputed'),
    'not raw')
})
