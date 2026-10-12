.public_pls_fixture <- function() {
  set.seed(812)
  X <- array(rnorm(12 * 4 * 5), c(12, 4, 5), dimnames = list(paste0('s',
    1:12), paste0('gene_', 1:4), paste0('time', 1:5)))
  Y <- matrix(rep(c(0, 1), 6), ncol = 1, dimnames = list(dimnames(X)[[1]],
    'Case'))
  X[, 1, ] <- X[, 1, ] + Y[, 1]
  list(X = X, Y = Y)
}
test_that('the public API calls the native core once with feature-fast unfolding', {
  d <- .public_pls_fixture()
  native <- .tensorpls_fit
  calls <- 0L
  captured <- NULL
  forbidden <- function(...)
    stop('unexpected extra pipeline')
  local_mocked_bindings(.tensorpls_fit = function(X, Y, ncomp, ...) {
    calls <<- calls + 1L
    captured <<- list(X = X, Y = Y, ncomp = ncomp)
    native(X, Y, ncomp, ...)
  }, nested_cv_tensorpls = forbidden, fit_tucker_imputer = forbidden, vip_feature_sets = forbidden,
    .nested_vips = forbidden, .package = 'TensorPLS')
  local_mocked_bindings(plsda = forbidden, pls = forbidden, .package = 'mixOmics')
  z <- fit_tensorpls_pls(d$X, d$Y, 3)
  expect_identical(calls, 1L)
  expect_identical(captured$Y, d$Y)
  expect_identical(captured$ncomp, 3L)
  for (ti in 1:5)
    for (f in 1:4)
    expect_identical(unname(captured$X[,
    f + 4 * (ti - 1)]), unname(d$X[, f, ti]))
  expect_identical(colnames(captured$X), as.vector(outer(dimnames(d$X)[[2]], dimnames(d$X)[[3]],
    paste, sep = '_')))
  expect_identical(z$x_scores, z$core_fit$x.scores)
  expect_null(z$core_fit$Q2)
  expect_null(z$core_fit$Q2cum)
  expect_false(any(c('VIP2D', 'VIP3Dmodel1', 'VIP3Dmodel2', 'accuracy', 'AUC') %in% names(z)))
})
test_that('the public object preserves the complete direct native fit and labels', {
  d <- .public_pls_fixture()
  Xm <- matrix(d$X, 12)
  colnames(Xm) <- as.vector(outer(dimnames(d$X)[[2]], dimnames(d$X)[[3]],
    paste, sep = '_'))
  direct <- .tensorpls_fit(Xm, d$Y, 3L)
  z <- fit_tensorpls_pls(d$X, d$Y, 3)
  expect_s3_class(z, 'tensorpls_model')
  expect_identical(z$core_fit, direct)
  expect_identical(z$x_scores, direct$x.scores)
  expect_identical(z$x_weights, direct$raw.wgs)
  expect_identical(z$x_loadings, direct$x.loads)
  expect_identical(z$coefficients, direct$reg.coefs)
  expect_identical(z$fitted_y, direct$y.pred)
  expect_identical(z$preprocessing, direct$preprocessing)
  expect_identical(z$convergence, direct$convergence)
  expect_identical(z$subject_ids, dimnames(d$X)[[1]])
  expect_identical(z$feature_names, dimnames(d$X)[[2]])
  expect_identical(z$time_names, dimnames(d$X)[[3]])
  expect_identical(z$ncomp, 3L)
  expect_equal(z$fitted_R2Y, 1 - sum((d$Y - direct$y.pred)^2) / sum((d$Y -
    mean(d$Y))^2))
  expect_identical(names(z$fitted_R2Y_per_response), 'Case')
  expect_identical(z$response_encoding, list(levels = c('0', '1'), binary = TRUE,
    response_names = 'Case'))
  expect_named(z, c('ncomp', 'feature_names', 'time_names', 'subject_ids',
    'training_class', 'x_scores', 'x_weights', 'x_loadings',
    'coefficients', 'fitted_y', 'fitted_R2Y', 'fitted_R2Y_per_response',
      'response_encoding', 'preprocessing', 'convergence', 'core_fit'))
})
test_that('classification responses follow existing TensorPLS encoding', {
  d <- .public_pls_fixture()
  ref <- fit_tensorpls_pls(d$X, d$Y, 2)
  f <- factor(ifelse(d$Y[, 1] == 1, 'Case', 'Control'), levels = c('Control',
    'Case'))
  names(f) <- rownames(d$Y)
  z <- fit_tensorpls_pls(d$X, f, 2)
  expect_identical(z$x_scores, ref$x_scores)
  expect_identical(z$fitted_y, ref$fitted_y)
  expect_identical(z$response_encoding$levels, c('Control', 'Case'))
  Y3 <- array(rep(d$Y, 5), c(12, 1, 5), dimnames = list(rownames(d$Y),
    'Case', dimnames(d$X)[[3]]))
  expect_identical(fit_tensorpls_pls(d$X, Y3, 2)$core_fit, ref$core_fit)
  Xm <- matrix(d$X, 12)
  colnames(Xm) <- as.vector(outer(dimnames(d$X)[[2]], dimnames(d$X)[[3]],
    paste, sep = '_'))
  group <- factor(rep(c('C', 'A', 'B'), 4), levels = c('C', 'A',
    'B'))
  onehot <- cbind(C = as.numeric(group == 'C'), A = as.numeric(group == 'A'),
    B = as.numeric(group == 'B'))
  zz <- fit_tensorpls_pls(d$X, group, 2)
  direct <- .tensorpls_fit(Xm, onehot, 2L)
  expect_identical(zz$core_fit, direct)
  expect_identical(fit_tensorpls_pls(d$X, onehot, 2)$core_fit, direct)
  expected <- 1 - colSums((onehot - direct$y.pred)^2) / colSums(sweep(onehot,
    2, colMeans(onehot), '-')^2)
  expect_equal(zz$fitted_R2Y_per_response, expected)
  expect_equal(zz$fitted_R2Y, mean(expected))
  expect_false(zz$response_encoding$binary)
  expect_identical(zz$response_encoding$response_names, c('C', 'A', 'B'))
  bad <- Y3
  bad[1, 1, 2] <- 1 - bad[1, 1, 2]
  expect_error(fit_tensorpls_pls(d$X, bad, 2), 'time-invariant')
  expect_error(fit_tensorpls_pls(d$X, matrix(0.2, 12, 1), 2), 'binary or one-hot')
  expect_error(fit_tensorpls_pls(d$X, rep('one class', 12), 2), 'At least two classes')
})
test_that('missingness, invalid dimensions and ncomp fail before any fit', {
  d <- .public_pls_fixture()
  local_mocked_bindings(.tensorpls_fit = function(...)
    stop('core should not run'), .package = 'TensorPLS')
  for (missing in c(NA_real_, NaN)) {
    X <- d$X
    X[1, 1, 1] <- missing
    expect_error(fit_tensorpls_pls(X, d$Y, 2), 'impute X first')
  }
  X <- d$X
  X[1, 1, 1] <- Inf
  expect_error(fit_tensorpls_pls(X, d$Y, 2), 'finite')
  expect_error(fit_tensorpls_pls(matrix(d$X, 12), d$Y, 2), '3D tensor')
  expect_error(fit_tensorpls_pls(array('a', c(12, 4, 5)), d$Y,
    2), 'real numeric')
  expect_error(fit_tensorpls_pls(array(1, c(1, 4, 5)), d$Y,
    1), 'at least two subjects')
  expect_error(fit_tensorpls_pls(array(numeric(), c(12, 0, 5)), d$Y,
    1), 'nonempty')
  for (nc in list(0, -1, 1.5, NA_real_, NaN, Inf, 1:2,
    '2', TRUE, NULL, 2 + 1i))
      expect_error(fit_tensorpls_pls(d$X,
    d$Y, nc), 'ncomp')
  expect_error(fit_tensorpls_pls(d$X, d$Y, 12), 'ncomp must not exceed')
  expect_error(fit_tensorpls_pls(d$X[, 1, 1:2, drop = FALSE], d$Y,
    3), 'ncomp must not exceed')
  expect_error(fit_tensorpls_pls(d$X, d$Y), 'ncomp.*missing')
})
test_that('labels align without reordering, with documented unnamed defaults', {
  d <- .public_pls_fixture()
  bad <- d$Y
  rownames(bad) <- rev(rownames(bad))
  expect_error(fit_tensorpls_pls(d$X, bad, 2), 'subject labels.*order')
  expect_error(fit_tensorpls_pls(d$X, d$Y[-1, , drop = FALSE], 2), 'subject')
  badx <- d$X
  dimnames(badx)[[2]][2] <- dimnames(badx)[[2]][1]
  expect_error(fit_tensorpls_pls(badx, d$Y, 2), 'unique, nonempty')
  y3 <- array(rep(d$Y, 5), c(12, 1, 5), dimnames = list(rownames(d$Y),
    'Case', rev(dimnames(d$X)[[3]])))
  expect_error(fit_tensorpls_pls(d$X, y3, 2), 'time axis')
  expect_error(fit_tensorpls_pls(d$X, y3[, , 1:4, drop = FALSE], 2),
    'time axis')
  x <- d$X
  dimnames(x) <- NULL
  z <- fit_tensorpls_pls(x, d$Y, 1)
  expect_identical(z$subject_ids, rownames(d$Y))
  expect_identical(z$feature_names, paste0('Var', 1:4))
  expect_identical(z$time_names, paste0('Time', 1:5))
  expect_identical(z$ncomp, 1L)
  yy <- d$Y
  rownames(yy) <- NULL
  expect_identical(fit_tensorpls_pls(x, yy, 1)$subject_ids, as.character(1:12))
  y3 <- array(rep(d$Y, 5), c(12, 1, 5), dimnames = list(rownames(d$Y),
    'Case', paste0('visit', 1:5)))
  expect_identical(fit_tensorpls_pls(x, y3, 1)$time_names, paste0('visit', 1:5))
})
test_that('native non-convergence is surfaced without a second fit', {
  d <- .public_pls_fixture()
  native <- .tensorpls_fit
  calls <- 0L
  local_mocked_bindings(.tensorpls_fit = function(...) {
    calls <<- calls + 1L
    z <- native(...)
    z$convergence$converged[1] <- FALSE
    z
  }, .package = 'TensorPLS')
  expect_warning(z <- fit_tensorpls_pls(d$X, d$Y, 2), 'inspect convergence')
  expect_identical(calls, 1L)
  expect_false(z$convergence$converged[1])
  expect_identical(z$x_scores, z$core_fit$x.scores)
})
test_that('training classes preserve original coding, levels and subject names', {
  d <- .public_pls_fixture()
  ids <- dimnames(d$X)[[1]]
  binary <- d$Y
  binary_factor <- factor(ifelse(binary[, 1] == 1, 'Case', 'Control'), levels = c('Case',
    'Control'))
  multi <- factor(rep(c('C', 'A', 'B'), 4), levels = c('B', 'C',
    'A'))
  onehot <- cbind(B = as.numeric(multi == 'B'), C = as.numeric(multi == 'C'),
    A = as.numeric(multi == 'A'))
  repeated <- array(rep(binary, 5), c(12, 1, 5))
  cases <- list(binary = binary, vector = as.vector(binary), factor = binary_factor,
    multiclass = onehot, multifactor = multi, repeated = repeated)
  expected <- list(binary = factor(binary[, 1], levels = c('0', '1')),
    vector = factor(binary[, 1], levels = c('0', '1')), factor = binary_factor,
    multiclass = multi, multifactor = multi, repeated = factor(binary[, 1],
      levels = c('0', '1')))
  for (nm in names(cases)) {
    if (nm %in% c('multiclass', 'multifactor')) {
      expect_warning(model <- fit_tensorpls_pls(d$X, cases[[nm]], 2),
        '^TensorPLS fit did not converge for every component; inspect convergence\\.$')
    } else {
      model <- fit_tensorpls_pls(d$X, cases[[nm]], 2)
    }
    expect_identical(model$training_class, stats::setNames(expected[[nm]], ids))
    expect_identical(names(model$training_class), model$subject_ids)
    expect_identical(levels(model$training_class), model$response_encoding$levels)
  }
})
