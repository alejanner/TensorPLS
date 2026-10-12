.public_vips_fixture <- function(p = 4L, k = 5L, H = 3L) {
  set.seed(471)
  X <- array(rnorm(20 * p * k), c(20, p, k), dimnames = list(paste0('s',
    1:20), paste0('gene', seq_len(p)), paste0('time', seq_len(k))))
  fit_tensorpls_pls(X, factor(rep(c('Control', 'Case'), 10)), H)
}
test_that('VIP representations never invoke fitting or change the supplied model', {
  model <- .public_vips_fixture()
  before <- serialize(model, NULL)
  forbidden <- function(...)
    stop('unexpected fitting call')
  local_mocked_bindings(.tensorpls_fit = forbidden, fit_tensorpls_pls = forbidden,
    nested_cv_tensorpls = forbidden,
    .package = 'TensorPLS')
  if (requireNamespace('mixOmics', quietly = TRUE))
    local_mocked_bindings(pls = forbidden, plsda = forbidden, .package = 'mixOmics')
  z <- tensorpls_vips(model)
  expect_s3_class(z, 'tensorpls_vips')
  expect_identical(serialize(model, NULL), before)
  expect_identical(unname(z$VIP2D), unname(model$core_fit$VIP))
  expect_named(z, c('V', 'VIP2D', 'VIP3Dmodel1', 'VIP3Dmodel2', 'M2_mean',
    'feature_names', 'time_names', 'component_names',
    'ncomp'))
})
test_that('coordinates, historical sums, all-time means and singleton axes are preserved', {
  for (dims in list(c(4L, 5L, 3L), c(4L, 3L, 3L), c(3L,
    2L, 4L), c(1L, 5L, 2L), c(4L, 1L, 2L),
    c(4L, 5L, 1L), c(1L, 1L, 1L))) {
    p <- dims[1]
    k <- dims[2]
    H <- dims[3]
    model <- .public_vips_fixture(p, k, H)
    z <- tensorpls_vips(model)
    native <- .nested_vips(model$core_fit, model$feature_names, model$time_names)
    expect_identical(z$VIP2D, native$VIP2D)
    expect_identical(z$VIP3Dmodel1, native$VIP3Dmodel1)
    expect_identical(z$M2_mean, native$M2_mean)
    expect_identical(z$VIP3Dmodel2, native$VIP3Dmodel2)
    expect_identical(dim(z$VIP3Dmodel2), c(p, k))
    expect_identical(dimnames(z$VIP3Dmodel2), list(model$feature_names, model$time_names))
    expect_equal(max(abs(z$VIP3Dmodel2 - H * z$M2_mean)), 0, tolerance = 1e-14)
    expect_equal(z$M2_mean, z$VIP3Dmodel2 / H, tolerance = 1e-14)
    expect_identical(max(abs(z$VIP2D - model$core_fit$VIP)), 0)
    expect_identical(dim(z$V), c(p, k, H))
    expect_identical(dim(z$VIP3Dmodel1), c(p, H))
    expect_identical(dim(z$M2_mean), c(p, k))
    expect_identical(z$feature_names, model$feature_names)
    expect_identical(z$time_names, model$time_names)
    expect_identical(z$component_names, paste0('Comp', seq_len(H)))
    expect_identical(z$ncomp, H)
    expect_identical(dimnames(z$V), list(model$feature_names, model$time_names, z$component_names))
    expect_identical(dimnames(z$VIP3Dmodel1), list(model$feature_names, paste0(z$component_names,
      '_meanTime')))
    expect_identical(dimnames(z$M2_mean), list(model$feature_names, model$time_names))
    for (f in seq_len(p))
      for (t in seq_len(k))
      for (h in seq_len(H))
      expect_identical(unname(z$V[f,
      t, h]), unname(model$core_fit$VIP[f + p * (t - 1), h]))
    for (f in seq_len(p))
      for (h in seq_len(H))
      expect_identical(unname(z$VIP3Dmodel1[f,
      h]), mean(z$V[f, seq_len(k), h]))
    for (f in seq_len(p))
      for (t in seq_len(k)) {
      expect_identical(unname(z$M2_mean[f, t]), mean(z$V[f, t, seq_len(H)]))
      expect_identical(unname(z$VIP3Dmodel2[f, t]), sum(z$V[f, t, seq_len(H)]))
    }
  }
})
test_that('late times and every cumulative path column contribute to the means', {
  model <- .public_vips_fixture()
  z <- tensorpls_vips(model)
  # Perturb a stored value solely to isolate the aggregation contract.
  for (t in 4:5) {
    changed <- model
    changed$core_fit$VIP[1 + 4 * (t - 1), 2] <- changed$core_fit$VIP[1 +
      4 * (t - 1), 2] + 5
    zz <- tensorpls_vips(changed)
    expect_equal(unname(zz$VIP3Dmodel1[1, 2] - z$VIP3Dmodel1[1,
      2]), 1)
  }
  for (h in 1:3) {
    changed <- model
    changed$core_fit$VIP[1, h] <- changed$core_fit$VIP[1, h] + 3
    zz <- tensorpls_vips(changed)
    expect_equal(unname(zz$M2_mean[1, 1] - z$M2_mean[1, 1]),
      1)
    expect_equal(unname(zz$VIP3Dmodel2[1, 1] - z$VIP3Dmodel2[1,
      1]), 3)
  }
})
test_that('invalid or inconsistent model objects fail with informative errors', {
  model <- .public_vips_fixture()
  for (x in list(NULL, list(), model$core_fit, unclass(model), 1))
    expect_error(tensorpls_vips(x),
    'tensorpls_model')
  for (nm in c('feature_names', 'time_names'))
    for (x in list(NULL, character(),
    NA_character_, c('a',
    'a'), '', 1:4)) {
    bad <- model
    bad[nm] <- list(x)
    expect_error(tensorpls_vips(bad), 'feature_names and time_names')
  }
  for (x in list(NULL, 0, -1, 1.5, NA_real_, Inf, c(1, 2),
    '3', TRUE, 1 + 1i)) {
    bad <- model
    bad['ncomp'] <- list(x)
    expect_error(tensorpls_vips(bad), 'ncomp')
  }
  bad <- model
  bad$core_fit <- NULL
  expect_error(tensorpls_vips(bad), 'core_fit')
  bad <- model
  bad$core_fit$ncomp <- 2L
  expect_error(tensorpls_vips(bad), 'must agree')
  for (x in list(NULL, 1:60, matrix(1, 19, 3), matrix(1,
    20, 2), matrix('a', 20, 3), matrix(NA_real_, 20,
    3), matrix(Inf, 20, 3), matrix(1 + 1i, 20,
    3))) {
    bad <- model
    bad$core_fit['VIP'] <- list(x)
    expect_error(tensorpls_vips(bad), 'finite numeric')
  }
})
