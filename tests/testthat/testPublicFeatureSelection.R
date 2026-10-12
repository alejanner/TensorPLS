.selection_vips <- function(p = 4L, k = 3L, H = 2L) {
  set.seed(701)
  X <- array(rnorm(20 * p * k), c(20, p, k), dimnames = list(NULL, paste0('gene_',
    seq_len(p)), paste0('t', seq_len(k))))
  tensorpls_vips(fit_tensorpls_pls(X, factor(rep(c('a', 'b'), 10)), H))
}
.selection_strategies <- c('Original', 'VIP2D', 'VIP3Dmodel1', 'VIP3Dmodel2',
  'VIP2D_intersect_VIP3Dmodel1', 'VIP2D_intersect_VIP3Dmodel2', 'VIP3Dmodel1_intersect_VIP3Dmodel2',
  'VIP2D_intersect_VIP3Dmodel1_intersect_VIP3Dmodel2')
test_that('single-strategy mode exactly preserves existing sets and tensor ordering', {
  for (d in list(c(4L, 5L, 3L), c(4L, 3L, 3L), c(3L,
    2L, 4L), c(1L, 3L, 2L), c(4L, 1L, 2L),
    c(4L, 3L, 1L), c(1L, 1L, 1L))) {
    v <- .selection_vips(d[1], d[2], d[3])
    scaled <- v
    scaled$VIP3Dmodel2 <- v$M2_mean
    for (thr in c(0, 50, 95, 99, 100)) {
      ref <- vip_feature_sets(v, thr = thr)
      for (s in .selection_strategies) {
        z <- vip_feature_sets(v, thr = thr, strategy = s)
        expected <- if (s == 'Original')
          v$feature_names else
            v$feature_names[v$feature_names %in% ref$sets[[s]]]
        expect_identical(z$selected_features, expected)
        expect_identical(z$selected_indices, match(expected, v$feature_names))
        expect_identical(z$percentile, thr / 100)
        zz <- vip_feature_sets(scaled, thr = thr, strategy = s)
        expect_identical(z$selected_features, zz$selected_features)
        expect_identical(z$selected_indices, zz$selected_indices)
        for (nm in names(z$thresholds))
          expect_identical(z$thresholds[[nm]],
          vapply(as.data.frame(v[[nm]]),
          function(x)
            quantile(x, thr / 100, names = FALSE), numeric(1)))
      }
    }
  }
})
test_that('selection uses all columns, inclusive boundaries and exact feature mapping', {
  v <- .selection_vips(4L, 2L, 2L)
  # Controlled score matrices isolate the selection rule from fitting.
  v$VIP2D[, ] <- 0
  v$VIP2D[1, 1] <- 10
  v$VIP2D[2, 2] <- 10
  v$VIP2D[5, 2] <- 10
  v$VIP3Dmodel2[, ] <- 0
  v$VIP3Dmodel2[1, 1] <- 10
  v$VIP3Dmodel2[3, 2] <- 10
  a <- vip_feature_sets(v, thr = 100, strategy = 'VIP2D')
  b <- vip_feature_sets(v, thr = 100, strategy = 'VIP3Dmodel2')
  both <- vip_feature_sets(v, thr = 100, strategy = 'VIP2D_intersect_VIP3Dmodel2')
  expect_identical(a$selected_features, c('gene_1', 'gene_2'))
  expect_identical(b$selected_features, c('gene_1', 'gene_3'))
  expect_identical(both$selected_features, 'gene_1')
  expect_identical(both$selected_indices, 1L)
  expect_identical(vip_feature_sets(v, thr = 0, strategy = 'VIP2D')$selected_features,
    v$feature_names)
  # Changing only an earlier cumulative column can change selection.
  changed <- v
  changed$VIP2D[4, 1] <- 20
  expect_identical(changed$VIP2D[, 2], v$VIP2D[, 2])
  expect_identical(vip_feature_sets(changed, thr = 100, strategy = 'VIP2D')$selected_features,
    c('gene_1', 'gene_2', 'gene_4'))
  v$VIP3Dmodel1[, ] <- c(0, 1, 1, 2)
  z <- vip_feature_sets(v, thr = 50, strategy = 'VIP3Dmodel1')
  expect_identical(z$selected_features, c('gene_2', 'gene_3', 'gene_4'))
  expect_identical(unname(z$thresholds$VIP3Dmodel1), c(1, 1))
})
test_that('selection does not fit, recompute VIPs, or mutate its input', {
  v <- .selection_vips()
  before <- serialize(v, NULL)
  forbidden <- function(...)
    stop('Unexpected computation')
  local_mocked_bindings(fit_tensorpls_pls = forbidden, .tensorpls_fit = forbidden,
    nested_cv_tensorpls = forbidden,
    tensorpls_vips = forbidden, fit_tucker_imputer = forbidden, .package = 'TensorPLS')
  local_mocked_bindings(pls = forbidden, plsda = forbidden, .package = 'mixOmics')
  for (s in .selection_strategies)
    expect_type(vip_feature_sets(v, thr = 99, strategy = s),
    'list')
  expect_identical(serialize(v, NULL), before)
})
test_that('malformed public objects and incompatible settings fail clearly', {
  v <- .selection_vips()
  select <- function(x, ...)
    vip_feature_sets(x, strategy = 'VIP2D', ...)
  expect_error(select(unclass(v)), 'tensorpls_vips')
  for (nm in c('V', 'VIP2D', 'VIP3Dmodel1', 'VIP3Dmodel2', 'M2_mean')) {
    bad <- v
    bad[[nm]] <- NULL
    expect_error(select(bad), 'missing/invalid')
    bad <- v
    bad[[nm]][1] <- Inf
    expect_error(select(bad), 'missing/invalid')
    bad <- v
    bad[[nm]][1] <- NA_real_
    expect_error(select(bad), 'missing/invalid')
    bad <- v
    dimnames(bad[[nm]]) <- NULL
    expect_error(select(bad), 'labels')
    bad <- v
    bad[[nm]] <- matrix(1, 2, 2)
    expect_error(select(bad), 'dimensions')
  }
  for (nm in c('feature_names', 'time_names')) {
    bad <- v
    bad[[nm]][1] <- ''
    expect_error(select(bad), 'labels')
    bad <- v
    bad[[nm]][2] <- bad[[nm]][1]
    expect_error(select(bad), 'labels')
  }
  for (x in list(-1, 101, NA, Inf, c(90, 99), '99',
    1 + 1i))
      expect_error(select(v, thr = x), 'percentile')
  for (x in list('bad', NA_character_, c('VIP2D', 'Original')))
    expect_error(vip_feature_sets(v,
    strategy = x), 'strategy')
  bad <- v
  bad$ncomp <- 3L
  expect_error(select(bad), 'labels|dimensions')
  expect_error(select(v, cols = 1), 'default')
  expect_error(select(v, aggregator = 'all'), 'default')
  expect_error(select(v, strip_time = FALSE), 'default')
})
test_that('all existing percentile aggregators and set combinations preserve positive scaling', {
  for (H in c(1L, 2L, 4L)) {
    v <- .selection_vips(5L, 3L, H)
    scaled <- v
    scaled$VIP3Dmodel2 <- v$M2_mean
    for (a in c('any', 'all', 'mean'))
      for (thr in c(0, 50, 99,
      100)) {
      x <- vip_feature_sets(v, thr = thr, aggregator = a, combinations = c('individual',
        'intersection', 'union'))
      y <- vip_feature_sets(scaled, thr = thr, aggregator = a, combinations = c('individual',
        'intersection', 'union'))
      expect_identical(x$sets, y$sets)
    }
  }
})
