test_that('nested and public VIP names agree on the same fit without changing percentile sets', {
  set.seed(123)
  features <- paste0('gene_', c(12, 1, 11, 2, 10, 3,
    9, 4, 8, 5, 7, 6))
  times <- c('-12', '-9', '-6', '-3', '0')
  X <- array(rnorm(36 * 12 * 5), c(36, 12, 5), dimnames = list(paste0('s',
    1:36), features, times))
  Y <- factor(rep(c('Control', 'Case'), 18))
  for (H in c(1L, 2L, 3L, 5L)) {
    model <- fit_tensorpls_pls(X, Y, H)
    public <- tensorpls_vips(model)
    nested <- .nested_vips(model$core_fit, features, times)
    expect_named(nested, c('VIP2D', 'VIP3Dmodel1', 'VIP3Dmodel2', 'M2_mean'))
    for (name in names(nested))
      expect_identical(public[[name]], nested[[name]])
    expect_equal(nested$VIP3Dmodel2, H * nested$M2_mean, tolerance = 1e-14)
    mean_path <- nested
    mean_path$VIP3Dmodel2 <- nested$M2_mean
    for (thr in c(95, 99)) {
      sums <- vip_feature_sets(nested, thr = thr, sources = c("VIP2D", "VIP3Dmodel2"))
      means <- vip_feature_sets(mean_path, thr = thr, sources = c("VIP2D", "VIP3Dmodel2"))
      for (strategy in c('VIP3Dmodel2', 'VIP2D_intersect_VIP3Dmodel2')) {
        expect_identical(sums$sets[[strategy]], means$sets[[strategy]])
        selected <- features[features %in% sums$sets[[strategy]]]
        expect_identical(selected, features[features %in% means$sets[[strategy]]])
        expect_identical(selected, vip_feature_sets(public, thr = thr,
          strategy = strategy)$selected_features)
      }
      expect_identical(feature_selection(nested$VIP3Dmodel2, thr = thr),
        feature_selection(nested$M2_mean,
        thr = thr))
    }
  }
})
