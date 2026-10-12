.table_fixture <- function() {
  m <- matrix(c(2, 1, 4, 2, 2, 5, 1, 2,
    1, 3, 4, 2), 4, 3, dimnames = list(c('b',
    'a', 'd', 'c'), c('0', '-12', '-3')))
  structure(list(M2_mean = m, VIP3Dmodel2 = m * 10, feature_names = rownames(m),
    time_names = colnames(m)), class = 'tensorpls_vips')
}
.table_nested <- function() {
  structure(list(selected_features_per_outer_fold = list('1' = c('b', 'a'),
    '2' = c('d', 'a'), '3' = NULL, '4' = c('c', 'a')),
    settings = list(outer_folds = 4), stability = list(frequency = data.frame(feature = c('b',
      'a', 'd', 'c'), count = 999, frequency = 999)),
        inner_tuning_results = list(selected_features = rep('b',
      100)), outer_fold_metrics = data.frame(fold = 1:4, status = c('ok',
      'failed', 'failed', 'ok'))), class = 'nested_tensorpls')
}
test_that('table values, schema and peak ties use canonical M2', {
  v <- .table_fixture()
  tab <- temporal_vip_table(v)
  m <- v$M2_mean[tab$feature, , drop = FALSE]
  expect_s3_class(tab, 'data.frame')
  expect_named(tab, c('feature', 'mean_M2', 'peak_M2', 'peak_time', 'M2_0', 'M2_-12', 'M2_-3'))
  expect_setequal(tab$feature, v$feature_names)
  expect_identical(tab$mean_M2, unname(rowMeans(m)))
  expect_identical(tab$peak_M2, unname(apply(m, 1, max)))
  expect_identical(tab$peak_time, c('-12', '0', '0', '0'))
  expect_equal(unname(as.matrix(tab[, 5:7])), unname(m), tolerance = 0)
  expect_identical(tab$feature, c('a', 'd', 'c', 'b'))
  expect_identical(temporal_vip_table(v, order_by = 'peak_M2')$feature, c('a', 'd',
    'b', 'c'))
  for (f in list(c('d', 'a'), c('a', 'd'))) {
    expect_identical(temporal_vip_table(v, f)$feature, f)
    expect_setequal(temporal_vip_table(v, f)$feature, f)
  }
  one <- v
  one$M2_mean <- v$M2_mean[, 2, drop = FALSE]
  one$time_names <- colnames(one$M2_mean)
  expect_identical(temporal_vip_table(one, 'b')$peak_time, '-12')
  expect_identical(temporal_vip_table(one, 'b')$mean_M2, 2)
})
test_that('all rows are supported and top 30 agrees with heatmap', {
  v <- .table_fixture()
  m <- matrix(seq_len(20001 * 3)%%29, 20001, 3, dimnames = list(paste0('gene',
    1:20001), v$time_names))
  v$M2_mean <- m
  v$feature_names <- rownames(m)
  expect_equal(nrow(temporal_vip_table(v)), 20001)
  f <- rev(v$feature_names[1:75])
  tab <- temporal_vip_table(v, f)
  p <- plot_temporal_vip_heatmap(v, f, top_n = 30)
  expect_identical(head(tab$feature, 30), rev(levels(p$data$feature)))
  for (t in seq_along(v$time_names))
    expect_identical(tab[[paste0('M2_', v$time_names[t])]][1:30],
    p$data$value[p$data$time == v$time_names[t]])
})
test_that('outer-refit selections define counts and eligible denominator', {
  v <- .table_fixture()
  n <- .table_nested()
  tab <- temporal_vip_table(v, nested = n)
  expect_named(tab, c('feature', 'mean_M2', 'peak_M2', 'peak_time', 'M2_0', 'M2_-12', 'M2_-3',
    'cv_selection_count', 'cv_selection_frequency'))
  expect_identical(tab$cv_selection_count, c(3L, 1L, 1L, 1L))
  expect_equal(tab$cv_selection_frequency, c(1, 1 / 3, 1 / 3,
    1 / 3))
  expect_identical(temporal_vip_table(v, nested = n,
    order_by = 'cv_selection_frequency')$feature, c('a',
    'b', 'd', 'c'))
  expect_identical(temporal_vip_table(v, c('c', 'b'), n, 'cv_selection_frequency')$feature, c('c',
    'b'))
  n$inner_tuning_results <- list(selected_features = 'c')
  n$stability$frequency$count <- 0
  expect_identical(tab, temporal_vip_table(v, nested = n))
  n$selected_features_per_outer_fold[['4']] <- 'a'
  zero <- temporal_vip_table(v, 'c', n)
  expect_identical(zero$cv_selection_count, 0L)
  expect_identical(zero$cv_selection_frequency, 0)
})
test_that('invalid feature and nested inputs fail clearly', {
  v <- .table_fixture()
  for (f in list(character(), 1, c('a', 'a'), NA_character_, ''))
    expect_error(temporal_vip_table(v,
    f), 'nonempty character vector')
  expect_error(temporal_vip_table(v, 'absent'), 'Unknown features')
  expect_error(temporal_vip_table(v, order_by = 'cv_selection_frequency'), 'requires nested')
  expect_error(temporal_vip_table(v, nested = list()), 'nested_cv_tensorpls')
  n <- .table_nested()
  n$selected_features_per_outer_fold <- NULL
  expect_error(temporal_vip_table(v, nested = n), 'selected_features_per_outer_fold')
  n <- .table_nested()
  n$selected_features_per_outer_fold <- n$selected_features_per_outer_fold[1:3]
  expect_error(temporal_vip_table(v, nested = n), 'every declared outer fold')
  n <- .table_nested()
  n$selected_features_per_outer_fold[['1']] <- c('a', 'a')
  expect_error(temporal_vip_table(v, nested = n), 'Invalid outer-refit')
  n <- .table_nested()
  n$selected_features_per_outer_fold[] <- list(NULL)
  expect_error(temporal_vip_table(v, nested = n), 'no eligible')
  n <- .table_nested()
  n$stability$frequency$feature[1] <- 'other'
  expect_error(temporal_vip_table(v, nested = n), 'absent from the nested')
})
test_that('export does not fit, select, recompute VIPs or plot', {
  v <- .table_fixture()
  n <- .table_nested()
  before <- serialize(list(v, n), NULL)
  forbidden <- function(...)
    stop('Unexpected computation')
  local_mocked_bindings(fit_tensorpls_pls = forbidden, .tensorpls_fit = forbidden,
    tensorpls_vips = forbidden,
    nested_cv_tensorpls = forbidden, vip_feature_sets = forbidden, feature_selection = forbidden,
    plot_temporal_vip = forbidden, plot_temporal_vip_heatmap = forbidden,
      plot_nplsda_scores = forbidden,
    .package = 'TensorPLS')
  local_mocked_bindings(ggplot = forbidden, .package = 'ggplot2')
  local_mocked_bindings(pls = forbidden, plsda = forbidden, block.plsda = forbidden,
    vip = forbidden, .package = 'mixOmics')
  expect_s3_class(temporal_vip_table(v), 'data.frame')
  expect_s3_class(temporal_vip_table(v, nested = n), 'data.frame')
  expect_identical(serialize(list(v, n), NULL), before)
})
