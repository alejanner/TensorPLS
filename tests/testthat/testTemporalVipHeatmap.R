.map_fixture <- function() {
  m <- matrix(c(1, 4, 2, 4, 3, 2, 1, 2,
    2, 3, 6, 3), 4, 3, dimnames = list(c('b',
    'a', 'd', 'c'), c('0', '-12', '-3')))
  structure(list(M2_mean = m, feature_names = rownames(m), time_names = colnames(m)),
    class = 'tensorpls_vips')
}
test_that('tiles preserve values and time order with deterministic ranked rows', {
  v <- .map_fixture()
  f <- c('c', 'd', 'a')
  p <- plot_temporal_vip_heatmap(v, f)
  expect_s3_class(p, 'ggplot')
  expect_null(p$labels$caption)
  expect_identical(levels(p$data$time), v$time_names)
  rows <- f[order(-rowMeans(v$M2_mean[f, , drop = FALSE]), seq_along(f))]
  expect_identical(rev(levels(p$data$feature)), rows)
  for (i in seq_len(nrow(p$data)))
    expect_identical(p$data$value[i],
    unname(v$M2_mean[as.character(p$data$feature[i]),
    as.character(p$data$time[i])]))
  b <- ggplot2::ggplot_build(p)
  expect_identical(b$layout$panel_params[[1]]$x$get_labels(), v$time_names)
  expect_identical(b$layout$panel_params[[1]]$y$get_labels(), rev(rows))
  expect_equal(nrow(b$data[[1]]), length(f) * length(v$time_names))
  p2 <- plot_temporal_vip_heatmap(v, f, top_n = 2)
  expect_identical(rev(levels(p2$data$feature)), rows[1:2])
  expect_identical(p2$data, plot_temporal_vip_heatmap(v, f, top_n = 2)$data)
  expect_identical(rev(levels(plot_temporal_vip_heatmap(v, c('c', 'a'))$data$feature)),
    c('c', 'a'))
  expect_length(p$layers, 1L)
  expect_no_warning(ggplot2::ggplotGrob(plot_temporal_vip_heatmap(v, 'a')))
})
test_that('external adjusted p-values are aligned and thresholded strictly', {
  v <- .map_fixture()
  q <- v$M2_mean
  q[] <- .7
  q['d', '-12'] <- .01
  q['a', '0'] <- .05
  q['c', '-3'] <- .049
  p <- plot_temporal_vip_heatmap(v, c('c', 'd', 'a'), padj = q)
  expect_length(p$layers, 2L)
  marked <- p$layers[[2]]$data
  expect_equal(nrow(marked), 2L)
  expect_setequal(paste(marked$feature, marked$time), c('d -12', 'c -3'))
  expect_equal(marked$padj, c(.01, .049))
  expect_identical(p$data$value, plot_temporal_vip_heatmap(v, c('c', 'd', 'a'))$data$value)
  expect_equal(nrow(plot_temporal_vip_heatmap(v, top_n = 1, padj = q)$layers[[2]]$data),
    0L)
  expect_null(p$labels$caption)
  for (bad in list(q[-1, , drop = FALSE], q[rev(rownames(q)), , drop = FALSE],
    q[, rev(colnames(q)), drop = FALSE], unname(q)))
      expect_error(plot_temporal_vip_heatmap(v,
    padj = bad), 'exactly match')
  for (value in list(NA_real_, Inf, -.1, 1.1, 1i, 'bad')) {
    bad <- q
    bad[1, 1] <- value
    expect_error(plot_temporal_vip_heatmap(v, padj = bad), 'finite numeric matrix')
  }
  for (a in list(NA_real_, Inf, -1, 2, c(.1, .2), 'bad'))
    expect_error(plot_temporal_vip_heatmap(v,
    alpha = a), 'alpha')
})
test_that('large sets require manageable top_n and invalid selections error', {
  v <- .map_fixture()
  for (n in list(0, -1, 1.2, NA_real_, Inf, c(1, 2),
    '2'))
      expect_error(plot_temporal_vip_heatmap(v, top_n = n), 'positive integer')
  expect_error(plot_temporal_vip_heatmap(v, c('a', 'a')), 'unique')
  expect_error(plot_temporal_vip_heatmap(v, character()), 'nonempty')
  expect_error(plot_temporal_vip_heatmap(v, 'unknown'), 'Unknown')
  m <- matrix(1, 51, 1, dimnames = list(paste0('f', 1:51),
    'only'))
  v$M2_mean <- m
  v$feature_names <- rownames(m)
  v$time_names <- colnames(m)
  expect_error(plot_temporal_vip_heatmap(v), 'top_n <= 50')
  expect_equal(nrow(plot_temporal_vip_heatmap(v, top_n = 20)$data), 20L)
  expect_no_warning(ggplot2::ggplotGrob(plot_temporal_vip_heatmap(v, top_n = 1)))
})
test_that('profiles and heatmaps do not fit, recompute VIP or invoke permutation logic', {
  v <- .map_fixture()
  before <- serialize(v, NULL)
  forbidden <- function(...)
    stop('Unexpected computation')
  local_mocked_bindings(fit_tensorpls_pls = forbidden, .tensorpls_fit = forbidden,
    tensorpls_vips = forbidden,
    vip_feature_sets = forbidden, feature_selection = forbidden, nested_cv_tensorpls = forbidden,
    .package = 'TensorPLS')
  local_mocked_bindings(pls = forbidden, plsda = forbidden, block.plsda = forbidden,
    vip = forbidden, .package = 'mixOmics')
  q <- v$M2_mean
  q[] <- .01
  expect_no_warning(ggplot2::ggplotGrob(plot_temporal_vip_heatmap(v, padj = q)))
  for (d in c('none', 'IQR'))
    for (r in c('none', 'all')) {
    p <- plot_temporal_vip(v, c('a', 'd'), reference = r, dispersion = d,
      caption = 'omitted')
    expect_null(p$labels$caption)
    expect_identical(p$data$importance[1:3], unname(colMeans(v$M2_mean[c('a',
      'd'), ])))
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(serialize(v, NULL), before)
})
