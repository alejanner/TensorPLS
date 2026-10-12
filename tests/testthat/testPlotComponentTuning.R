.tuning_plot_fixture <- function(metric = 'Q2Y') {
  structure(list(tuning = data.frame(ncomp = c(3L, 1L, 2L), selection_score = c(.6,
    .2, .5)), selected_ncomp = 3L, settings = list(selection_metric = metric),
    fold_results = data.frame(ncomp = rep(c(3L, 1L, 2L), 2),
      inner_fold = rep(1:2, each = 3), Q2Y = c(.3, -.2,
        .1, .4, -.1, .2), BER = c(.1, .4, .2,
        .2, .5, .3))), class = 'tensorpls_component_tuning')
}
test_that('pooled Q2Y is stored data in H order and plotting performs no computation', {
  withr::local_pdf(tempfile(fileext = '.pdf'))
  x <- .tuning_plot_fixture()
  before <- serialize(x, NULL)
  forbidden <- function(...)
    stop('unexpected computation')
  local_mocked_bindings(tune_tensorpls_components = forbidden, .nested_develop = forbidden,
    .nested_fit_pipeline = forbidden, .nested_predict_pipeline = forbidden,
      .nested_metrics = forbidden,
    .tensorpls_fit = forbidden, fit_tensorpls_pls = forbidden, fit_tucker_imputer = forbidden,
    predict_tucker_imputer = forbidden, tensorpls_vips = forbidden, vip_feature_sets = forbidden,
    .package = 'TensorPLS')
  z <- withVisible(plot(x, show_folds = FALSE))
  expect_false(z$visible)
  expect_identical(z$value$ncomp, 1:3)
  expect_identical(z$value$pooled_score, c(.2, .5, .6))
  expect_identical(z$value$selected, c(FALSE, FALSE, TRUE))
  expect_identical(z$value$metric, rep('Q2Y', 3))
  expect_identical(serialize(x, NULL), before)
})
test_that('BER uses positive display scale and keeps stored selection', {
  withr::local_pdf(tempfile(fileext = '.pdf'))
  x <- .tuning_plot_fixture('BER')
  x$tuning$selection_score <- c(-.1, -.4, -.2)
  z <- plot(x, show_folds = FALSE)
  expect_identical(z$pooled_score, c(.4, .2, .1))
  expect_identical(z$metric, rep('BER', 3))
  expect_identical(z$ncomp[z$selected], 3L)
})
test_that('fold trajectories use stored metrics without replacing pooled scores', {
  withr::local_pdf(tempfile(fileext = '.pdf'))
  calls <- list()
  original <- graphics::lines
  local_mocked_bindings(lines = function(x, y, ...) {
    calls[[length(calls) + 1L]] <<- list(x = x, y = y)
    original(x, y, ...)
  }, .package = 'graphics')
  for (metric in c('Q2Y', 'BER')) {
    calls <- list()
    x <- .tuning_plot_fixture(metric)
    if (metric == 'BER')
      x$tuning$selection_score <- -x$tuning$selection_score
    z <- plot(x, show_folds = TRUE)
    expect_length(calls, 3)
    for (f in 1:2) {
      a <- x$fold_results[x$fold_results$inner_fold == f, ]
      a <- a[order(a$ncomp), ]
      expect_identical(calls[[f]]$y, a[[metric]])
      expect_identical(calls[[f]]$x, a$ncomp)
    }
    expect_identical(calls[[3]]$y, z$pooled_score)
    expect_identical(z, plot(x, show_folds = FALSE))
  }
})
test_that('missing selection and nonfinite scores do not invent a component', {
  withr::local_pdf(tempfile(fileext = '.pdf'))
  x <- .tuning_plot_fixture()
  x$selected_ncomp <- NA_integer_
  local_mocked_bindings(abline = function(...)
    stop('unexpected marker'), .package = 'graphics')
  z <- plot(x)
  expect_false(any(z$selected))
  expect_identical(z$pooled_score, c(.2, .5, .6))
  x$tuning$selection_score[2] <- NA_real_
  expect_identical(plot(x, show_folds = FALSE)$pooled_score, c(NA_real_, .5, .6))
  x$tuning$selection_score[] <- NA_real_
  expect_true(all(is.na(plot(x, show_folds = FALSE)$pooled_score)))
  x$selected_ncomp <- 3L
  x$tuning$selection_score <- c(.6, .2, .5)
  expect_true(plot(x, mark_selected = FALSE)$selected[3])
})
test_that('malformed plot inputs fail clearly', {
  method <- getS3method('plot', 'tensorpls_component_tuning')
  expect_error(method(list()), 'tensorpls_component_tuning')
  x <- .tuning_plot_fixture()
  x$settings$selection_metric <- 'AUC'
  expect_error(plot(x), 'selection metric')
  x <- .tuning_plot_fixture()
  x$tuning$selection_score <- NULL
  expect_error(plot(x), 'tuning must')
  for (h in list(c(1, NA, 3), c(1, Inf, 3), c(1,
    1, 3), c(1, 2.5, 3))) {
    x <- .tuning_plot_fixture()
    x$tuning$ncomp <- h
    expect_error(plot(x), 'tuning must')
  }
  x <- .tuning_plot_fixture()
  x$fold_results$Q2Y <- NULL
  expect_error(plot(x, show_folds = TRUE), 'fold_results')
})
test_that('default plot draws only all ten pooled candidates without annotations', {
  withr::local_pdf(tempfile(fileext = '.pdf'))
  forbidden <- function(...)
    stop('unexpected computation or annotation')
  local_mocked_bindings(tune_tensorpls_components = forbidden, .nested_develop = forbidden,
    .nested_fit_pipeline = forbidden, .nested_predict_pipeline = forbidden,
      .nested_metrics = forbidden,
    .tensorpls_fit = forbidden, fit_tensorpls_pls = forbidden, .nested_imputer = forbidden,
    fit_tucker_imputer = forbidden, predict_tucker_imputer = forbidden, tensorpls_vips = forbidden,
    vip_feature_sets = forbidden, feature_selection = forbidden, .package = 'TensorPLS')
  lines_seen <- list()
  ticks <- NULL
  args <- NULL
  original_lines <- graphics::lines
  original_axis <- graphics::axis
  original_modify <- utils::modifyList
  local_mocked_bindings(modifyList = function(x, val, ...) {
    args <<- x
    original_modify(x, val, ...)
  }, .package = "utils")
  local_mocked_bindings(legend = forbidden, text = forbidden, mtext = forbidden, lines = function(x,
    y, ...) {
    lines_seen[[length(lines_seen) + 1L]] <<- list(x = x, y = y)
    original_lines(x, y, ...)
  }, axis = function(side, at, labels, ...) {
    ticks <<- at
    original_axis(side, at, labels, ...)
  }, .package = 'graphics')
  for (metric in c('Q2Y', 'BER')) {
    x <- .tuning_plot_fixture(metric)
    x$tuning <- data.frame(ncomp = 10:1, selection_score = seq(.1,
      1, length.out = 10))
    if (metric == 'BER')
      x$tuning$selection_score <- -x$tuning$selection_score
    x$selected_ncomp <- 7L
    before <- serialize(x, NULL)
    lines_seen <- list()
    z <- plot(x)
    expect_identical(z$ncomp, 1:10)
    expect_identical(ticks, 1:10)
    expect_length(lines_seen, 1)
    expect_identical(lines_seen[[1]]$x, 1:10)
    expect_identical(lines_seen[[1]]$y, z$pooled_score)
    expect_identical(z$pooled_score, rev(if (metric == 'BER')
      -x$tuning$selection_score else
        x$tuning$selection_score))
    expect_identical(z$ncomp[z$selected], 7L)
    expect_identical(serialize(x, NULL), before)
    expect_identical(args$xlab, 'Number of components')
    expect_identical(args$ylab, metric)
    expect_identical(args$main, 'Component tuning')
    expect_null(args$sub)
  }
})
