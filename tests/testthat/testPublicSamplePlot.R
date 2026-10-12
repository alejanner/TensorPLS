.sample_model <- function(H = 3L, multiclass = FALSE) {
  set.seed(615)
  X <- array(rnorm(24 * 4 * 3), c(24, 4, 3), dimnames = list(paste0('s',
    1:24), NULL, NULL))
  Y <- if (multiclass)
    factor(rep(c('z', 'a', 'm'), 8), levels = c('m', 'z',
      'a')) else
        factor(rep(c('z', 'a'), 12), levels = c('z', 'a'))
  list(model = fit_tensorpls_pls(X, Y, H), Y = Y)
}
test_that('model plots preserve exact native scores, class coding and order', {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (multi in c(FALSE, TRUE)) {
    d <- .sample_model(multiclass = multi)
    m <- d$model
    before <- serialize(m, NULL)
    for (comp in list(c(1, 2), c(1, 3), c(2, 3),
      c(3, 1))) {
      p <- plot_nplsda_scores(m, comp = comp)
      expect_s3_class(p, 'data.frame')
      expect_identical(p$x, unname(m$x_scores[, comp[1]]))
      expect_identical(p$y, unname(m$x_scores[, comp[2]]))
      expect_identical(p$subject_id, m$subject_ids)
      expect_identical(p$class, d$Y)
    }
    expect_identical(serialize(m, NULL), before)
  }
})
test_that('plotting never invokes fitting or selection helpers', {
  m <- .sample_model()$model
  forbidden <- function(...)
    stop('Unexpected computation')
  local_mocked_bindings(fit_tensorpls_pls = forbidden, .tensorpls_fit = forbidden,
    nested_cv_tensorpls = forbidden,
    tensorpls_vips = forbidden, fit_tucker_imputer = forbidden, predict_tucker_imputer = forbidden,
    vip_feature_sets = forbidden, .package = 'TensorPLS')
  local_mocked_bindings(pls = forbidden, plsda = forbidden, block.plsda = forbidden,
    .package = 'mixOmics')
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot_nplsda_scores(m), 'data.frame')
})
test_that('invalid components and malformed model data are rejected', {
  m <- .sample_model()$model
  for (comp in list(c(0, 1), c(-1, 2), c(1, 4),
    c(1, 1), 1, 1:3, c(1, 1.5), c(NA,
      2), c(Inf, 2), c('1', '2'), c(TRUE, FALSE)))
        expect_error(plot_nplsda_scores(m,
    comp = comp), 'two distinct')
  expect_error(plot_nplsda_scores(.sample_model(1L)$model), 'two distinct')
  expect_error(plot_nplsda_scores(m, class_vec = c('a', 'b')), 'class_vec')
  expect_error(plot_nplsda_scores(m, pc1 = 1), 'only the model')
  expect_error(plot_nplsda_scores(m, variance = TRUE), 'only the model')
  bad <- m
  bad$x_scores[1, 1] <- NA
  expect_error(plot_nplsda_scores(bad), 'x_scores')
  bad <- m
  bad$subject_ids <- rev(m$subject_ids[-1])
  expect_error(plot_nplsda_scores(bad), 'subject_ids')
  for (value in list(NULL, as.character(m$training_class), unname(m$training_class))) {
    bad <- m
    bad$training_class <- value
    expect_error(plot_nplsda_scores(bad), 'training_class')
  }
  bad <- m
  names(bad$training_class) <- rev(m$subject_ids)
  expect_error(plot_nplsda_scores(bad), 'training_class')
  bad <- m
  bad$training_class[1] <- NA
  expect_error(plot_nplsda_scores(bad), 'training_class')
})
test_that('model plotting uses the base renderer and zero-reference lines exactly', {
  m <- .sample_model()$model
  captured <- lines <- NULL
  e <- new.env(parent = environment(plot_nplsda_scores))
  e$plot_nplsda_scores <- plot_nplsda_scores
  e$.tensorpls_sample_plot <- .tensorpls_sample_plot
  environment(e$plot_nplsda_scores) <- e
  environment(e$.tensorpls_sample_plot) <- e
  e$plot <- function(x, y, ...) {
    captured <<- c(list(x = x, y = y), list(...))
  }
  e$abline <- function(...) {
    lines <<- list(...)
  }
  e$legend <- function(...)
    NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  labels <- factor(rep(c('Control', 'Case'), 12), levels = c('Control', 'Case'))
  before <- serialize(m, NULL)
  d <- e$plot_nplsda_scores(m, class_vec = labels)
  expect_identical(captured$x, m$x_scores[, 1])
  expect_identical(captured$y, m$x_scores[, 2])
  expect_identical(captured$pch, 21)
  expect_identical(captured$col, 'black')
  expect_identical(unname(captured$bg), rep(c('dodgerblue2', 'orange2'), 12))
  expect_identical(lines, list(h = 0, v = 0, lty = 2, col = 'black'))
  expect_identical(d$class, labels)
  expect_identical(serialize(m, NULL), before)
  expect_error(plot_nplsda_scores(m, class_vec = factor(labels, levels = rev(levels(labels)))),
    'class_vec')
})
test_that('plotting uses stored classes without reading predictions or residuals', {
  m <- .sample_model()$model
  poison <- new.env(parent = emptyenv())
  reads <- 0L
  for (nm in c('y.pred', 'resid'))
    makeActiveBinding(nm, function(value) {
    reads <<- reads + 1L
    stop('Training predictions/residuals must not be read')
  }, poison)
  m$core_fit <- poison
  m$training_class[] <- rev(as.character(m$training_class))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- plot_nplsda_scores(m)
  expect_identical(d$class, unname(m$training_class))
  expect_identical(reads, 0L)
  expect_identical(d$x, m$x_scores[, 1])
  expect_identical(d$y, m$x_scores[, 2])
  bad_labels <- factor(rev(as.character(m$training_class)), levels = levels(m$training_class))
  expect_error(plot_nplsda_scores(m, class_vec = bad_labels), 'class_vec')
})
