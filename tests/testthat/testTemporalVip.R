.temporal_fixture <- function() {
  m <- matrix(c(1, 5, 2, 8, 3, 7, 4, 9,
    2, 6, 1, 10), 4, 3, dimnames = list(c("geneB",
    "geneA", "geneD", "geneC"), c("0", "-12", "-3")))
  structure(list(M2_mean = m, VIP3Dmodel2 = 99 * m, feature_names = rownames(m),
    time_names = colnames(m)), class = "tensorpls_vips")
}
test_that("temporal summaries use canonical means and preserve axis order", {
  v <- .temporal_fixture()
  m <- v$M2_mean
  selected <- c("geneC", "geneB")
  s <- temporal_vip_summary(v, selected, dispersion = "IQR")
  expect_identical(s$importance, unname(colMeans(m[selected, , drop = FALSE])))
  expect_identical(s$time, colnames(m))
  expect_identical(attr(s, "feature_names"), selected)
  expect_identical(s$n_features, rep(2L, 3))
  expect_identical(temporal_vip_summary(v)$importance, unname(colMeans(m)))
  expect_identical(attr(temporal_vip_summary(v), "feature_names"), rownames(m))
  expect_equal(temporal_vip_summary(v, selected, "median")$importance, unname(apply(m[selected,
    , drop = FALSE], 2, median)))
  expect_equal(s$q25, unname(apply(m[selected, , drop = FALSE], 2, quantile,
    .25)))
  expect_equal(s$q75, unname(apply(m[selected, , drop = FALSE], 2, quantile,
    .75)))
  expect_true(all(is.na(temporal_vip_summary(v)$q25)))
})
test_that("plot reference and ribbons use the same summaries", {
  v <- .temporal_fixture()
  f <- c("geneC", "geneB")
  p <- plot_temporal_vip(v, f, reference = "all", dispersion = "IQR")
  expect_s3_class(p, "ggplot")
  expect_identical(p$data$importance[p$data$profile == "All features"],
    temporal_vip_summary(v)$importance)
  expect_identical(p$data$importance[p$data$profile == "Selected feature signature"], temporal_vip_summary(v,
    f)$importance)
  expect_identical(p$data$time, rep(v$time_names, 2))
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$data[[1]]$ymin, temporal_vip_summary(v, dispersion = "IQR")$q25)
  expect_equal(b$data[[4]]$ymax, temporal_vip_summary(v, f, dispersion = "IQR")$q75)
  expect_identical(b$layout$panel_params[[1]]$x$get_labels(), v$time_names)
  expect_null(p$labels$caption)
  expect_match(p$labels$subtitle, "n = 2", fixed = TRUE)
  expect_identical(plot_temporal_vip(v, summary = "median")$labels$y, "Median temporal VIP")
  expect_equal(nrow(plot_temporal_vip(v, reference = "all")$data), length(v$time_names))
  expect_equal(nrow(plot_temporal_vip(v, rev(v$feature_names), reference = "all")$data),
    length(v$time_names))
})
test_that("single-feature and single-time inputs retain dimensions and draw cleanly", {
  v <- .temporal_fixture()
  for (single_time in c(FALSE, TRUE)) {
    if (single_time) {
      v$M2_mean <- v$M2_mean[, 2, drop = FALSE]
      v$time_names <- colnames(v$M2_mean)
    }
    s <- temporal_vip_summary(v, "geneA", dispersion = "IQR")
    expect_identical(s$importance, unname(as.vector(v$M2_mean["geneA", ])))
    expect_equal(s$q25, s$importance)
    expect_equal(s$q75, s$importance)
    expect_no_warning(ggplot2::ggplotGrob(plot_temporal_vip(v, "geneA", reference = "all",
      dispersion = "IQR")))
  }
})
test_that("invalid temporal inputs fail clearly", {
  v <- .temporal_fixture()
  expect_error(temporal_vip_summary(list()), "tensorpls_vips")
  for (m in list(NULL, matrix(NA_real_, 4, 3), matrix(Inf, 4, 3),
    matrix("1", 4, 3), matrix(1i, 4, 3), matrix(numeric(),
      0, 3))) {
    bad <- v
    bad$M2_mean <- m
    expect_error(temporal_vip_summary(bad), "M2_mean")
  }
  for (f in list(character(), c("geneA", "geneA"), NA_character_, "",
    1L))
      expect_error(plot_temporal_vip(v,
    f), "nonempty character vector")
  expect_error(temporal_vip_summary(v, "unknown"), "Unknown features: unknown", fixed = TRUE)
  for (axis in 1:2) {
    for (labels in list(NULL, rep("duplicate", dim(v$M2_mean)[axis]))) {
      bad <- v
      dimnames(bad$M2_mean)[axis] <- list(labels)
      expect_error(temporal_vip_summary(bad), "feature and time names")
    }
  }
  bad <- v
  bad$time_names <- rev(v$time_names)
  expect_error(temporal_vip_summary(bad), "names and order")
  bad <- v
  bad$feature_names <- rev(v$feature_names)
  expect_error(temporal_vip_summary(bad), "names and order")
})
test_that("temporal plotting invokes no fitting, VIP, selection or factor routines", {
  v <- .temporal_fixture()
  before <- serialize(v, NULL)
  forbidden <- function(...)
    stop("Unexpected computation")
  local_mocked_bindings(fit_tensorpls_pls = forbidden, .tensorpls_fit = forbidden,
    tensorpls_vips = forbidden,
    nested_cv_tensorpls = forbidden, vip_feature_sets = forbidden, feature_selection = forbidden,
    fit_tucker_imputer = forbidden, predict_tucker_imputer = forbidden, .package = "TensorPLS")
  if (requireNamespace("mixOmics", quietly = TRUE))
    local_mocked_bindings(pls = forbidden, plsda = forbidden, block.plsda = forbidden,
      vip = forbidden, .package = "mixOmics")
  for (ref in c("none", "all"))
    for (disp in c("none", "IQR")) {
    p <- plot_temporal_vip(v, c("geneD", "geneB"), reference = ref, dispersion = disp)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(serialize(v, NULL), before)
})
test_that("public fitted VIP summaries yield the direct V average", {
  set.seed(123)
  X <- array(rnorm(24 * 4 * 5), c(24, 4, 5))
  model <- fit_tensorpls_pls(X, factor(rep(c("Control", "Case"), 12)), 3)
  v <- tensorpls_vips(model)
  expected <- apply(v$V, 2, mean)
  expect_equal(temporal_vip_summary(v)$importance, unname(expected))
  expect_identical(temporal_vip_summary(v)$importance, unname(colMeans(v$M2_mean)))
})
