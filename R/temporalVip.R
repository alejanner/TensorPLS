#' Summarize temporal VIP importance across features
#'
#' @param vips A `tensorpls_vips` object returned by [tensorpls_vips()].
#' @param features Unique feature names in the desired order, or NULL for all.
#' @param summary Aggregate features by their mean or median.
#' @param dispersion Include no dispersion or descriptive feature-wise IQR.
#' @details Uses `M2_mean`, the mean VIP across the cumulative component path.
#'   For the mean summary, `TI[t]` = mean over selected features of `M2_mean[f,t]`.
#'   Time order is preserved. IQR uses type-7 quantiles across features;
#'   it is not a confidence interval. This is descriptive fitted-model
#'   importance, not predictive performance or a Q2 contribution.
#' @return A data frame in model time order with time, importance, q25, q75,
#'   n_features, summary and dispersion. Quantiles are NA when not requested.
#'   The feature_names attribute records the features used, in supplied order.
#' @export
temporal_vip_summary <- function(vips, features = NULL, summary = c("mean", "median"),
  dispersion = c("none", "IQR")) {
  summary <- match.arg(summary)
  dispersion <- match.arg(dispersion)
  if (!inherits(vips, "tensorpls_vips") || !is.list(vips))
    stop("vips must be an object returned by tensorpls_vips().", call. = FALSE)
  m <- vips[["M2_mean", exact = TRUE]]
  if (!is.matrix(m) || !is.numeric(m) || is.complex(m) || any(dim(m) == 0L) ||
    any(!is.finite(m)))
    stop("vips$M2_mean must be a nonempty finite numeric matrix.", call. = FALSE)
  for (axis in seq_len(2L)) {
    labels <- dimnames(m)[[axis]]
    if (is.null(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
      stop("M2_mean must have unique, nonempty feature and time names.", call. = FALSE)
  }
  if (!identical(rownames(m), vips$feature_names) || !identical(colnames(m), vips$time_names))
    stop("M2_mean names and order must match vips feature_names and time_names.", call. = FALSE)
  if (is.null(features))
    features <- rownames(m)
  if (!is.character(features) || !length(features) || anyNA(features) || any(!nzchar(features)) ||
    anyDuplicated(features))
    stop("features must be a nonempty character vector of unique feature names.", call. = FALSE)
  missing <- setdiff(features, rownames(m))
  if (length(missing))
    stop("Unknown features: ", paste(missing, collapse = ", "), call. = FALSE)
  m <- m[features, , drop = FALSE]
  values <- if (summary == "mean")
    colMeans(m) else
      apply(m, 2L, stats::median)
  q25 <- q75 <- rep(NA_real_, ncol(m))
  if (dispersion == "IQR") {
    q25 <- apply(m, 2L, stats::quantile, probs = .25, names = FALSE, type = 7)
    q75 <- apply(m, 2L, stats::quantile, probs = .75, names = FALSE, type = 7)
  }
  result <- data.frame(time = colnames(m), importance = unname(values), q25 = unname(q25),
    q75 = unname(q75), n_features = nrow(m), summary = summary, dispersion = dispersion,
    stringsAsFactors = FALSE)
  attr(result, "feature_names") <- features
  result
}
#' Plot a descriptive temporal VIP profile
#'
#' @inheritParams temporal_vip_summary
#' @param reference Optionally overlay the all-feature profile.
#' @param ... Additional labels passed to [ggplot2::labs()], such as title. Captions are omitted.
#' @details Returns a line-and-point plot in model time order, with equally
#'   spaced time positions. A requested all-feature reference is omitted when
#'   the primary profile already contains all features. The IQR ribbon shows
#'   feature-to-feature dispersion, not uncertainty of the mean or median.
#'   No fitting, VIP recomputation or feature selection is performed.
#' @return A ggplot object. Its data contains the plotted temporal summaries,
#'   profile labels and time_index in model order.
#' @examples
#' set.seed(123)
#' X <- array(rnorm(24 * 4 * 5), c(24, 4, 5))
#' model <- fit_tensorpls_pls(X, factor(rep(c("Control", "Case"), 12)), 3)
#' vips <- tensorpls_vips(model)
#' plot_temporal_vip(vips, features = vips$feature_names[1:2], reference = "all")
#' @export
plot_temporal_vip <- function(vips, features = NULL, reference = c("none", "all"),
  summary = c("mean", "median"), dispersion = c("none", "IQR"), ...) {
  reference <- match.arg(reference)
  summary <- match.arg(summary)
  dispersion <- match.arg(dispersion)
  primary <- temporal_vip_summary(vips, features, summary, dispersion)
  all_features <- setequal(attr(primary, "feature_names"), vips$feature_names)
  label <- if (all_features)
    "All features" else
      "Selected feature signature"
  primary$profile <- label
  primary$time_index <- seq_len(nrow(primary))
  ref <- NULL
  if (reference == "all" && !all_features) {
    ref <- temporal_vip_summary(vips, NULL, summary, dispersion)
    ref$profile <- "All features"
    ref$time_index <- seq_len(nrow(ref))
  }
  data <- if (is.null(ref))
    primary else
      rbind(primary, ref)
  data$profile <- factor(data$profile, levels = unique(c(label, "All features")))
  colors <- stats::setNames(if (is.null(ref))
    "#176B87" else
      c("#176B87", "#9299A1"), if (is.null(ref))
    label else
      c(label, "All features"))
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$time_index, y = .data$importance,
    color = .data$profile, group = .data$profile))
  for (d in list(ref, primary)) {
    if (is.null(d))
      next
    is_primary <- identical(d$profile[1], label)
    if (dispersion == "IQR")
      p <- p + ggplot2::geom_ribbon(data = d, ggplot2::aes(ymin = .data$q25, ymax = .data$q75,
        fill = .data$profile), alpha = if (is_primary)
        .16 else
          .08, colour = NA, show.legend = FALSE)
    if (nrow(d) > 1L)
      p <- p + ggplot2::geom_line(data = d, linewidth = if (is_primary)
        .9 else
          .6, linetype = if (is_primary)
        "solid" else
          "dashed")
    p <- p + ggplot2::geom_point(data = d, size = if (is_primary)
      2.8 else
        1.8)
  }
  if (dispersion == "IQR")
    p <- p + ggplot2::scale_fill_manual(values = colors)
  p + ggplot2::scale_color_manual(values = colors,
    breaks = names(colors)) + ggplot2::scale_x_continuous(breaks = primary$time_index,
    labels = primary$time,
      expand = ggplot2::expansion(mult = .06)) + ggplot2::labs(title = "Temporal VIP profile",
    subtitle = paste0(label, " (n = ", format(primary$n_features[1], big.mark = ",",
      trim = TRUE), ")"), x = "Time point", y = paste0(if (summary == "mean")
      "Mean" else
        "Median", " temporal VIP"), color = NULL,
        fill = NULL) + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 12, 12, 10),
      panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"), legend.position = if (is.null(ref))
        "none" else
          "bottom") + ggplot2::labs(...) + ggplot2::labs(caption = NULL)
}
