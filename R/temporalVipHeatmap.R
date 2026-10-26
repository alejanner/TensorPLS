#' Plot a temporal VIP heatmap
#'
#' @param vips A `tensorpls_vips` object returned by [tensorpls_vips()].
#' @param features Unique feature names, or NULL for all features.
#' @param top_n Optional positive integer limiting displayed rows. At most
#'   50 rows may be displayed; use top_n for larger sets.
#' @param padj Optional finite numeric matrix in `[0,1]`, with exactly the same
#'   dimensions, feature/time names and order as the full `M2_mean` matrix.
#' @param alpha Marker threshold; cells with padj strictly below alpha are marked.
#' @details Fill is `M2_mean[f,t]`, the mean VIP across the cumulative component
#'   path. Rows are ranked by descending mean across time within the supplied
#'   feature set; ties retain input order. top_n keeps the first ranked rows.
#'   Time order is unchanged. No clustering or statistical calculation occurs.
#'   Optional markers denote externally supplied adjusted p-values, separate
#'   from the descriptive VIP magnitude. No VIP2D permutation statistics are used.
#' @return A ggplot object with feature, time and value columns in its data;
#'   padj is included only when supplied. No caption is produced.
#' @importFrom rlang .data
#' @importFrom utils head
#' @export
plot_temporal_vip_heatmap <- function(vips, features = NULL, top_n = NULL, padj = NULL,
  alpha = 0.05) {
  checked <- temporal_vip_summary(vips, features)
  selected <- attr(checked, "feature_names")
  m <- vips[["M2_mean", exact = TRUE]]
  if (
    !is.numeric(alpha) ||
    is.complex(alpha) ||
    length(alpha) != 1L ||
    !is.finite(alpha) ||
    alpha < 0 ||
    alpha > 1
  )
    stop("alpha must be one finite number in [0,1].", call. = FALSE)
  if (!is.null(padj)) {
    if (
      !is.matrix(padj) ||
      !is.numeric(padj) ||
      is.complex(padj) ||
      any(!is.finite(padj)) ||
      any(padj < 0 | padj > 1)
    )
      stop("padj must be a finite numeric matrix in [0,1].", call. = FALSE)
    if (!identical(dim(padj), dim(m)) || !identical(dimnames(padj), dimnames(m)))
      stop("padj dimensions, feature/time names and order must exactly match the full M2_mean matrix.", call. = FALSE)
  }
  if (
    !is.null(top_n) && (
      !is.numeric(top_n) ||
      is.complex(top_n) ||
      length(top_n) != 1L ||
      !is.finite(top_n) ||
      top_n < 1 ||
      top_n != floor(top_n)
    )
  )
    stop("top_n must be one positive integer.", call. = FALSE)
  count <- if (is.null(top_n))
    length(selected) else
      min(length(selected), top_n)
  if (count > 50L)
    stop("More than 50 features would be displayed; supply top_n <= 50 or a smaller feature set.", call. = FALSE)
  ranked <- order(-rowMeans(m[selected, , drop = FALSE]), seq_along(selected))
  rows <- selected[head(ranked, count)]
  d <- data.frame(feature = factor(rep(rows, length(vips$time_names)), levels = rev(rows)),
    time = factor(rep(vips$time_names, each = length(rows)), levels = vips$time_names),
    value = as.vector(m[rows, , drop = FALSE]))
  if (!is.null(padj))
    d$padj <- as.vector(padj[rows, , drop = FALSE])
  label <- if (is.null(features))
    "All features" else
      "Selected feature signature"
  subtitle <- paste0(label, " (n = ", length(rows), ")")
  if (length(rows) < length(selected))
    subtitle <- paste0(subtitle, " \u00b7 Top ", length(rows), " of ", length(selected))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$feature, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = .25) +
    ggplot2::scale_fill_viridis_c(name = "M2 VIP") +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::labs(title = "Temporal VIP map", subtitle = subtitle,
      x = "Time point", y = NULL, caption = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 12, 12, 10),
      panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold"))
  if (!is.null(padj))
    p <- p + ggplot2::geom_point(data = d[d$padj < alpha, , drop = FALSE], shape = 21,
      fill = "white", colour = "black", stroke = .4, size = 1.8, show.legend = FALSE)
  p
}
