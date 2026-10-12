#' Plot stored development component-tuning scores
#' @param x A tensorpls_component_tuning object.
#' @param show_folds Optionally overlay stored fold-specific diagnostic trajectories.
#' @param mark_selected Mark the stored selected component, when available.
#' @param main Optional title.
#' @param ... Additional base plot arguments, such as ylim or xlab.
#' @details The pooled curve uses tuning$selection_score directly. Q2Y is
#'   unchanged; negative BER scores are negated for display. Fold curves are
#'   not averaged. Scores are development tuning diagnostics, not unbiased
#'   predictive performance. No component selection or model fitting occurs.
#'   Nonfinite scores leave gaps. With no finite scores an empty frame is drawn.
#' @return Invisibly, a data.frame with ncomp, metric, pooled_score (display
#'   scale), and selected, ordered by H. Selection reflects the stored value
#'   even when mark_selected is FALSE. The input is unchanged.
#' @export
plot.tensorpls_component_tuning <- function(x, show_folds = FALSE, mark_selected = TRUE,
  main = NULL, ...) {
  if (!inherits(x, "tensorpls_component_tuning") || !is.list(x))
    stop("x must be a tensorpls_component_tuning object.", call. = FALSE)
  .vip_flag(show_folds, "show_folds")
  .vip_flag(mark_selected, "mark_selected")
  tab <- x$tuning
  metric <- x$settings$selection_metric
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) || !metric %in%
    c("Q2Y", "BER"))
    stop("Unsupported selection metric; use Q2Y or BER.", call. = FALSE)
  if (!is.data.frame(tab) || !nrow(tab) || !all(c("ncomp", "selection_score") %in% names(tab)) ||
    !is.numeric(tab$ncomp) || is.complex(tab$ncomp) || any(!is.finite(tab$ncomp)) ||
    any(tab$ncomp < 1 | tab$ncomp != floor(tab$ncomp)) || anyDuplicated(tab$ncomp) ||
    !is.numeric(tab$selection_score) || is.complex(tab$selection_score))
    stop("tuning must contain unique finite integer ncomp and numeric selection_score.", call. = FALSE)
  selected <- x$selected_ncomp
  if (!is.numeric(selected) || is.complex(selected) || length(selected) != 1L)
    stop("selected_ncomp must be one numeric value or NA.", call. = FALSE)
  d <- data.frame(ncomp = tab$ncomp, metric = metric, pooled_score = if (metric == "BER")
    -tab$selection_score else
      tab$selection_score,
      selected = is.finite(selected) & tab$ncomp == selected)
  d <- d[order(d$ncomp), , drop = FALSE]
  rownames(d) <- NULL
  folds <- NULL
  if (show_folds && !is.null(x$fold_results)) {
    folds <- x$fold_results
    if (!is.data.frame(folds) || !all(c("ncomp", "inner_fold", metric) %in% names(folds)) ||
      !is.numeric(folds$ncomp) || any(!is.finite(folds$ncomp)) || !all(folds$ncomp %in%
      d$ncomp) || anyNA(folds$inner_fold) || !is.numeric(folds[[metric]]) || is.complex(folds[[metric]]))
      stop("fold_results must contain ncomp, inner_fold and numeric metric values.", call. = FALSE)
  }
  values <- c(d$pooled_score, if (!is.null(folds))
    folds[[metric]])
  finite <- values[is.finite(values)]
  args <- utils::modifyList(list(x = d$ncomp, y = d$pooled_score, type = "n", xaxt = "n",
    ylim = if (length(finite))
      range(finite) else
        c(0, 1), xlab = "Number of components",
    ylab = metric, main = if (is.null(main))
      "Component tuning" else
        main), list(...))
  do.call(graphics::plot, args)
  graphics::axis(1, at = d$ncomp, labels = d$ncomp)
  if (!is.null(folds))
    for (f in unique(folds$inner_fold)) {
      a <- folds[folds$inner_fold == f, , drop = FALSE]
      a <- a[order(a$ncomp), , drop = FALSE]
      y <- a[[metric]]
      y[!is.finite(y)] <- NA_real_
      graphics::lines(a$ncomp, y, col = "grey75", lwd = 1)
    }
  y <- d$pooled_score
  y[!is.finite(y)] <- NA_real_
  graphics::lines(d$ncomp, y, type = "b", pch = 21, bg = "#287D8E", col = "black",
    lwd = 2)
  idx <- which(d$selected & is.finite(d$pooled_score))
  if (mark_selected && length(idx)) {
    graphics::abline(v = d$ncomp[idx], lty = 2, col = "#287D8E")
    graphics::points(d$ncomp[idx], d$pooled_score[idx], pch = 21, bg = "#E69F00",
      cex = 1.3)
  }
  invisible(d)
}
