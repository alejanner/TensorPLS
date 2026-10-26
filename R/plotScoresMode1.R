#' Plot TensorPLS scores
#'
#' @param scores_matrix A tensorpls_model from [fit_tensorpls_pls()], or a
#'   matrix of scores for the matrix-input interface.
#' @param expl_var_matrix Matrix with variance explained (with R2.Y column)
#' @param class_vec Factor with class labels. For model input, optional display
#'   labels with the same subject assignments and level order as model classes.
#' @param pc1 First component to plot (default = 1)
#' @param pc2 Second component to plot (default = 2)
#' @param comp Two distinct positive fitted component indices for model input.
#' @param variance Whether matrix-input axes display supplied variance values.
#'   Model input requires FALSE and uses plain component labels.
#' @param main_title Plot title (default = "TensorPLS Score Plot")
#'
#' @details For model input, coordinates are taken directly from x_scores.
#'   Subject classes come directly from model$training_class, a factor named
#'   by subject_ids. Its class-level order is preserved. Model input requires
#'   this field; missing training classes produce an error.
#'   No model fitting, score transformation or scaling is performed.
#'   This is descriptive visualization, not predictive validation.
#'   Model input accepts comp, main_title and optional class_vec.
#'   Points have class-dependent fills, black borders, and dashed zero-reference
#'   lines. These lines are not classification boundaries.
#' @return Draws a base-R score plot. For model input, invisibly returns a
#'   data frame with subject_id, x, y and class in model subject order.
#'   For matrix input, invisibly returns NULL.
#' @examples
#' set.seed(123)
#' X <- array(rnorm(20 * 4 * 3), c(20, 4, 3))
#' model <- fit_tensorpls_pls(X, factor(rep(c("A", "B"), 10)), 3)
#' plot_nplsda_scores(model, comp = c(1, 3))
#' @importFrom graphics abline legend par
#' @importFrom stats setNames
#' @export
plot_nplsda_scores <- function(scores_matrix, expl_var_matrix, class_vec, pc1 = 1,
  pc2 = 2, variance = FALSE, main_title = "TensorPLS Score Plot", comp = c(1, 2)) {
  if (inherits(scores_matrix, "tensorpls_model")) {
    if (!missing(expl_var_matrix) || !missing(pc1) || !missing(pc2) || !identical(variance,
      FALSE))
      stop("For model input, supply only the model, comp, main_title and optional class_vec.", call. = FALSE)
    return(.tensorpls_sample_plot(scores_matrix, comp, main_title, if (missing(class_vec))
      NULL else
        class_vec))
  }
  if (!missing(comp))
    stop("comp requires a tensorpls_model input; use pc1/pc2 for a score matrix.", call. = FALSE)
  if (pc1 > ncol(scores_matrix) || pc2 > ncol(scores_matrix)) {
    stop(sprintf("Component indices (%d, %d) exceed available components (%d)", pc1,
      pc2, ncol(scores_matrix)))
  }
  # Extract coordinates
  x <- scores_matrix[, pc1]
  y <- scores_matrix[, pc2]
  if (length(class_vec) != length(x)) {
    stop(sprintf("Length mismatch: class_vec has %d elements, scores have %d observations", length(class_vec), length(x)))
  }
  # Axis labels with R2.Y percentages
  xlab <- sprintf("Comp %d: %.1f%%", pc1, expl_var_matrix[pc1, 3] * 100)
  ylab <- sprintf("Comp %d: %.1f%%", pc2, expl_var_matrix[pc2, 3] * 100)
  # Set colors based on classes
  class_levels <- levels(class_vec)
  n_classes <- length(class_levels)
  color_palette <- c("dodgerblue2", "orange2", "green3", "red3", "purple3", "gold2")
  class_colors <- setNames(color_palette[1:n_classes], class_levels)
  cols <- class_colors[as.character(class_vec)]
  # Plot with external legend
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4, 4, 8))
  if (variance) {
    plot(x, y, pch = 21, bg = cols, col = "black", xlab = xlab, ylab = ylab,
      main = main_title)
    abline(h = 0, v = 0, lty = 2, col = "black")
  } else {
    plot(x, y, pch = 21, bg = cols, col = "black", xlab = sprintf("Comp %d ", pc1),
      ylab = sprintf("Comp %d ", pc2), main = main_title)
    abline(h = 0, v = 0, lty = 2, col = "black")
  }
  # External legend - set xpd=TRUE only for legend
  par(xpd = TRUE)
  legend(x = par("usr")[2] + diff(par("usr")[1:2]) * 0.02,
    y = mean(par("usr")[3:4]), legend = names(class_colors), pch = 21,
    pt.bg = class_colors, col = "black", bty = "n", title = "Classes")
  invisible(NULL)
}
.tensorpls_sample_plot <- function(model, comp, main_title, display_class = NULL) {
  H <- model$ncomp
  T <- model$x_scores
  ids <- model$subject_ids
  if (!is.numeric(H) || is.complex(H) || length(H) != 1L || !is.finite(H) || H <
    1 || H != floor(H) || !is.matrix(T) || !is.numeric(T) || is.complex(T) ||
    !nrow(T) || ncol(T) != H || any(!is.finite(T)))
    stop("Model must contain finite native x_scores matching ncomp.", call. = FALSE)
  if (!is.numeric(comp) || is.complex(comp) || length(comp) != 2L || any(!is.finite(comp)) ||
    any(comp != floor(comp)) || any(comp < 1 | comp > H) || anyDuplicated(comp))
    stop("comp must contain two distinct positive integer indices <= model$ncomp.", call. = FALSE)
  if (!is.character(ids) || length(ids) != nrow(T) || anyNA(ids) || any(!nzchar(ids)) ||
    anyDuplicated(ids))
    stop("Model subject_ids must be unique and match score row order.", call. = FALSE)
  classes <- model$training_class
  if (!is.factor(classes) || length(classes) != nrow(T) || anyNA(classes) || !identical(names(classes),
    ids) || nlevels(classes) < 2L || anyNA(levels(classes)) || any(!nzchar(levels(classes))) ||
    anyDuplicated(levels(classes)))
    stop("Model training_class must be a valid factor named in subject_ids order.", call. = FALSE)
  d <- data.frame(subject_id = ids, x = T[, comp[1]], y = T[, comp[2]],
    class = unname(classes))
  if (!is.null(display_class)) {
    if (!is.factor(display_class) || length(display_class) != nrow(d) || anyNA(display_class) ||
      nlevels(display_class) != nlevels(classes) || anyNA(levels(display_class)) ||
      any(!nzchar(levels(display_class))) || !identical(as.integer(display_class),
      as.integer(d$class)))
      stop("class_vec must preserve model class assignments and factor level order.", call. = FALSE)
    d$class <- display_class
  }
  plot_nplsda_scores(T, matrix(NA_real_, H, 3L), d$class, pc1 = comp[1], pc2 = comp[2],
    variance = FALSE, main_title = main_title)
  invisible(d)
}
