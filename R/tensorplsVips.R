#' Expose cumulative VIP representations from a fitted TensorPLS model
#'
#' Reshape and average the base cumulative VIP values already stored by
#' [fit_tensorpls_pls()]. This function never fits or refits a PLS model.
#'
#' @param model An unmodified \code{tensorpls_model} returned by
#'   [fit_tensorpls_pls()].
#' @details For p features, k times, and H fitted components, the existing
#'   feature-fast unfolding is preserved exactly:
#'   V[f, t, h] = VIP2D[f + p * (t - 1), h]. The h axis indexes the
#'   cumulative PLS VIP path, not separate component contributions.
#'
#'   These are TensorPLS structured/derived summaries of the base cumulative
#'   VIP representation. VIP3Dmodel1[f, h] averages over ALL k time points,
#'   independently of H. VIP3Dmodel2[f, t] is the historical global sum over
#'   all H cumulative VIP columns; M2_mean[f, t] is their arithmetic mean.
#'   Thus VIP3Dmodel2 = H * M2_mean up to floating-point precision for a fixed
#'   fitted model. Neither is a sum of independent component contributions.
#'   Column h is the VIP for the cumulative h-component model.
#'
#'   No fitting, feature selection, thresholding, normalization of underlying
#'   VIP values, or recomputation of the base VIP formula is performed.
#'   M2_mean provides an explicitly mean-scaled alternative for comparison.
#'   The former duplicate M1/M2 fields are removed: use VIP3Dmodel1/M2_mean.
#'   Nested CV uses these same sum and mean field meanings.
#'
#'   VIP2D preserves the base matrix values, with feature_time row labels and
#'   Comp1, Comp2, etc. column labels. V uses feature, time, and Comp labels.
#'   VIP3Dmodel1 uses Comp1_meanTime, Comp2_meanTime, etc., consistently with the
#'   all-time summaries; VIP3Dmodel2 and M2_mean use time labels. Singleton
#'   dimensions are retained. The supplied model is not modified.
#'
#' @return A list of class \code{tensorpls_vips}, containing:
#' \describe{
#'   \item{V}{Numeric array, p x k x H, with feature/time/component dimnames.}
#'   \item{VIP2D}{The original (p * k) x H base cumulative VIP values, labelled.}
#'   \item{VIP3Dmodel1}{Numeric matrix, p x H, averaging all k times.}
#'   \item{VIP3Dmodel2}{Numeric matrix, p x k, summing all H cumulative VIP columns.}
#'   \item{M2_mean}{Numeric matrix, p x k, averaging the same H columns.}
#'   \item{feature_names, time_names}{Axis labels from model metadata.}
#'   \item{component_names}{Comp1 through CompH (cumulative-component indices).}
#'   \item{ncomp}{The fixed fitted component count H.}
#' }
#' @examples
#' set.seed(123)
#' X <- array(rnorm(24 * 4 * 5), c(24, 4, 5))
#' model <- fit_tensorpls_pls(X, factor(rep(c("Control", "Case"), 12)), 3)
#' vips <- tensorpls_vips(model)
#' dim(vips$V)
#' vips$VIP3Dmodel1
#' @export
tensorpls_vips <- function(model) {
  if (!inherits(model, "tensorpls_model") || !is.list(model))
    stop("model must be a tensorpls_model produced by fit_tensorpls_pls().", call. = FALSE)
  features <- model$feature_names
  times <- model$time_names
  for (labels in list(features, times)) {
    if (!is.character(labels) || !length(labels) || anyNA(labels) || any(!nzchar(labels)) ||
      anyDuplicated(labels))
      stop("model must contain valid unique feature_names and time_names.", call. = FALSE)
  }
  h <- model$ncomp
  if (!is.numeric(h) || is.complex(h) || length(h) != 1L || is.na(h) || !is.finite(h) ||
    h < 1 || h > .Machine$integer.max || h != floor(h))
    stop("model$ncomp must be one positive integer.", call. = FALSE)
  if (!is.list(model$core_fit))
    stop("model must contain core_fit from fit_tensorpls_pls().", call. = FALSE)
  h <- as.integer(h)
  if (!identical(as.numeric(model$core_fit$ncomp), as.numeric(h)))
    stop("model$ncomp and model$core_fit$ncomp must agree.", call. = FALSE)
  vip <- model$core_fit$VIP
  p <- length(features)
  k <- length(times)
  if (!is.matrix(vip) || !is.numeric(vip) || is.complex(vip) || !all(dim(vip) == c(p *
    k, h)) || any(!is.finite(vip)))
    stop("model$core_fit$VIP must be a finite numeric (features * times) by ncomp matrix.", call. = FALSE)
  components <- paste0("Comp", seq_len(h))
  dimnames(vip) <- list(as.vector(outer(features, times, paste, sep = "_")), components)
  V <- array(vip, c(p, k, h), dimnames = list(features, times, components))
  M1 <- matrix(apply(V, c(1, 3), mean), p, h, dimnames = list(features,
    paste0(components, "_meanTime")))
  M2_sum <- matrix(apply(V, c(1, 2), sum), p, k, dimnames = list(features,
    times))
  M2_mean <- matrix(apply(V, c(1, 2), mean), p, k, dimnames = list(features,
    times))
  structure(list(V = V, VIP2D = vip, VIP3Dmodel1 = M1, VIP3Dmodel2 = M2_sum, M2_mean = M2_mean,
    feature_names = features, time_names = times, component_names = components, ncomp = h),
    class = "tensorpls_vips")
}
