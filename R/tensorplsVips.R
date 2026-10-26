#' Extract structured VIP scores
#'
#' Extract cumulative VIP scores from a fitted TensorPLS model and compute
#' feature- and time-specific summaries without refitting.
#'
#' @param model A `tensorpls_model` returned by [fit_tensorpls_pls()].
#'
#' @details `V[f,t,h]` contains VIP scores from the cumulative
#'   h-component model. The summaries are:
#'   - `VIP3Dmodel1` (M1): mean across time for each component.
#'   - `M2_mean` (M2): mean across the cumulative component path.
#'   - `VIP3Dmodel2`: sum across the cumulative component path.
#'
#' @return A `tensorpls_vips` list containing:
#' \describe{
#'   \item{V}{VIP array (features x time x components).}
#'   \item{VIP2D}{Unfolded VIP matrix (feature-time x components).}
#'   \item{VIP3Dmodel1}{M1 matrix (features x components).}
#'   \item{M2_mean}{M2 matrix (features x time).}
#'   \item{VIP3Dmodel2}{Cumulative VIP sum (features x time).}
#'   \item{feature_names, time_names, component_names, ncomp}{Dimension metadata.}
#' }
#'
#' @examples
#' set.seed(123)
#' X <- array(rnorm(24 * 4 * 5), c(24, 4, 5))
#' model <- fit_tensorpls_pls(X, factor(rep(c("Control", "Case"), 12)), 3)
#' vips <- tensorpls_vips(model)
#' dim(vips$V)
#' vips$M2_mean
#'
#' @export
tensorpls_vips <- function(model) {
  if (!inherits(model, "tensorpls_model") || !is.list(model))
    stop("model must be a tensorpls_model produced by fit_tensorpls_pls().", call. = FALSE)
  features <- model$feature_names
  times <- model$time_names
  for (labels in list(features, times)) {
    if (
      !is.character(labels) ||
      !length(labels) ||
      anyNA(labels) ||
      any(!nzchar(labels)) ||
      anyDuplicated(labels)
    )
      stop("model must contain valid unique feature_names and time_names.", call. = FALSE)
  }
  h <- model$ncomp
  if (
    !is.numeric(h) ||
    is.complex(h) ||
    length(h) != 1L ||
    is.na(h) ||
    !is.finite(h) ||
    h < 1 ||
    h > .Machine$integer.max ||
    h != floor(h)
  )
    stop("model$ncomp must be one positive integer.", call. = FALSE)
  if (!is.list(model$core_fit))
    stop("model must contain core_fit from fit_tensorpls_pls().", call. = FALSE)
  h <- as.integer(h)
  if (!identical(as.numeric(model$core_fit$ncomp), as.numeric(h)))
    stop("model$ncomp and model$core_fit$ncomp must agree.", call. = FALSE)
  vip <- model$core_fit$VIP
  p <- length(features)
  k <- length(times)
  if (
    !is.matrix(vip) ||
    !is.numeric(vip) ||
    is.complex(vip) ||
    !all(dim(vip) == c(p * k, h)) ||
    any(!is.finite(vip))
  )
    stop("model$core_fit$VIP must be a finite numeric (features * times) by ncomp matrix.", call. = FALSE)
  components <- paste0("Comp", seq_len(h))
  dimnames(vip) <- list(as.vector(outer(features, times, paste, sep = "_")), components)
  V <- array(vip, c(p, k, h), dimnames = list(features, times, components))
  M1 <- matrix(apply(V, c(1, 3), mean), p, h, dimnames = list(features,
    paste0(components, "_meanTime")))
  M2_sum <- matrix(apply(V, c(1, 2), sum), p, k,
    dimnames = list(features, times))
  M2_mean <- matrix(apply(V, c(1, 2), mean), p, k,
    dimnames = list(features, times))
  structure(list(V = V, VIP2D = vip, VIP3Dmodel1 = M1, VIP3Dmodel2 = M2_sum, M2_mean = M2_mean,
    feature_names = features, time_names = times, component_names = components, ncomp = h),
    class = "tensorpls_vips")
}
