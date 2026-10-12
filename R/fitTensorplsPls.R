#' Fit one TensorPLS PLS-DA model with fixed components
#'
#' Fit the same native PLS core used by [nested_cv_tensorpls()] to a complete
#' tensor. This interface fits one model only: no imputation, cross-validation,
#' component tuning, feature selection, or alternative PLS implementation.
#'
#' @param X Complete, finite, real numeric array with dimensions subjects x
#'   features x time. Missing values must be imputed before calling this function.
#' @param Y Classification response using the nested TensorPLS conventions:
#'   a factor or class-label vector, a one-column numeric 0/1 matrix, a numeric
#'   one-hot matrix, or a three-way response array with identical class coding
#'   at every time. Both binary classes / all declared one-hot classes must occur.
#'   A two-level factor encodes its second level as 1; other class-label vectors
#'   use factor ordering. More than two levels are one-hot encoded. Repeated
#'   response arrays must have the same time dimension/order as X and collapse
#'   to one response block. Continuous-response regression and time-varying
#'   response coding are not supported by this classification interface.
#' @param ncomp One fixed positive integer, at most min(n - 1, p * k), where
#'   dim(X) = c(n, p, k). The requested value is never tuned or truncated.
#'   Numerical rank exhaustion is reported by the native core.
#'
#' @details Unfolding is subjects x (features * time), with column
#'   j = f + p * (t - 1). Centering and sample-SD scaling of X and encoded Y
#'   are estimated on all supplied subjects by the native core.
#'
#'   Supplied labels are preserved. When X has no subject labels, Y subject
#'   labels are used if present; otherwise subject IDs are "1", "2", etc.
#'   Unlabelled features use Var1, Var2, etc.; unlabelled times use Y time
#'   labels if available, otherwise Time1, Time2, etc. Present labels must be
#'   unique and nonempty; shared X/Y axes must match in order. No reordering
#'   is performed. Native matrices are not relabelled: use the metadata to
#'   identify their rows/columns while preserving native scores exactly.
#'
#'   For each encoded response, fitted R2Y is 1 - SSE/SST, where SSE is the
#'   native training residual sum of squares and SST is the sum of squares
#'   about that response's training mean. fitted_R2Y is the unweighted mean
#'   of these response-specific values. It is an in-sample descriptive fit
#'   statistic, not cross-validated or predictive performance. No accuracy
#'   or AUC is calculated. A non-converged core fit is returned with a warning
#'   and its native convergence diagnostics; it is not automatically refitted.
#'
#'   The core's base cumulative VIP matrix is retained only inside core_fit.
#'   This function does not reshape or aggregate VIPs into TensorPLS summaries.
#'
#' @return A list of class \code{tensorpls_model}, containing:
#' \describe{
#'   \item{ncomp}{The fixed fitted component count (integer).}
#'   \item{feature_names, time_names, subject_ids}{Axis labels in input order.}
#'   \item{training_class}{Factor of original training class assignments in
#'     subject order, with levels from the response encoder and names equal
#'     to subject_ids.}
#'   \item{x_scores}{Native core x.scores (n x ncomp), unchanged.}
#'   \item{x_weights}{Native raw.wgs ((p * k) x ncomp), unchanged.}
#'   \item{x_loadings}{Native x.loads ((p * k) x ncomp), unchanged.}
#'   \item{coefficients}{Native reg.coefs on the original response/predictor
#'     scales, ((p * k + 1) x q). The last row is INTERCEPT.}
#'   \item{fitted_y}{Native y.pred (n x q) for the encoded response.}
#'   \item{fitted_R2Y}{Mean response-specific in-sample R2Y, a scalar.}
#'   \item{fitted_R2Y_per_response}{Named vector of in-sample R2Y values.}
#'   \item{response_encoding}{List with levels, binary, and response_names,
#'     identifying the class coding and fitted-response column order.}
#'   \item{preprocessing}{Native x_mean, x_sd, y_mean, and y_sd.}
#'   \item{convergence}{Native per-component convergence data frame.}
#'   \item{core_fit}{The native TensorPLS fit, including its base VIP.}
#' }
#' @examples
#' set.seed(123)
#' X <- array(rnorm(24 * 4 * 3), c(24, 4, 3))
#' Y <- factor(rep(c("Control", "Case"), each = 12),
#'             levels = c("Control", "Case"))
#' model <- fit_tensorpls_pls(X, Y, ncomp = 2)
#' model$fitted_R2Y # descriptive in-sample fit only
#' head(model$x_scores)
#' @export
fit_tensorpls_pls <- function(X, Y, ncomp) {
  if (!is.numeric(X) || is.complex(X) || length(dim(X)) != 3L || any(dim(X) < 1L) ||
    dim(X)[1L] < 2L)
    stop("X must be a real numeric 3D tensor with at least two subjects and nonempty feature/time axes.", call. = FALSE)
  if (anyNA(X))
    stop("X contains missing values; impute X first before calling fit_tensorpls_pls().", call. = FALSE)
  if (any(!is.finite(X)))
    stop("X must contain only finite values.", call. = FALSE)
  if (is.complex(ncomp))
    stop("ncomp must be one positive integer.", call. = FALSE)
  .nested_integer(ncomp, "ncomp")
  n <- dim(X)[1L]
  p <- dim(X)[2L]
  k <- dim(X)[3L]
  if (ncomp > min(n - 1L, p * k))
    stop("ncomp must not exceed min(n_subjects - 1, n_features * n_times); it is not truncated.", call. = FALSE)
  ncomp <- as.integer(ncomp)
  ids <- dimnames(X)[[1L]]
  features <- dimnames(X)[[2L]]
  times <- dimnames(X)[[3L]]
  yids <- if (is.null(dim(Y)))
    names(Y) else
      dimnames(Y)[[1L]]
  ytimes <- if (length(dim(Y)) == 3L)
    dimnames(Y)[[3L]] else
      NULL
  if (!is.null(ids) && !is.null(yids) && !identical(ids, yids))
    stop("X and Y subject labels must match in identical order.", call. = FALSE)
  if (length(dim(Y)) == 3L && (any(dim(Y) < 1L) || dim(Y)[3L] != k ||
    (!is.null(times) && !is.null(ytimes) && !identical(times, ytimes))))
    stop("The Y response tensor must match the X time axis in length and order.", call. = FALSE)
  if (is.null(ids))
    ids <- if (!is.null(yids))
      yids else
        as.character(seq_len(n))
  if (is.null(features))
    features <- paste0("Var", seq_len(p))
  if (is.null(times))
    times <- if (!is.null(ytimes))
      ytimes else
        paste0("Time", seq_len(k))
  for (axis in list(subject_ids = ids, feature_names = features, time_names = times)) {
    if (!is.character(axis) || anyNA(axis) || any(!nzchar(axis)) || anyDuplicated(axis))
      stop("Subject, feature and time labels must be unique, nonempty strings.", call. = FALSE)
  }
  if (length(ids) != n)
    stop("X and Y must have the same number of subjects.", call. = FALSE)
  if (is.complex(Y))
    stop("Y must use real binary or one-hot class coding.", call. = FALSE)
  encoded <- .nested_response(Y)
  if (nrow(encoded$Y) != n)
    stop("X and Y must have the same number of subjects.", call. = FALSE)
  # Use the nested pipeline's feature-fast unfolding and preserve native scores.
  Xmat <- matrix(X, n)
  colnames(Xmat) <- as.vector(outer(features, times, paste, sep = "_"))
  fit <- .tensorpls_fit(Xmat, encoded$Y, ncomp)
  if (!all(fit$convergence$converged))
    warning("TensorPLS fit did not converge for every component; inspect convergence.", call. = FALSE)
  sst <- colSums(sweep(encoded$Y, 2L, fit$preprocessing$y_mean, "-")^2)
  r2 <- stats::setNames(1 - colSums(fit$resid^2) / sst, colnames(encoded$Y))
  structure(list(ncomp = fit$ncomp, feature_names = features, time_names = times, subject_ids = ids,
    training_class = stats::setNames(factor(encoded$class, levels = encoded$levels),
      ids), x_scores = fit$x.scores, x_weights = fit$raw.wgs, x_loadings = fit$x.loads,
    coefficients = fit$reg.coefs, fitted_y = fit$y.pred, fitted_R2Y = mean(r2),
      fitted_R2Y_per_response = r2,
    response_encoding = list(levels = encoded$levels, binary = encoded$binary,
      response_names = colnames(encoded$Y)),
    preprocessing = fit$preprocessing, convergence = fit$convergence, core_fit = fit),
    class = "tensorpls_model")
}
