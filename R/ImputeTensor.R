#' Complete a tensor for full-data descriptive analysis
#'
#' @param X Numeric subjects x features x time tensor, possibly incomplete.
#' @param method Completion method: feature-time means or Tucker imputation.
#' @param ranks Three Tucker ranks in subject, feature, time order. NULL uses
#'   the fitter's reference ranks, c(4L, 4L, 3L). Ignored for mean completion.
#' @param tucker_control Named controls passed to [fit_tucker_imputer()],
#'   excluding X_train and ranks. Defaults are those of the fitter.
#' @return A numeric three-dimensional array with the dimensions and dimnames
#'   of X. Observed values are preserved.
#' @details Mean completion learns one mean per feature-time predictor from
#'   the supplied subjects. An entirely missing predictor is filled with zero.
#'   Tucker completion returns the fitted imputer's completed training tensor;
#'   the reference ranks are not claimed to be optimal for every dataset.
#'
#'   This function uses all supplied subjects and is intended for descriptive
#'   full-data completion. For predictive validation, pass incomplete X to
#'   [nested_cv_tensorpls()] and choose its fold-local imputation method.
#'   Completing X globally before cross-validation can leak information.
#' @seealso [tune_tucker_ranks()], [fit_tucker_imputer()]
#' @examples
#' X <- array(seq_len(120), c(6, 5, 4))
#' X[3] <- NA_real_
#' X_complete <- ImputeTensor(X, method = "mean")
#' @export
ImputeTensor <- function(X, method = c("mean", "tucker"), ranks = NULL, tucker_control = list()) {
  method <- match.arg(method)
  .tk_validate(X)
  if (is.complex(X))
    stop("X must be real numeric.", call. = FALSE)
  if (method == "mean") {
    if (length(tucker_control))
      stop("tucker_control applies only to method = 'tucker'.", call. = FALSE)
    predictors <- matrix(X, dim(X)[1L])
    imputer <- .nested_imputer(predictors, "training_mean")
    completed <- .nested_impute(imputer, predictors)
    return(array(completed, dim(X), dimnames(X)))
  }
  if (is.null(ranks)) {
    ranks <- c(4L, 4L, 3L)
    if (any(ranks > dim(X)))
      stop("Default Tucker ranks c(4, 4, 3) exceed tensor dimensions; supply valid ranks or run tune_tucker_ranks().", call. = FALSE)
  }
  allowed <- setdiff(names(formals(fit_tucker_imputer)), c("X_train", "ranks"))
  if (!is.list(tucker_control) || (length(tucker_control) && (is.null(names(tucker_control)) ||
    anyNA(names(tucker_control)) || anyDuplicated(names(tucker_control)) || !all(names(tucker_control) %in%
    allowed))))
    stop("tucker_control must contain named fitter controls, excluding X_train and ranks.", call. = FALSE)
  fit <- do.call(fit_tucker_imputer, c(list(X_train = X, ranks = ranks), tucker_control))
  fit$X_train_imputed
}
