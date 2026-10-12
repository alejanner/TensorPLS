\name{ImputeTensor}
\alias{ImputeTensor}
\title{Complete a tensor for full-data descriptive analysis}
\usage{
ImputeTensor(X, method = c("mean", "tucker"), ranks = NULL,
  tucker_control = list())
}
\arguments{
\item{X}{Numeric subjects by features by time tensor, possibly incomplete.}
\item{method}{Completion method: feature-time means or Tucker imputation.}
\item{ranks}{Three Tucker ranks in subject, feature, time order. NULL uses the reference ranks c(4L, 4L, 3L). Ignored for mean completion.}
\item{tucker_control}{Named controls passed unchanged to \code{fit_tucker_imputer()}, excluding X_train and ranks. Defaults are those of the fitter. Used only for Tucker completion.}
}
\description{Complete missing entries using the existing mean or Tucker imputation implementation.}
\details{Mean completion estimates one mean per feature-time predictor from the supplied subjects. Entirely missing predictors are filled with zero. Tucker completion returns the fitter's completed training tensor. The reference ranks c(4, 4, 3) are not universally optimal; invalid default ranks require explicit valid ranks or rank tuning. Warnings and errors from the fitter propagate unchanged.

This function uses all supplied subjects for descriptive full-data completion. For predictive validation, pass incomplete X to \code{nested_cv_tensorpls()} and choose its fold-local imputation method. Globally completing X before cross-validation can leak information.

Use \code{tune_tucker_ranks(X)} to recommend ranks by masked reconstruction error. Its \code{selected_ranks} can be passed to this function, or the user can deliberately choose another valid tuple after inspecting \code{results}.}
\value{A numeric three-dimensional array with the dimensions and dimnames of X. Observed values are preserved.}
\seealso{\code{\link{fit_tucker_imputer}}, \code{\link{tune_tucker_ranks}}, \code{\link{nested_cv_tensorpls}}}
\examples{
X <- array(seq_len(120), c(6, 5, 4))
X[3] <- NA_real_
X_complete <- ImputeTensor(X, method = "mean")
\dontrun{
X_complete <- ImputeTensor(X, method = "tucker")
rank_tuning <- tune_tucker_ranks(X)
rank_tuning$selected_ranks
X_complete <- ImputeTensor(X, method = "tucker",
                          ranks = rank_tuning$selected_ranks)
rank_tuning$results
X_complete <- ImputeTensor(X, method = "tucker", ranks = c(3L, 3L, 2L))
}
}
