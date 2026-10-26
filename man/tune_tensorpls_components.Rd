\name{tune_tensorpls_components}
\alias{tune_tensorpls_components}
\title{Choose components by full-development-data cross-validation}
\usage{
tune_tensorpls_components(X, Y, ncomp = 1:10,
  feature_strategy = "Original", vip_thr = 99, cv_folds = 5L,
  seed = 123L, folds = NULL,
  imputation = c("none", "training_mean", "tucker"),
  tucker_control = list(ranks = c(4L, 4L, 3L)),
  input_provenance = "raw", selection_metric = c("Q2Y", "BER"),
  verbose = TRUE)
}
\arguments{
\item{X}{Numeric subjects by features by time tensor with unique axis labels.}
\item{Y}{Classification response accepted by \code{nested_cv_tensorpls()}.}
\item{ncomp}{Candidate integers from 1 to 10; duplicates removed.}
\item{feature_strategy}{One fixed strategy accepted by \code{vip_feature_sets()},
including M1 (all-time mean), M2 (cumulative-path mean), and their
ordered intersections. VIP3Dmodel2 uses cumulative-path sum selection.}
\item{vip_thr}{One percentile in [0,100]; ignored for Original.}
\item{cv_folds}{Number of stratified folds, at least two.}
\item{seed}{Nonnegative integer seed; caller RNG state is preserved.}
\item{folds}{Optional numeric fold vector named in exact X subject order. Labels must be consecutive integers 1 through cv_folds. Every training and validation split must contain all classes.}
\item{imputation}{Training-only imputation: none, training_mean, or tucker.}
\item{tucker_control}{Named fixed Tucker controls, as in \code{nested_cv_tensorpls()}.}
\item{input_provenance}{raw, unknown, or preimputed. Non-raw input warns about earlier preprocessing leakage.}
\item{selection_metric}{Q2Y (maximize pooled Q2Y) or BER (minimize pooled BER).}
\item{verbose}{Report development progress.}
}
\description{Choose H using one stratified cross-validation on the supplied development dataset and the existing nested development engine.}
\details{Scores are tuning diagnostics, not unbiased predictive performance; use \code{nested_cv_tensorpls()} for predictive assessment. The selected component count is intended for subsequent final full-data fitting and interpretation. Imputation and feature selection are fitted within each training fold. Q2Y uses pooled held-out errors and corresponding training response means. BER selection scores are negative pooled BER. Ties prefer fewer components, then candidate order. A candidate must complete every validation fold without diagnostics. The shared engine then refits the chosen pipeline on all supplied data. Original stores its irrelevant threshold as NA.}
\value{A \code{tensorpls_component_tuning} list with selected_ncomp, tuning (candidate scores/configurations), fold_results (shared-engine diagnostics), folds, settings, diagnostics, and refit (internal full-data pipeline). selected_ncomp is NA and refit is NULL if development or final refitting fails; diagnostics$error explains the failure. Scores remain available.}
