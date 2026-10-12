\name{temporal_vip_table}
\alias{temporal_vip_table}
\title{Export feature-level temporal VIP information}
\usage{temporal_vip_table(vips, features = NULL, nested = NULL,
  order_by = c("mean_M2", "peak_M2", "cv_selection_frequency"))}
\arguments{
\item{vips}{A \code{tensorpls_vips} object from \code{\link{tensorpls_vips}}.}
\item{features}{Unique feature names, or NULL for all features.}
\item{nested}{Optional \code{nested_tensorpls} result from \code{\link{nested_cv_tensorpls}}.}
\item{order_by}{Descending mean M2, peak M2, or outer-refit selection frequency.}
}
\description{Return an ordered feature-level data frame of temporal VIP values.}
\details{Uses only \code{M2_mean}. mean_M2 averages a feature's temporal VIP
across time; peak_M2 is its maximum and peak_time is the first maximum in
model time order. Ranking ties retain input feature order. No row limit applies.

Nested counts use \code{selected_features_per_outer_fold}, populated from
final outer-training refits. The denominator is the number of nonempty refit
selections, matching the nested object's stability calculation. A completed
refit remains eligible if subsequent test prediction fails; failed development
with no refit is excluded. No inner candidate sets are used.
\code{stability$frequency$feature} validates the nested feature universe.
Frequencies describe selection stability, not p-values or biological validation.
The nested result must concern compatible features.}
\value{A standard data.frame with feature, mean_M2, peak_M2, peak_time, then
M2_<time> columns in model time order. With nested, appends cv_selection_count
and cv_selection_frequency. No model-level metrics or selection-rule evidence
are included. Suitable for write.csv().}
