\name{plot_temporal_vip_heatmap}
\alias{plot_temporal_vip_heatmap}
\title{Plot a temporal VIP heatmap}
\usage{plot_temporal_vip_heatmap(vips, features = NULL, top_n = NULL,
  padj = NULL, alpha = 0.05)}
\arguments{
\item{vips}{A \code{tensorpls_vips} object from \code{\link{tensorpls_vips}}.}
\item{features}{Unique feature names, or NULL for all features.}
\item{top_n}{Optional positive integer limiting displayed rows to at most 50.}
\item{padj}{Optional finite numeric matrix in [0,1], with the exact dimensions,
feature/time names and order of the full \code{M2_mean} matrix.}
\item{alpha}{Marker threshold; cells with padj strictly below alpha are marked.}
}
\description{Display feature-by-time mean cumulative VIP values.}
\details{Fill represents \code{M2_mean[f,t]}, the mean VIP across the cumulative
component path. Rows rank by descending mean across time within the supplied
feature set. Ties retain input order; top_n keeps the first ranked rows.
Time order is unchanged. At most 50 rows may be displayed; use top_n or a
smaller feature set for larger inputs. No clustering or statistical calculation
occurs. Optional markers denote externally supplied adjusted p-values, separate
from descriptive VIP magnitude. No VIP2D permutation statistics are used.}
\value{A ggplot with feature, time and value columns in its data, and padj only
when supplied. No caption is produced.}
