\name{plot.tensorpls_component_tuning}
\alias{plot.tensorpls_component_tuning}
\title{Plot stored development component-tuning scores}
\usage{
\method{plot}{tensorpls_component_tuning}(x, show_folds = FALSE,
  mark_selected = TRUE, main = NULL, ...)
}
\arguments{
\item{x}{A tensorpls_component_tuning object.}
\item{show_folds}{Optionally overlay stored fold-specific diagnostic trajectories; omitted by default.}
\item{mark_selected}{Mark the stored selected component, when available.}
\item{main}{Optional title.}
\item{...}{Additional base plot arguments, such as ylim or xlab.}
}
\description{Display the stored component tuning curve without fitting or selecting a model.}
\details{The pooled curve uses tuning$selection_score directly. Q2Y is unchanged; negative BER scores are negated for display, so lower BER is better. Fold curves use stored Q2Y or BER values and are not averaged. These are development tuning diagnostics, not unbiased predictive performance. Nonfinite scores leave gaps. With no finite scores an empty frame is drawn. A missing selected component does not prevent plotting finite scores; no selection is invented.}
\value{Invisibly, a data.frame with ncomp, metric, pooled_score on the displayed scale, and selected, ordered by H. Selection reflects the stored component even when mark_selected is FALSE. The input is unchanged.}
