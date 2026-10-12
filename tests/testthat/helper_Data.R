# Synthetic fixtures used across tests
# Basic long-format table:
#  - 2 subjects: P1, P2
#  - 2 time points: 0, 1
#  - 2 features: FeatA, FeatB
make_long_df <- function() {
  expand.grid(SubjectID = c("P1", "P2"), Visit = c(0, 1), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE) |>
    transform(FeatA = c(1, 2, 3, 4), FeatB = c(10, 20,
      30, 40))
}
# Same as above but SubjectID is numeric-like
make_long_numeric_id <- function() {
  df <- make_long_df()
  df$SubjectID <- c(101, 202, 101, 202)
  df
}
# Time as character labels instead of numeric
make_time_as_labels <- function() {
  df <- make_long_df()
  df$Visit <- ifelse(df$Visit == 0, "t0", "t1")
  names(df)[names(df) == "Visit"] <- "VisitLabel"
  df
}
# Minimal cohort table:
#  - CohortKey matches the ID column in the data
#  - SetTag marks Train/Valid/etc.
make_cohort_df <- function(ids) {
  data.frame(CohortKey = ids, SetTag = rep(c("Train", "Train", "Valid", "Train"),
    length.out = length(ids)),
    stringsAsFactors = FALSE)
}
# A table that is "wide" across samples and requires transpose to become long.
# The first column is a non-numeric unique label → transpose(auto) should trigger
# given your .transpose_if_needed() heuristic.
make_transposed_like <- function() {
  # Structure (rows = features, columns = subject+time "cells"):
  #   FeatureKey  P1_v0 P1_v1 P2_v0 P2_v1
  #   FeatA         1     2     3     4
  #   FeatB        10    20    30    40
  data.frame(FeatureKey = c("FeatA", "FeatB"), P1_v0 = c(1, 10), P1_v1 = c(2,
    20), P2_v0 = c(3, 30), P2_v1 = c(4, 40), stringsAsFactors = FALSE)
}
