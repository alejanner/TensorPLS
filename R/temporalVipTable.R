#' Export feature-level temporal VIP information
#'
#' @param vips A `tensorpls_vips` object from [tensorpls_vips()].
#' @param features Unique feature names, or NULL for all features.
#' @param nested Optional `nested_tensorpls` result from [nested_cv_tensorpls()].
#' @param order_by Descending mean M2, peak M2, or outer-refit selection frequency.
#' @details Uses only `M2_mean`. mean_M2 averages a feature's temporal VIP across
#'   time; peak_M2 is its maximum and peak_time is the first maximum in model
#'   time order. Ranking ties retain input feature order. No row limit applies.
#'
#'   Nested counts use `selected_features_per_outer_fold`, populated from the
#'   final outer-training refits. The denominator is the number of nonempty
#'   refit selections, matching the nested object's stability calculation.
#'   A completed refit remains eligible if subsequent test prediction fails;
#'   failed development with no refit is excluded. No inner candidate sets are
#'   used. `stability$frequency$feature` validates the nested feature universe.
#'   Frequencies describe selection stability, not p-values or biological
#'   validation. The supplied nested result must concern compatible features.
#' @return A standard data.frame with feature, mean_M2, peak_M2, peak_time,
#'   then M2_<time> columns in model time order. With nested, appends
#'   cv_selection_count and cv_selection_frequency. No model-level metrics
#'   or selection-rule evidence are included. Suitable for write.csv().
#' @export
temporal_vip_table <- function(vips, features = NULL, nested = NULL, order_by = c("mean_M2",
  "peak_M2", "cv_selection_frequency")) {
  order_by <- match.arg(order_by)
  if (order_by == "cv_selection_frequency" && is.null(nested))
    stop("order_by = 'cv_selection_frequency' requires nested.", call. = FALSE)
  checked <- temporal_vip_summary(vips, features)
  features <- attr(checked, "feature_names")
  m <- vips[["M2_mean", exact = TRUE]][features, , drop = FALSE]
  peak <- max.col(m, ties.method = "first")
  result <- data.frame(feature = features, mean_M2 = rowMeans(m),
    peak_M2 = m[cbind(seq_len(nrow(m)),
    peak)], peak_time = vips$time_names[peak], stringsAsFactors = FALSE, row.names = NULL,
    check.names = FALSE)
  times <- as.data.frame(m, optional = TRUE)
  names(times) <- paste0("M2_", vips$time_names)
  rownames(times) <- NULL
  result <- cbind(result, times)
  if (!is.null(nested)) {
    if (!inherits(nested, "nested_tensorpls") || !is.list(nested))
      stop("nested must be a nested_cv_tensorpls() result.", call. = FALSE)
    sets <- nested$selected_features_per_outer_fold
    ids <- names(sets)
    if (!is.list(sets) || !length(sets) || is.null(ids) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids))
      stop("nested must contain named selected_features_per_outer_fold.", call. = FALSE)
    k <- nested$settings$outer_folds
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1 ||
      k != floor(k) || length(sets) != k || !identical(ids, as.character(seq_len(k))))
      stop("nested must retain every declared outer fold in order.", call. = FALSE)
    universe <- nested$stability$frequency$feature
    if (!is.character(universe) || !length(universe) || anyNA(universe) || any(!nzchar(universe)) ||
      anyDuplicated(universe))
      stop("nested must contain a valid stability$frequency$feature universe.", call. = FALSE)
    if (!all(features %in% universe))
      stop("Requested features are absent from the nested feature universe.", call. = FALSE)
    valid <- vapply(sets, function(s)
      is.null(s) || (is.character(s) && !anyNA(s) &&
      all(nzchar(s)) && !anyDuplicated(s) && all(s %in% universe)), logical(1))
    if (!all(valid))
      stop("Invalid outer-refit selected feature sets in nested.", call. = FALSE)
    eligible <- lengths(sets) > 0L
    denominator <- sum(eligible)
    if (!denominator)
      stop("nested contains no eligible outer-training refit selections.", call. = FALSE)
    counts <- table(unlist(sets[eligible], use.names = FALSE))
    result$cv_selection_count <- as.integer(counts[match(features, names(counts))])
    result$cv_selection_count[is.na(result$cv_selection_count)] <- 0L
    result$cv_selection_frequency <- result$cv_selection_count / denominator
  }
  result <- result[order(-result[[order_by]], seq_len(nrow(result))), , drop = FALSE]
  rownames(result) <- NULL
  result
}
