#' Wrapper selecting 2-D or 3-D tensor representation. For future extensions on 
#' 2D matrix. 
#'
#' Dispatches to `.create_array3d()` when a valid `time_col` is provided and
#' present in `df`; 
#'
#' @param time_col Optional; if `NULL` or not found in `df`, a 2-D matrix is
#'   produced.
#' @param ... Additional arguments passed on to the chosen helper.
#'
#' @return Either a 3-D numeric array or a 2-D numeric matrix.
#'
#' @keywords internal
#' @rdname helpers-internal
#' 
#'
#'
#'
#'
#' @keywords internal
.build_tensor <- function(df, id_col, time_col, ..., deduplicate = c("first", "last",
  "mean")) {
  if (is.null(time_col) || !(time_col %in% names(df))) {
    stop("`time_col` is missing or not present in the data; required for a 3-D tensor.")
  }
  .create_array3d(df, id_col = id_col, time_col = time_col, ...,
    deduplicate = match.arg(deduplicate))
}
