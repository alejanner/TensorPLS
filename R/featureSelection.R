#The function is agnostic to the VIP object name: pass `VIP3Dmodel1`,
#' `VIP3Dmodel2`, `VIP2D`, or a mixOmics `vip()` matrix.
#' @param vip         Matrix/data‑frame of VIP scores (variables × columns).
#' @param thr         Percentile threshold (0–100). Default 95.
#' @param cols        Columns to use (numeric indices or names). `NULL` keeps
#'                    all columns.
#' @param aggregator  How to summarise across selected columns:
#'                    * "any"  – keep if passes **one** column (logical OR)
#'                    * "mean" – percentile filter on the **row‑wise mean**
#' @param strip_time  Logical. Remove the suffix matched by `pattern`.
#' @param pattern     Regex used when `strip_time = TRUE` (default "\|t\\d+$").
#' @param na.rm       Logical passed to `quantile()` and row operations.
#'
#' @return Character vector of selected variables (deduplicated)
#' @export
feature_selection <- function(vip, thr = 95, cols = NULL, aggregator = c("any",
  "mean", "all"), strip_time = TRUE, pattern = "\\|t\\d+$", na.rm = TRUE) {
  aggregator <- match.arg(aggregator)
  vip_names <- rownames(vip)
  vip_df <- as.data.frame(vip)
  if (is.null(vip_names))
    vip_names <- rownames(vip_df)
  if (!is.null(cols)) {
    if (is.numeric(cols)) {
      bad <- cols[cols < 1 | cols > ncol(vip_df)]
      if (length(bad))
        stop(sprintf("`cols` index out of range 1..%d (provided: %s)", ncol(vip_df),
          paste(bad, collapse = ", ")), call. = FALSE)
    } else {
      missing <- setdiff(cols, colnames(vip_df))
      if (length(missing))
        stop(sprintf("`cols` has unknown column name(s): %s. Available: %s",
          paste(missing, collapse = ", "), paste(colnames(vip_df),
          collapse = ", ")), call. = FALSE)
    }
    vip_df <- vip_df[, cols, drop = FALSE]
  }
  if (ncol(vip_df) == 0)
    stop("No columns selected in vip matrix.", call. = FALSE)
  # Compute a separate percentile threshold for each column.
  bool_matrix <- do.call(cbind, lapply(vip_df, function(col) {
    thr_val <- stats::quantile(col, thr / 100, na.rm = na.rm, names = FALSE)
    col >= thr_val
  }))
  keep <- switch(aggregator, any = apply(bool_matrix, 1, any), all = apply(bool_matrix,
    1, all), mean = {
    row_mean <- rowMeans(vip_df, na.rm = na.rm)
    row_mean >= stats::quantile(row_mean, thr / 100, names = FALSE, na.rm = na.rm)
  })
  vars <- vip_names[keep]
  if (strip_time)
    vars <- sub(pattern, "", vars)
  vars <- unique(vars)
  vars
}
