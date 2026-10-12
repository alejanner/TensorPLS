#' Prepare omics data for downstream analysis (3D tensor)
#'
#' Reads tabular data (long or wide), applies preprocessing steps (optional
#' transpose, type/NA normalization, cohort filtering, duplicate handling),
#' and returns a 3-D numeric tensor \emph{(Subject × Feature × Time)}.
#'
#' The function always builds a 3-D tensor; if `time_col` is missing or results
#' in zero valid time points after filtering/normalization, an explicit error is
#' raised.
#'
#' @section Workflow (brief):
#' \enumerate{
#'   \item Input loading via \code{.read_or_pass()}.
#'   \item Optional transpose via \code{.transpose_if_needed()}.
#'   \item Type/NA normalization via \code{normalize_types_by_spec()} controlled by
#'         \code{coercion_mode} (and, for \code{"custom"}, \code{id_as}/\code{time_as}/\code{features_numeric}).
#'   \item Cohort filter (\code{cohort}, \code{cohort_filter}); ID type is aligned across tables.
#'   \item Minimum time-points filter.
#'   \item Guardrail: explicit error if the time axis would be empty.
#'   \item Tensor creation via \code{.build_tensor()}.
#' }
#'
#' @param data Data.frame, matrix/array, or a readable file path.
#' @param id_col Character; subject identifier column name.
#' @param time_col Character; time-point column name (required).
#' @param transpose One of \code{"auto"}, \code{"never"}, \code{"always"}; passed to
#'   \code{.transpose_if_needed()}.
#' @param header_in_row Logical; used only when transposing.
#' @param cohort Optional data.frame or file path with metadata.
#' @param cohort_id_col Column in \code{cohort} matching \code{id_col}.
#' @param cohort_filter Character scalar; an R expression evaluated within
#'   \code{cohort} to subset rows.
#' @param min_timepoints Integer; keep only subjects observed at least this
#'   many times.
#' @param exclude_cols Columns to drop entirely when building features.
#' @param feature_cols Character vector; explicit feature columns to keep.
#'   Overrides \code{exclude_cols}.
#' @param subjects Optional character vector defining order (and subset) of
#'   subjects in the tensor.
#' @param time_points Optional vector defining order (and subset) of time points
#'   in the tensor.
#' @param coercion_mode One of \code{"force_numeric"}, \code{"force_character"}, \code{"custom"}.
#'   \itemize{
#'     \item \code{"force_numeric"}: ID/TIME coerced to numeric; features numeric; forces \code{deduplicate="first"}.
#'     \item \code{"force_character"}: ID/TIME as character; features not forced.
#'     \item \code{"custom"}: use \code{id_as}, \code{time_as}, \code{features_numeric}.
#'   }
#' @param id_as (custom only) One of \code{"numeric"}, \code{"character"}.
#' @param time_as (custom only) One of \code{"numeric"}, \code{"character"}, \code{"date"}.
#' @param features_numeric Logical; if \code{TRUE} (default), coerce features to numeric.
#' @param deduplicate Strategy when duplicate subject/time pairs are found:
#'   \code{"first"}, \code{"last"}, or \code{"mean"}.
#'
#' @return 3-D numeric array (Subject × Feature × Time).
#'
#' @examples
#' data <- data.frame(subject = rep(1:4, each = 3),
#'                    time = rep(1:3, 4),
#'                    feature1 = 1:12, feature2 = 13:24)
#' arr <- prepare_omics(data, id_col = "subject", time_col = "time",
#'                      transpose = "never")
#' dim(arr)
#'
#' @export
prepare_omics <- function(data, id_col, time_col, transpose = "auto", header_in_row = TRUE,
  cohort = NULL, cohort_id_col = id_col, cohort_filter = NULL, min_timepoints = NULL,
  exclude_cols = NULL, feature_cols = NULL, subjects = NULL, time_points = NULL,
    coercion_mode = c("force_numeric",
    "force_character", "custom"), id_as = c("numeric", "character"),
      time_as = c("numeric", "character",
    "date"), features_numeric = TRUE, deduplicate = c("first", "last", "mean")) {
  coercion_mode <- match.arg(coercion_mode)
  id_as <- match.arg(id_as)
  time_as <- match.arg(time_as)
  deduplicate <- match.arg(deduplicate)
  if (coercion_mode == "force_numeric")
    deduplicate <- "first"
  # I/O + transpose
  obj <- .read_or_pass(data)
  if (is.matrix(obj) || (is.array(obj) && length(dim(obj)) == 3L))
    return(obj)
  if (!is.data.frame(obj))
    stop("Unsupported input type: ", paste(class(obj), collapse = ", "))
  obj <- .transpose_if_needed(obj, transpose, header_in_row)
  # Require a time column before normalization.
  if (is.null(time_col) || !(time_col %in% names(obj))) {
    stop("`time_col` is missing or not present in the data; required for a 3-D tensor.")
  }
  #  Type/NA normalization
  obj <- normalize_types_by_spec(df = obj, id_col = id_col, time_col = time_col,
    feature_cols = feature_cols,
    coercion_mode = coercion_mode, id_as = id_as, time_as = time_as,
      features_numeric = features_numeric)
  #  Cohort filter (align ID types consistently across tables)
  if (!is.null(cohort)) {
    cdf <- .read_or_pass(cohort)
    if (!is.data.frame(cdf))
      stop("Cohort must be a data.frame or readable file.")
    if (cohort_id_col %in% names(cdf))
      names(cdf)[names(cdf) == cohort_id_col] <- id_col
    if (!is.null(cohort_filter))
      cdf <- cdf[with(cdf, eval(parse(text = cohort_filter))), ]
    # Align ID type for the %in% join
    join_as_numeric <- FALSE
    if (coercion_mode == "force_numeric") {
      join_as_numeric <- TRUE
    } else
      if (coercion_mode == "custom")
      {
        join_as_numeric <- (id_as == "numeric")
      }  # force_character -> FALSE
    if (join_as_numeric) {
      obj[[id_col]] <- suppressWarnings(as.numeric(as.character(obj[[id_col]])))
      cdf[[id_col]] <- suppressWarnings(as.numeric(as.character(cdf[[id_col]])))
    } else {
      obj[[id_col]] <- as.character(obj[[id_col]])
      cdf[[id_col]] <- as.character(cdf[[id_col]])
    }
    idx <- match(obj[[id_col]], cdf[[id_col]])
    obj <- obj[!is.na(idx), , drop = FALSE]
    if (nrow(obj) == 0L || !any(!is.na(obj[[time_col]]))) {
      stop("No valid time points after filtering/normalization; check ID/TIME types, cohort, and coercion settings.")
    }
  }
  #  Minimum time-points
  if (!is.null(min_timepoints) && !is.null(time_col) && time_col %in% names(obj)) {
    distinct_t_by_id <- tapply(obj[[time_col]], obj[[id_col]], function(x)
      length(unique(x)))
    keep <- names(distinct_t_by_id)[distinct_t_by_id >= min_timepoints]
    obj <- obj[obj[[id_col]] %in% keep, ]
  }
  # Validate time points after filtering.
  if (is.null(time_col) || !(time_col %in% names(obj))) {
    stop("`time_col` is missing or not present in the data; required for a 3-D tensor.")
  }
  tp_vals <- unique(obj[[time_col]][!is.na(obj[[time_col]])])
  if (nrow(obj) == 0L || length(tp_vals) == 0L) {
    stop("No valid time points after filtering/normalization; check ID/TIME types, cohort, and coercion settings.")
  }
  # (5) Build tensor
  excl <- unique(c(id_col, time_col, if (is.null(exclude_cols))
    character() else
      exclude_cols))
  .build_tensor(obj, id_col = id_col, time_col = time_col, exclude_cols = excl,
    feature_cols = feature_cols,
    subjects = subjects, time_points = time_points, deduplicate = deduplicate)
}
