#' Build a 3-D numeric array (subjects × features × time)
#'
#' Robust 3-D tensor construction with tolerant numeric conversion for features
#' and flexible time handling (Date/POSIXct, numeric-like, or character labels).
#'
#' @importFrom stats aggregate
#' @keywords internal
#' @noRd
.create_array3d <- function(df, id_col, time_col, feature_cols = NULL, exclude_cols = c(id_col,
  time_col), subjects = NULL, time_points = NULL, deduplicate = c("first", "last",
  "mean")) {
  deduplicate <- match.arg(deduplicate)
  #  helpers
  .num <- function(x)
    suppressWarnings(as.numeric(gsub(",", ".", as.character(x),
    fixed = FALSE)))
  .is_date_time <- function(x)
    inherits(x, c("Date", "POSIXct", "POSIXt"))
  # Validate input columns.
  if (!is.data.frame(df)) {
    stop("'.create_array3d' expects a data.frame.")
  }
  if (!all(c(id_col, time_col) %in% names(df))) {
    stop("Missing key columns: expected '", id_col, "' and '", time_col, "'.")
  }
  #  ID as character (stable key)
  df[[id_col]] <- as.character(df[[id_col]])
  #  feature set
  if (is.null(feature_cols)) {
    feature_cols <- setdiff(names(df), exclude_cols)
  }
  if (length(feature_cols) == 0L) {
    stop("No features available after 'exclude_cols' / 'feature_cols'.")
  }
  if (!all(feature_cols %in% names(df))) {
    miss <- setdiff(feature_cols, names(df))
    stop("Unknown feature columns: ", paste(miss, collapse = ", "))
  }
  #  time handling
  tc <- df[[time_col]]
  time_is_date <- .is_date_time(tc)
  if (time_is_date) {
    # no NaN in time points
    if (is.null(time_points)) {
      time_points <- sort(unique(tc[!is.na(tc)]))
    } else {
      if (inherits(tc, "Date")) {
        time_points <- as.Date(time_points)
      } else {
        tz <- attr(tc, "tzone")
        if (is.null(tz))
          tz <- ""
        time_points <- as.POSIXct(time_points, tz = tz)
      }
      time_points <- time_points[!is.na(time_points)]
    }
  } else {
    tc_as_char <- as.character(tc)
    is_num_like <- grepl("^\\s*[-+]?\\d+(?:[\\.,]\\d+)?\\s*$", tc_as_char)
    if (all(is_num_like | is.na(tc_as_char), na.rm = TRUE)) {
      # numeric-like
      df[[time_col]] <- .num(tc_as_char)
      if (is.null(time_points)) {
        time_points <- sort(unique(df[[time_col]][!is.na(df[[time_col]])]))
      } else {
        time_points <- .num(time_points)
        time_points <- time_points[!is.na(time_points)]
      }
      time_is_date <- FALSE
    } else {
      # character: no NaN
      df[[time_col]] <- tc_as_char
      if (is.null(time_points)) {
        tp0 <- unique(df[[time_col]])
        time_points <- tp0[!is.na(tp0) & tp0 != "NA" & tp0 != ""]
      } else {
        time_points <- as.character(time_points)
        time_points <- time_points[!is.na(time_points) & time_points != "NA" &
          time_points != ""]
      }
    }
  }
  #  subjects order/set
  if (is.null(subjects)) {
    subjects <- unique(df[[id_col]])
  }
  subjects <- as.character(subjects)
  #  destination array
  arr <- array(NA_real_, dim = c(length(subjects), length(feature_cols), length(time_points)),
    dimnames = list(Subject = subjects, Feature = feature_cols, Time = as.character(time_points)))
  #  fill slices
  for (k in seq_along(time_points)) {
    tp <- time_points[k]
    # time match by type
    if (time_is_date) {
      slice <- df[df[[time_col]] == tp, c(id_col, feature_cols), drop = FALSE]
    } else
      if (is.numeric(df[[time_col]])) {
      slice <- df[df[[time_col]] == as.numeric(tp), c(id_col, feature_cols), drop = FALSE]
    } else {
      slice <- df[df[[time_col]] == as.character(tp), c(id_col, feature_cols), drop = FALSE]
    }
    if (nrow(slice) > 0L) {
      slice[[id_col]] <- as.character(slice[[id_col]])
    }
    # deduplicate by subject
    if (nrow(slice) > 0L && any(duplicated(slice[[id_col]]))) {
      if (deduplicate == "first") {
        slice <- slice[!duplicated(slice[[id_col]]), , drop = FALSE]
      } else
        if (deduplicate == "last") {
        slice <- slice[!duplicated(slice[[id_col]], fromLast = TRUE), , drop = FALSE]
      } else
        if (deduplicate == "mean") {
        # force numeric before to aggregate
        for (cc in feature_cols)
          slice[[cc]] <- .num(slice[[cc]])
        agg <- aggregate(slice[feature_cols], by = list(.id = slice[[id_col]]),
          FUN = function(z) {
          m <- mean(z, na.rm = TRUE)
          if (is.nan(m))
            NA_real_ else
              m
          })
        names(agg)[1] <- id_col
        slice <- agg
      }
    }
    # anchor to preserve subject order
    anchor <- data.frame(tmp = subjects, stringsAsFactors = FALSE)
    names(anchor)[1] <- id_col
    anchor[[id_col]] <- as.character(anchor[[id_col]])
    merged <- merge(anchor, slice, by = id_col, all.x = TRUE, sort = FALSE)
    # reorder strictly to anchor
    ord <- match(anchor[[id_col]], merged[[id_col]])
    if (anyNA(ord)) {
      # Handle missing rows properly
      missing_rows <- is.na(ord)
      fill <- matrix(NA_real_, nrow = sum(missing_rows), ncol = length(feature_cols))
      colnames(fill) <- feature_cols
      fill_df <- as.data.frame(fill)
      # Create ID column dataframe
      id_df <- data.frame(subjects[missing_rows], stringsAsFactors = FALSE)
      names(id_df) <- id_col
      # Get valid merged rows
      merged_valid <- merged[!is.na(match(merged[[id_col]], anchor[[id_col]])),
        , drop = FALSE]
      # Combine valid and missing rows
      if (nrow(merged_valid) > 0) {
        valid_ord <- match(anchor[[id_col]][!missing_rows], merged_valid[[id_col]])
        merged <- rbind(merged_valid[valid_ord, , drop = FALSE], cbind(id_df,
          fill_df))
      } else {
        merged <- cbind(id_df, fill_df)
      }
      # Final reorder
      ord2 <- match(anchor[[id_col]], merged[[id_col]])
      merged <- merged[ord2, , drop = FALSE]
    } else {
      merged <- merged[ord, , drop = FALSE]
    }
    # force features to numeric (tolerant)
    if (nrow(merged) > 0L) {
      for (cc in feature_cols)
        merged[[cc]] <- .num(merged[[cc]])
    }
    mat <- as.matrix(merged[, feature_cols, drop = FALSE])
    storage.mode(mat) <- "double"
    arr[, , k] <- mat
  }
  arr
}
