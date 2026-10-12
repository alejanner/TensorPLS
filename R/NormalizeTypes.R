#' Normalize types and textual NAs according to a concise spec
#'
#' Internal helper to consistently normalize identifiers, time, and feature
#' columns prior to tensor construction.
#'
#' Behavior is controlled by `coercion_mode` and, when set to `"custom"`,
#' by `id_as`, `time_as`, and `features_numeric`.
#'
#' @param df A data.frame.
#' @param id_col Character; subject identifier column name.
#' @param time_col Character; time-point column name.
#' @param feature_cols Character vector of feature columns. If `NULL`, all
#'   columns except `id_col` and `time_col` are treated as features.
#' @param coercion_mode One of `"force_numeric"`, `"force_character"`, `"custom"`.
#'   - `"force_numeric"`: ID and TIME to numeric; features to numeric.
#'   - `"force_character"`: ID and TIME to character; features not forced.
#'   - `"custom"`: use `id_as`, `time_as`, `features_numeric`.
#' @param id_as (custom only) One of `"numeric"`, `"character"`.
#' @param time_as (custom only) One of `"numeric"`, `"character"`, `"date"`.
#' @param features_numeric Logical; if `TRUE`, coerce features to numeric.
#' @param na_strings Character values to treat as missing (default includes
#'   `""`, `"NA"`, `"N/A"`, `"n/a"`, `"."`, `"null"`).
#'
#' @return The normalized data.frame.
#' @keywords internal
#' @noRd
normalize_types_by_spec <- function(df, id_col, time_col, feature_cols,
  coercion_mode = c("force_numeric",
  "force_character", "custom"), id_as = c("numeric", "character"),
    time_as = c("numeric", "character", "date"),
  features_numeric = TRUE, na_strings = c("", "NA", "N/A", "n/a", ".",
    "null")) {
  stopifnot(is.data.frame(df))
  coercion_mode <- match.arg(coercion_mode)
  id_as <- match.arg(id_as)
  time_as <- match.arg(time_as)
  # Macro → specific settings
  if (coercion_mode != "custom") {
    if (coercion_mode == "force_numeric") {
      id_as <- "numeric"
      time_as <- "numeric"
      features_numeric <- TRUE
    } else
      if (coercion_mode == "force_character") {
      id_as <- "character"
      time_as <- "character"
      features_numeric <- FALSE
    }
  }
  .num <- function(x)
    suppressWarnings(as.numeric(gsub(",", ".", as.character(x),
    fixed = FALSE)))
  # Clean textual NAs and trim whitespace
  for (nm in names(df)) {
    x <- df[[nm]]
    if (is.character(x) || is.factor(x)) {
      x <- trimws(as.character(x))
      x[x %in% na_strings] <- NA
      df[[nm]] <- x
    }
  }
  # ID
  if (!is.null(id_col) && id_col %in% names(df)) {
    if (id_as == "numeric") {
      df[[id_col]] <- .num(df[[id_col]])
    } else {
      # character
      df[[id_col]] <- as.character(df[[id_col]])
    }
  }
  # TIME
  if (!is.null(time_col) && time_col %in% names(df)) {
    tc <- df[[time_col]]
    if (inherits(tc, c("Date", "POSIXct", "POSIXt"))) {
      if (time_as == "numeric")
        df[[time_col]] <- as.numeric(tc) else
          if (time_as == "character")
        df[[time_col]] <- as.character(tc)
      # "date": leave as is
    } else {
      if (time_as == "numeric") {
        df[[time_col]] <- .num(tc)
      } else
        if (time_as == "character") {
        df[[time_col]] <- as.character(tc)
      } else
        if (time_as == "date") {
        tc_try <- suppressWarnings(as.Date(tc))
        tc_try2 <- suppressWarnings(as.POSIXct(tc))
        if (!all(is.na(tc_try)))
          df[[time_col]] <- tc_try else
            if (!all(is.na(tc_try2)))
          df[[time_col]] <- tc_try2 else
            df[[time_col]] <- as.character(tc)
      }
    }
  }
  # FEATURES → numeric (optional)
  if (isTRUE(features_numeric)) {
    feats <- if (is.null(feature_cols))
      setdiff(names(df), c(id_col, time_col)) else
        feature_cols
    for (nm in feats)
      if (nm %in% names(df))
      df[[nm]] <- .num(df[[nm]])
  }
  df
}
