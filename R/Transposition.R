#' Conditionally transpose a data.frame
#'
#' Decides—according to `mode` and the first column contents—whether to
#' transpose the input `df`. If `header_in_row = TRUE`, the first row after
#' transposition becomes the column names.
#'
#' @param df A data.frame to be examined and possibly transposed.
#' @param mode Character; `"auto"`, `"never"`, or `"always"`.
#' @param header_in_row Logical; if `TRUE`, first transposed row ⇒ column
#'   names.
#'
#' @return A data.frame, transposed if criteria met, else the original `df`.
#'
#' @keywords internal
#' @rdname helpers-internal
#' 
#' 
.transpose_if_needed <- function(df, mode = c("auto", "never", "always"), header_in_row = TRUE) {
  mode <- match.arg(mode)
  if (!is.data.frame(df))
    return(df)
  do_it <- switch(mode, never = FALSE, always = TRUE, auto = {
    if (ncol(df) == 0 || nrow(df) == 0) {
      FALSE
    } else {
      # Inspect the first column for feature identifiers.
      fc <- df[[1]]
      # Transpose when the first column contains unique nonnumeric identifiers.
      !is.numeric(fc) && length(unique(fc)) == nrow(df)
    }
  })
  if (!do_it)
    return(df)
  tmp <- t(df)
  if (header_in_row) {
    colnames(tmp) <- tmp[1, ]
    tmp <- tmp[-1, , drop = FALSE]
  }
  as.data.frame(tmp, stringsAsFactors = FALSE)
}
