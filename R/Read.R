#' Read a path or return the object unchanged
#'
#' Utility that detects whether `x` is a single-length character string
#' pointing to an existing file. If so, it loads the file using a reader
#' appropriate to the extension (`csv`, `tsv`, `rds`, `feather`, `parquet`);
#' otherwise returns `x` untouched.
#'
#' @param x Either an object already in memory or a string path to a file.
#' @return A data.frame, matrix, list, or the original object, depending on
#'   the input and file type.
#' @details CSV/TSV are read with `read.csv()` / `read.delim()`.  
#'   RDS via `readRDS()`. Feather / Parquet require the **arrow** package.
#' @importFrom utils read.csv read.delim
#' @importFrom tools file_ext
#' @importFrom arrow read_feather read_parquet
#' @keywords internal
#' @rdname helpers-internal
.read_or_pass <- function(x) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    ext <- tolower(file_ext(x))
    switch(ext, csv = read.csv(x, header = TRUE, check.names = FALSE), tsv = read.delim(x,
      header = TRUE, check.names = FALSE), rds = readRDS(x), feather = arrow::read_feather(x),
      parquet = arrow::read_parquet(x), stop("Unsupported file extension: ", ext))
  } else {
    x
  }
}
