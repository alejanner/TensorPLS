#' Build named candidate feature sets from VIP scores
#'
#' Apply [feature_selection()] to each requested VIP definition, then construct
#' individual sets and all pairwise/triple intersections or unions.
#'
#' @param vip_object A named list containing the
#'   requested numeric VIP matrices/data frames with explicit feature row names.
#' @param thr Percentile in `[0,100]`, passed to [feature_selection()].
#' @param strategy Optional single strategy for a `tensorpls_vips` object:
#'   `Original`, `VIP2D`, `M1`, `M2`, `VIP3Dmodel2`, ordered intersections
#'   of VIP2D, M1 and M2, or `VIP2D_intersect_VIP3Dmodel2`.
#'   NULL returns all requested sets.
#'   Single-strategy mode uses all columns, `aggregator = "any"`, exact
#'   feature/time mapping, and the default sources/combinations/options.
#' @param sources Distinct names from `VIP2D`, `M1`, `M2`, `VIP3Dmodel2`,
#'   in output order. Defaults to `c("VIP2D", "M1", "M2")`.
#' @param combinations Any of `"individual"`, `"intersection"`, `"union"`.
#' @param cols,aggregator,na.rm Passed to [feature_selection()]. The same column
#'   selection is applied to each source; use column indices when names differ.
#' @param strip_time Whether to resolve feature/time identifiers to features.
#' @param pattern Optional suffix regex passed to [feature_selection()] for
#'   every source. With `NULL`, only VIP2D names are mapped using the exact
#'   feature/time names in VIP3Dmodel2 (`feature_time`). This avoids stripping
#'   meaningful suffixes from gene names. Nonstandard names require a pattern.
#'
#' @return An S3 `vip_feature_sets` list: `sets` (named character vectors),
#'   `summary` (feature_set, n_features, operation and list-column sources),
#'   `selection` (source-level counts, including collapsed identifiers), and
#'   `settings`. Empty intersections are retained. Repeated identifiers are
#'   deduplicated after percentile selection, including repeated time rows.
#'   With `strategy` supplied, returns only `selected_features`,
#'   `selected_indices` (in input feature order), `strategy`, `percentile`
#'   (thr/100), and `thresholds` (named per-column vectors for used sources).
#'   `Original` returns every feature and an empty threshold list.
#' @details M1 reads `VIP3Dmodel1`, the mean across all times for each
#'   feature/component pair. M2 reads `M2_mean`, the mean across the
#'   cumulative VIP path for each feature/time pair. `VIP3Dmodel2` reads
#'   the sum across that path. Sum and mean selection can differ at
#'   floating-point quantile boundaries.
#'
#'   Single-strategy selection uses all cumulative VIP2D and M1 columns.
#'   Each column has its own type-7 quantile at thr/100; a row passes if its
#'   value is greater than or equal to that threshold in any column.
#'   VIP2D quantiles span all feature/time rows within each component column;
#'   selected rows map to features and are deduplicated. M1 quantiles span
#'   features per cumulative column; M2 and VIP3Dmodel2 quantiles span
#'   features per time. Intersections combine these feature-level sets.
#'   No fitting or VIP recomputation occurs, and the input is unchanged.
#'   Percentiles include all ties and need not select a fixed number of genes.
#'
#'   This is exploratory selection on the supplied VIP scores, not
#'   independent validation. With `strip_time = FALSE`, composite VIP2D names
#'   remain unchanged and may not match a tensor's feature names.
#' @examples
#' vip <- matrix(c(1, 4, 3, 2), 4, 1,
#'               dimnames = list(LETTERS[1:4], "Comp1"))
#' sets <- vip_feature_sets(list(VIP3Dmodel1 = vip), thr = 50,
#'                         sources = "M1", combinations = "individual")
#' sets$sets
#' @export
vip_feature_sets <- function(vip_object, thr = 95, sources = c("VIP2D", "M1", "M2"),
  combinations = c("individual", "intersection"), cols = NULL,
  aggregator = c("any", "mean", "all"), strip_time = TRUE, pattern = NULL,
  na.rm = TRUE, strategy = NULL) {
  aggregator <- match.arg(aggregator)
  if (!is.null(strategy)) {
    if (
      !identical(sources, c("VIP2D", "M1", "M2")) ||
      !identical(combinations, c("individual", "intersection")) ||
      !is.null(cols) ||
      aggregator != "any" ||
      !identical(strip_time, TRUE) ||
      !is.null(pattern) ||
      !identical(na.rm, TRUE)
    )
      stop("Single-strategy mode requires default sources, combinations, cols, aggregator, strip_time, pattern and na.rm.", call. = FALSE)
    return(.vip_select_strategy(vip_object, strategy, thr))
  }
  if (!is.numeric(thr) || length(thr) != 1L || !is.finite(thr) || thr < 0 ||
    thr > 100)
    stop("thr must be one finite percentile between 0 and 100.", call. = FALSE)
  .vip_flag(strip_time, "strip_time")
  .vip_flag(na.rm, "na.rm")
  fields <- .vip_source_fields()
  allowed <- names(fields)
  if (!is.character(sources) || !length(sources) || anyNA(sources) || anyDuplicated(sources) ||
    !all(sources %in% allowed))
    stop("sources must be distinct VIP2D, M1, M2 or VIP3Dmodel2 names.", call. = FALSE)
  backing <- unname(fields[sources])
  if (!is.list(vip_object) || !all(backing %in% names(vip_object)))
    stop("vip_object is missing requested VIP sources: ", paste(setdiff(backing,
      names(vip_object)), collapse = ", "),
      call. = FALSE)
  if (!is.character(combinations) || !length(combinations) || anyNA(combinations) ||
    !all(combinations %in% c("individual", "intersection", "union")))
    stop("combinations must contain individual, intersection or union.", call. = FALSE)
  combinations <- unique(combinations)
  if (!is.null(pattern)) {
    if (!is.character(pattern) || length(pattern) != 1L || is.na(pattern))
      stop("pattern must be NULL or one regular expression.", call. = FALSE)
    tryCatch(sub(pattern, "", "test"), error = function(e)
      stop("Invalid pattern: ", conditionMessage(e),
      call. = FALSE))
  }
  base <- stats::setNames(vector("list", length(sources)), sources)
  audit <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    src <- sources[i]
    v <- vip_object[[fields[[src]]]]
    if (
      !(is.matrix(v) || is.data.frame(v)) ||
      !nrow(v) ||
      !ncol(v) ||
      !all(vapply(as.data.frame(v), is.numeric, logical(1)))
    )
      stop(src, " must be a nonempty numeric matrix or data frame.", call. = FALSE)
    ids <- rownames(v)
    if (
      is.null(ids) ||
      anyNA(ids) ||
      any(!nzchar(ids)) ||
      (is.data.frame(v) && .row_names_info(v, 1L) < 0L)
    )
      stop(src, " requires explicit, nonempty feature row names.", call. = FALSE)
    if (any(is.infinite(as.matrix(v))))
      stop(src, " contains infinite VIP scores.", call. = FALSE)
    native <- strip_time && is.null(pattern) && src == "VIP2D"
    if (native) {
      ref <- vip_object$VIP3Dmodel2
      features <- rownames(ref)
      times <- colnames(ref)
      if (is.null(features) || is.null(times) || anyDuplicated(features) || anyNA(features) ||
        anyNA(times) || anyDuplicated(times))
        stop("Automatic VIP2D mapping needs unique VIP3Dmodel2 feature/time names; supply pattern.", call. = FALSE)
      keys <- as.vector(outer(features, times, paste, sep = "_"))
      if (anyDuplicated(keys) || anyNA(match(ids, keys)))
        stop("VIP2D feature/time names are inconsistent or ambiguous; supply an explicit pattern.", call. = FALSE)
      v <- as.matrix(v)
      rownames(v) <- rep(features, length(times))[match(ids, keys)]
    }
    args <- list(vip = v, thr = thr, cols = cols, aggregator = aggregator,
      strip_time = strip_time && !is.null(pattern), na.rm = na.rm)
    if (!is.null(pattern))
      args$pattern <- pattern
    selected <- tryCatch(do.call(feature_selection, args), error = function(e)
      stop(src,
      ": ", conditionMessage(e), call. = FALSE))
    if (anyNA(selected) || any(!nzchar(selected)))
      stop(src, " selection produced missing/empty identifiers; check scores, na.rm and pattern.", call. = FALSE)
    base[[src]] <- unique(selected)
    resolved <- rownames(v)
    if (strip_time && !is.null(pattern))
      resolved <- sub(pattern, "", resolved)
    audit[[i]] <- data.frame(source = src, n_vip_rows = nrow(v),
      n_resolved_identifiers = length(unique(resolved)),
      n_collapsed_rows = nrow(v) - length(unique(resolved)), n_selected = length(selected),
      stringsAsFactors = FALSE)
  }
  sets <- list()
  operations <- character()
  used <- list()
  add <- function(name, values, operation, src) {
    sets[[name]] <<- unique(values)
    operations[name] <<- operation
    used[[name]] <<- src
  }
  if ("individual" %in% combinations)
    for (src in sources)
      add(src, base[[src]], "individual", src)
  for (operation in intersect(c("intersection", "union"), combinations)) {
    if (length(sources) < 2L)
      next
    for (size in 2:length(sources)) {
      for (src in utils::combn(sources, size, simplify = FALSE)) {
        fun <- if (operation == "intersection")
          intersect else
            union
        sep <- if (operation == "intersection")
          "_intersect_" else
            "_union_"
        add(paste(src, collapse = sep), Reduce(fun, base[src]), operation, src)
      }
    }
  }
  if (!length(sets))
    stop("The requested sources/combinations produce no candidate sets.", call. = FALSE)
  meta <- data.frame(feature_set = names(sets), n_features = lengths(sets),
    operation = unname(operations),
    row.names = NULL)
  meta$sources <- unname(used)
  structure(list(sets = sets, summary = meta, selection = do.call(rbind, audit),
    settings = list(thr = thr,
    sources = sources, combinations = combinations, cols = cols, aggregator = aggregator,
    strip_time = strip_time, pattern = pattern, na.rm = na.rm)), class = "vip_feature_sets")
}
.vip_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
}
.vip_with_seed <- function(seed, fun) {
  if (is.null(seed))
    return(fun())
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
    seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be NULL or a nonnegative integer.", call. = FALSE)
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed)
    old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed)
    assign(".Random.seed", old, envir = .GlobalEnv) else
      if (exists(".Random.seed",
    envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  fun()
}
#' @export
print.vip_feature_sets <- function(x, ...) {
  print(x$summary, row.names = FALSE)
  invisible(x)
}
.vip_select_strategy <- function(vips, strategy, thr) {
  strategies <- .vip_strategies()
  if (
    !is.character(strategy) ||
    length(strategy) != 1L ||
    is.na(strategy) ||
    !strategy %in% strategies
  )
    stop("Unsupported strategy. Use one of: ", paste(strategies, collapse = ", "), call. = FALSE)
  if (!is.numeric(thr) || is.complex(thr) || length(thr) != 1L || !is.finite(thr) ||
    thr < 0 || thr > 100)
    stop("thr must be one finite percentile between 0 and 100.", call. = FALSE)
  if (!inherits(vips, "tensorpls_vips") || !is.list(vips))
    stop("Single-strategy selection requires a tensorpls_vips object.", call. = FALSE)
  features <- vips$feature_names
  times <- vips$time_names
  H <- vips$ncomp
  for (labels in list(features, times)) {
    if (
      !is.character(labels) ||
      !length(labels) ||
      anyNA(labels) ||
      any(!nzchar(labels)) ||
      anyDuplicated(labels)
    )
      stop("VIP feature/time labels must be unique nonempty strings.", call. = FALSE)
  }
  if (
    !is.numeric(H) ||
    is.complex(H) ||
    length(H) != 1L ||
    !is.finite(H) ||
    H < 1 ||
    H > .Machine$integer.max ||
    H != floor(H)
  )
    stop("VIP ncomp must be one positive integer.", call. = FALSE)
  H <- as.integer(H)
  p <- length(features)
  k <- length(times)
  components <- paste0("Comp", seq_len(H))
  if (!identical(vips$component_names, components))
    stop("VIP component labels do not match ncomp.", call. = FALSE)
  keys <- as.vector(outer(features, times, paste, sep = "_"))
  if (anyDuplicated(keys))
    stop("Ambiguous feature/time labels.", call. = FALSE)
  dims <- list(V = c(p, k, H), VIP2D = c(p * k, H), VIP3Dmodel1 = c(p, H), VIP3Dmodel2 = c(p,
    k), M2_mean = c(p, k))
  labels <- list(V = list(features, times, components), VIP2D = list(keys, components),
    VIP3Dmodel1 = list(features, paste0(components, "_meanTime")), VIP3Dmodel2 = list(features,
      times), M2_mean = list(features, times))
  for (nm in names(dims)) {
    a <- vips[[nm]]
    if (
      !is.numeric(a) ||
      is.complex(a) ||
      is.null(dim(a)) ||
      length(dim(a)) != length(dims[[nm]]) ||
      !all(dim(a) == dims[[nm]]) ||
      any(!is.finite(a)) ||
      !identical(dimnames(a), labels[[nm]])
    )
      stop(nm, " has missing/invalid values, dimensions or labels.", call. = FALSE)
  }
  thresholds <- list()
  if (strategy == "Original")
    selected <- features else {
    sets <- vip_feature_sets(vips, thr = thr, sources = .vip_strategy_sources(strategy))
    selected <- features[features %in% sets$sets[[strategy]]]
    used <- strsplit(strategy, "_intersect_", fixed = TRUE)[[1L]]
    thresholds <- stats::setNames(lapply(used, function(nm) {
      a <- vips[[.vip_source_fields()[[nm]]]]
      stats::setNames(vapply(seq_len(ncol(a)), function(j)
        stats::quantile(a[, j],
        thr / 100, na.rm = TRUE, names = FALSE), numeric(1)), colnames(a))
    }), used)
  }
  list(selected_features = selected, selected_indices = match(selected, features),
    strategy = strategy,
    percentile = thr / 100, thresholds = thresholds)
}

.vip_source_fields <- function() {
  c(VIP2D = "VIP2D", M1 = "VIP3Dmodel1", M2 = "M2_mean", VIP3Dmodel2 = "VIP3Dmodel2")
}
.vip_strategies <- function() {
  c("Original", "VIP2D", "M1", "M2", "VIP3Dmodel2",
    "VIP2D_intersect_M1", "VIP2D_intersect_M2", "M1_intersect_M2",
    "VIP2D_intersect_M1_intersect_M2", "VIP2D_intersect_VIP3Dmodel2")
}
.vip_strategy_sources <- function(strategy) {
  strsplit(strategy, "_intersect_", fixed = TRUE)[[1L]]
}
