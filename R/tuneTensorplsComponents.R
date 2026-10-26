#' Choose components by full-development-data cross-validation
#'
#' @param X Numeric subjects x features x time tensor with unique axis labels.
#' @param Y Classification response accepted by [nested_cv_tensorpls()].
#' @param ncomp Candidate integers from 1 to 10; duplicates are removed.
#' @param feature_strategy One fixed strategy accepted by [vip_feature_sets()],
#'   including M1 (all-time mean), M2 (cumulative-path mean), and their
#'   ordered intersections. VIP3Dmodel2 uses cumulative-path sum selection.
#' @param vip_thr One percentile in [0,100]; ignored for Original.
#' @param cv_folds Number of stratified folds, at least two.
#' @param seed Nonnegative integer seed; caller RNG state is preserved.
#' @param folds Optional numeric fold vector named in exact X subject order.
#'   Labels must be consecutive integers 1 through cv_folds. Every training
#'   and validation split must contain all classes.
#' @param imputation Training-only imputation: none, training_mean, or tucker.
#' @param tucker_control Named fixed Tucker controls, as in [nested_cv_tensorpls()].
#' @param input_provenance raw, unknown, or preimputed. Non-raw input warns
#'   because earlier preprocessing may have used validation information.
#' @param selection_metric Q2Y (maximize pooled Q2Y) or BER (minimize pooled BER).
#' @param verbose Report development progress.
#' @details This selects H for subsequent final full-data fitting/interpretation.
#'   CV scores are tuning diagnostics, not unbiased predictive performance;
#'   use [nested_cv_tensorpls()] for predictive assessment. Feature selection
#'   and imputation are fitted separately inside each training fold. Q2Y uses
#'   pooled held-out errors and the corresponding training response means.
#'   Score ties prefer fewer components, then candidate order. The shared
#'   development engine also refits the chosen pipeline on all supplied data.
#'   A candidate must complete every validation fold without diagnostics.
#' @return A tensorpls_component_tuning list with selected_ncomp, tuning
#'   (candidate scores/configurations), fold_results (shared-engine diagnostics),
#'   folds, settings, diagnostics, and refit (the internal full-data pipeline).
#'   selected_ncomp is NA and refit is NULL if development or the final refit
#'   fails; diagnostics$error explains the failure. Scores remain available.
#' @export
tune_tensorpls_components <- function(X, Y, ncomp = 1:10, feature_strategy = "Original",
  vip_thr = 99, cv_folds = 5L, seed = 123L, folds = NULL, imputation = c("none",
    "training_mean", "tucker"), tucker_control = list(ranks = c(4L, 4L, 3L)),
  input_provenance = "raw", selection_metric = c("Q2Y", "BER"), verbose = TRUE) {
  imputation <- match.arg(imputation)
  input_provenance <- match.arg(input_provenance, c("raw", "unknown", "preimputed"))
  selection_metric <- match.arg(selection_metric)
  .nested_integer(seed, "seed", 0L)
  .nested_integer(cv_folds, "cv_folds", 2L)
  .vip_flag(verbose, "verbose")
  if (!is.numeric(X) || is.complex(X) || length(dim(X)) != 3L || any(dim(X) < 1L) ||
    any(is.infinite(X)))
    stop("X must be a real numeric three-way tensor without infinite values.", call. = FALSE)
  ids <- dimnames(X)[[1]]
  features <- dimnames(X)[[2]]
  times <- dimnames(X)[[3]]
  .nested_ids(ids, dim(X)[1])
  .nested_ids(features, dim(X)[2])
  .nested_ids(times, dim(X)[3])
  encoded <- .nested_response(Y)
  if (nrow(encoded$Y) != dim(X)[1])
    stop("X and Y subject counts differ.", call. = FALSE)
  yids <- if (is.null(dim(Y)))
    names(Y) else
      dimnames(Y)[[1]]
  if (!is.null(yids) && !identical(yids, ids))
    stop("Y subjects must match X in identical order.", call. = FALSE)
  if (length(dim(Y)) == 3L && (dim(Y)[3] != length(times) || (!is.null(dimnames(Y)[[3]]) &&
    !identical(dimnames(Y)[[3]], times))))
    stop("Y time axis must match X.", call. = FALSE)
  strategies <- .vip_strategies()
  if (!is.character(feature_strategy) || length(feature_strategy) != 1L || is.na(feature_strategy) ||
    !feature_strategy %in% strategies)
    stop("feature_strategy must be one supported strategy.", call. = FALSE)
  if (!is.numeric(vip_thr) || is.complex(vip_thr) || length(vip_thr) != 1L || !is.finite(vip_thr) ||
    vip_thr < 0 || vip_thr > 100)
    stop("vip_thr must be one percentile in [0,100].", call. = FALSE)
  if (!is.numeric(ncomp) || is.complex(ncomp) || !length(ncomp) || any(!is.finite(ncomp)) ||
    any(ncomp < 1 | ncomp > 10 | ncomp != floor(ncomp)))
    stop("ncomp must contain integers from 1 to 10.", call. = FALSE)
  if (is.null(folds))
    folds <- .vip_with_seed(seed, function()
      stats::setNames(.nested_stratify(encoded$class,
      cv_folds), ids))
  if (
    !is.numeric(folds) ||
    is.complex(folds) ||
    length(folds) != length(ids) ||
    !identical(names(folds), ids) ||
    any(!is.finite(folds)) ||
    any(folds < 1 | folds > cv_folds | folds != floor(folds)) ||
    !identical(sort(unique(as.integer(folds))), seq_len(cv_folds))
  )
    stop("folds must be numbered 1 through cv_folds in exact subject order.", call. = FALSE)
  for (f in seq_len(cv_folds)) {
    if (
      !setequal(encoded$class[folds != f], encoded$levels) ||
      !setequal(encoded$class[folds == f], encoded$levels)
    )
      stop("Every training and validation fold must contain all classes.", call. = FALSE)
  }
  method <- imputation
  if (method == "tucker") {
    if (!is.list(tucker_control) || is.null(names(tucker_control)) || anyDuplicated(names(tucker_control)) ||
      !all(names(tucker_control) %in% c("ranks", "tolerance", "max_iter", "seed",
        "projection_tolerance",
        "initialization", "ridge")))
      stop("Invalid tucker_control; use named fitting controls, with fixed ranks.", call. = FALSE)
    if (is.null(tucker_control$seed))
      tucker_control$seed <- seed
    imputation <- list(method = "tucker", control = tucker_control)
  }
  if (method == "none" && anyNA(X))
    stop("Missing X requires training_mean or tucker imputation.", call. = FALSE)
  Xm <- if (method == "tucker")
    X else
      matrix(X, dim(X)[1])
  if (method != "tucker")
    colnames(Xm) <- as.vector(outer(features, times, paste, sep = "_"))
  note <- character()
  if (input_provenance != "raw") {
    note <- "Input provenance is not raw: prior preprocessing may leak across subjects."
    warning(note, call. = FALSE)
  }
  configurations <- data.frame(strategy = feature_strategy, thr = if (feature_strategy ==
    "Original")
    NA_real_ else
      vip_thr, ncomp = sort(unique(as.integer(ncomp))), stringsAsFactors = FALSE)
  if (verbose)
    message("Development CV: evaluating components, then refitting the chosen pipeline.")
  development <- .nested_develop(Xm, encoded$Y, encoded$class, encoded[c("levels", "binary")],
    features, times, folds, configurations, imputation, selection_metric)
  structure(
    list(
      selected_ncomp = if (is.null(development$model))
        NA_integer_ else
          as.integer(development$model$config$ncomp),
      tuning = development$tuning,
      fold_results = development$inner_results,
      folds = folds,
      settings = list(
        feature_strategy = feature_strategy,
        vip_thr = configurations$thr[1],
        ncomp = configurations$ncomp,
        cv_folds = cv_folds,
        seed = seed,
        imputation = method,
        tucker_control = if (method == "tucker") tucker_control else NULL,
        input_provenance = input_provenance,
        selection_metric = selection_metric,
        score_interpretation = "Development tuning diagnostic, not unbiased predictive performance"
      ),
      diagnostics = list(
        error = development$error,
        provenance = note,
        imputation = development$imputation_diagnostics
      ),
      refit = development$model
    ),
    class = "tensorpls_component_tuning"
  )
}
