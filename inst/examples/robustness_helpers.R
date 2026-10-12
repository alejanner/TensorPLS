# Analysis helpers only: all fitting delegates to the unchanged nested-CV API.
robustness_permutations <- function(Y, R = 100L, seed = 20260923L) {
  stopifnot(is.matrix(Y), nrow(Y) > 1L, length(R) == 1L, R >= 1L, R ==
    as.integer(R))
  TensorPLS:::.vip_with_seed(seed, function()
    lapply(seq_len(R), function(i)
    sample.int(nrow(Y))))
}
robustness_response <- function(Y, index) {
  stopifnot(identical(sort(as.integer(index)), seq_len(nrow(Y))))
  answer <- Y[index, , drop = FALSE]
  dimnames(answer) <- dimnames(Y)  # subject positions stay fixed; only labels move
  answer
}
robustness_predictions <- function(result, ids) {
  p <- result$outer_predictions
  stopifnot(length(unique(p$response)) == 1L, nrow(p) == length(ids), !anyDuplicated(p$subject_id),
    setequal(p$subject_id, ids), identical(as.integer(p$outer_fold),
      as.integer(result$outer_fold_assignments[p$subject_id])))
  p$selected_features <- NULL
  p[match(ids, p$subject_id), , drop = FALSE]
}
robustness_case <- function(X, Y, fold_seed, ncomp = 2:5, permutation = NULL,
  outer_folds = 5L, inner_folds = 5L, vip_thr = 99) {
  original_Y <- Y
  local_Y <- if (is.null(permutation))
    Y else
      robustness_response(Y, permutation)
  ids <- dimnames(X)[[1]]
  folds <- TensorPLS::make_tensorpls_folds(local_Y, outer_folds, inner_folds, fold_seed,
    ids)
  results <- list()
  elapsed <- numeric()
  captured <- character()
  for (strategy in c('Original', 'VIP2D_intersect_VIP3Dmodel2')) {
    start <- proc.time()[['elapsed']]
    r <- withCallingHandlers(TensorPLS::nested_cv_tensorpls(X, local_Y, folds = folds,
      seed = fold_seed, ncomp = ncomp, feature_strategy = strategy, vip_thr = vip_thr,
      strategy_mode = 'fixed', imputation = 'training_mean', input_provenance = 'raw',
      keep_models = FALSE, verbose = FALSE), warning = function(w)
        captured <<- c(captured,
      paste(strategy, conditionMessage(w))))
    stopifnot(identical(Y, original_Y), identical(r$folds, folds))
    # Check that chosen ncomp is the inner-CV optimum, never an outer metric.
    for (f in seq_along(folds$inner)) {
      tuning <- r$inner_tuning_results[[f]]$tuning
      good <- which(is.finite(tuning$selection_score))
      if (length(good) && is.finite(r$selected_ncomp_per_outer_fold[f])) {
        best <- good[order(-tuning$selection_score[good], tuning$ncomp[good],
          good)[1]]
        stopifnot(r$selected_ncomp_per_outer_fold[f] == tuning$ncomp[best])
      }
    }
    r$outer_predictions <- robustness_predictions(r, ids)
    # Retain convergence, tuning, folds, genes and metrics; avoid hundreds of
    # duplicate full-size scaling vectors in these analysis archives.
    r$preprocessing <- NULL
    results[[strategy]] <- r
    elapsed[strategy] <- proc.time()[['elapsed']] - start
  }
  list(results = results, folds = folds, permutation = permutation, fold_seed = fold_seed,
    ncomp = ncomp, elapsed_seconds = elapsed, warnings = captured, labels_restored = identical(Y,
      original_Y))
}
