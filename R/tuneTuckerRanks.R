#' Tune Tucker ranks by masked-entry reconstruction
#' @param X Numeric subjects x features x time tensor, possibly incomplete.
#' @param rank_candidates Nonempty matrix or data.frame with three positive
#'   integer rank columns in subject, feature, time order. Rows must be unique
#'   and satisfy the fitter's dimension and product-of-other-ranks constraints.
#'   NULL generates ranks 1:4 per mode, filtered by these constraints.
#' @param mask_fraction Fraction of eligible observed entries hidden, in (0,1).
#' @param repeats Positive integer number of independent artificial masks.
#' @param seed Nonnegative integer seed; caller RNG state is preserved.
#' @param tucker_control Named controls passed to [fit_tucker_imputer()],
#'   excluding X_train and ranks.
#' @param verbose Report repetition progress.
#' @param candidates Alias for rank_candidates. Supply only one of the two.
#' @details Selects ranks by unsupervised masked reconstruction, not downstream
#'   PLS performance or predictive validation. Within a validation workflow,
#'   supply only the appropriate training subjects. Full-data rank tuning is
#'   suitable for a final descriptive/deployment model, not a claim of
#'   untouched-test predictive performance.
#'
#'   Each repetition hides a random subset of observed entries while retaining
#'   observations for every subject, feature and time. All candidates use the
#'   same mask within a repetition. Selection minimizes mean masked-entry RMSE
#'   among candidates valid in every repetition; exact ties retain input order.
#'   Observed RMSE is a training diagnostic, not the selection criterion.
#'   A fit is valid only if it converges, preserves retained observations and
#'   yields finite masked predictions. Naturally missing entries are not scored.
#'   Fit errors and warnings are retained in the results.
#' @return A tensorpls_tucker_tuning list: selected_ranks (named integer vector,
#'   all NA if none eligible), tuning (candidate summaries), repetition_results
#'   (metrics, validity and diagnostics), masks (linear indices and repetition
#'   seeds), settings, and diagnostics. Summaries use valid repetitions only;
#'   SD is NA with fewer than two. No completed tensors are retained.
#'   results is an alias of the candidate-level tuning table. The default grid
#'   varies subject rank fastest, then feature rank, then time rank. Users may
#'   supply any other valid ranks to [ImputeTensor()] instead of selected_ranks.
#' @export
tune_tucker_ranks <- function(X, rank_candidates = NULL, mask_fraction = 0.10, repeats = 3L,
  seed = 123L, tucker_control = list(), verbose = TRUE, candidates = NULL) {
  .tk_validate(X)
  if (is.complex(X))
    stop("X must be real numeric.", call. = FALSE)
  if (!is.null(candidates)) {
    if (!is.null(rank_candidates))
      stop("Supply only one of candidates and rank_candidates.", call. = FALSE)
    rank_candidates <- candidates
  }
  if (is.null(rank_candidates))
    rank_candidates <- .tucker_rank_grid(dim(X))
  .nested_integer(repeats, "repeats")
  .nested_integer(seed, "seed", 0L)
  .vip_flag(verbose, "verbose")
  if (!is.numeric(mask_fraction) || is.complex(mask_fraction) || length(mask_fraction) !=
    1L || !is.finite(mask_fraction) || mask_fraction <= 0 || mask_fraction >=
    1)
    stop("mask_fraction must be between zero and one.", call. = FALSE)
  if (!(is.matrix(rank_candidates) || is.data.frame(rank_candidates)) || ncol(rank_candidates) !=
    3L || !nrow(rank_candidates) || !all(vapply(as.data.frame(rank_candidates),
    function(z)
      is.numeric(z) && !is.complex(z), logical(1))))
    stop("rank_candidates must be a nonempty numeric matrix/data.frame with three columns.", call. = FALSE)
  ranks <- as.matrix(rank_candidates)
  d <- dim(X)
  if (any(!is.finite(ranks)) || any(ranks < 1 | ranks != floor(ranks)) || any(sweep(ranks,
    2, d, ">")) || any(ranks > .Machine$integer.max))
    stop("Ranks must be positive integers no larger than tensor dimensions.", call. = FALSE)
  if (anyDuplicated(as.data.frame(ranks)))
    stop("Duplicate rank candidates.", call. = FALSE)
  for (j in 1:3)
    if (any(ranks[, j] > ranks[, setdiff(1:3, j)[1]] *
    ranks[, setdiff(1:3, j)[2]]))
    stop("Each Tucker rank must be no larger than the product of the other two.", call. = FALSE)
  storage.mode(ranks) <- "integer"
  colnames(ranks) <- c("rank_subject", "rank_feature", "rank_time")
  candidates <- as.data.frame(ranks)
  rownames(candidates) <- NULL
  allowed <- setdiff(names(formals(fit_tucker_imputer)), c("X_train", "ranks"))
  if (!is.list(tucker_control) || (length(tucker_control) && (is.null(names(tucker_control)) ||
    anyNA(names(tucker_control)) || anyDuplicated(names(tucker_control)) || !all(names(tucker_control) %in%
    allowed))))
    stop("tucker_control must contain named fitter controls, excluding X_train and ranks.", call. = FALSE)
  observed <- which(!is.na(X))
  coords <- arrayInd(observed, d)
  if (!length(observed) || any(vapply(1:3, function(j)
    length(unique(coords[,
    j])) != d[j], logical(1))))
    stop("Every subject, feature and time must have an observed entry.", call. = FALSE)
  protected <- unique(unlist(lapply(1:3, function(j)
    observed[!duplicated(coords[,
    j])])))
  eligible <- setdiff(observed, protected)
  if (!length(eligible))
    stop("No observed entries remain eligible for masking after protection.", call. = FALSE)
  nmask <- max(1, round(mask_fraction * length(eligible)))
  seeds <- as.integer((as.double(seed) + seq_len(repeats) - 1)%%.Machine$integer.max)
  masks <- vector("list", repeats)
  rows <- vector("list", repeats * nrow(ranks))
  row <- 0L
  scalar <- function(z)
    if (is.numeric(z) && length(z) == 1L)
    as.numeric(z) else
      NA_real_
  for (rep in seq_len(repeats)) {
    hidden <- .vip_with_seed(seeds[rep], function()
      sort(eligible[sample.int(length(eligible),
      nmask)]))
    masks[[rep]] <- list(repetition = rep, seed = seeds[rep], indices = hidden)
    masked <- X
    masked[hidden] <- NA_real_
    retained <- setdiff(observed, hidden)
    if (verbose)
      message("Tucker rank tuning: repetition ", rep, "/", repeats)
    for (i in seq_len(nrow(ranks))) {
      diagnostic <- character()
      fit <- NULL
      fit <- tryCatch(withCallingHandlers(.vip_with_seed(seeds[rep], function()
        do.call(fit_tucker_imputer,
        c(list(X_train = masked, ranks = unname(ranks[i, ])), tucker_control))),
        warning = function(w) {
          diagnostic <<- c(diagnostic, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) {
        diagnostic <<- c(diagnostic, conditionMessage(e))
        NULL
      })
      converged <- !is.null(fit) && isTRUE(fit$convergence$converged)
      valid <- FALSE
      rmse <- mae <- NA_real_
      unchanged <- NA
      if (!is.null(fit)) {
        completed <- fit$X_train_imputed
        if (is.numeric(completed) && identical(dim(completed), d)) {
          unchanged <- identical(as.numeric(completed[retained]), as.numeric(X[retained]))
          if (!unchanged)
          diagnostic <- c(diagnostic, "Observed non-masked entries changed.")
          valid <- converged && unchanged && all(is.finite(completed[hidden]))
          if (valid) {
          error <- completed[hidden] - X[hidden]
          rmse <- sqrt(mean(error^2))
          mae <- mean(abs(error))
          valid <- is.finite(rmse) && is.finite(mae)
          } else
            if (any(!is.finite(completed[hidden])))
          diagnostic <- c(diagnostic, "Nonfinite masked predictions.")
        } else
          diagnostic <- c(diagnostic, "Invalid completed tensor.")
        if (!converged)
          diagnostic <- c(diagnostic, "Tucker fit did not converge.")
      }
      if (!valid) {
        rmse <- mae <- NA_real_
      }
      row <- row + 1L
      rows[[row]] <- cbind(data.frame(candidate = i, repetition = rep), candidates[i,
        , drop = FALSE], data.frame(RMSE_masked = rmse, MAE_masked = mae, valid = valid,
        converged = converged, iterations = scalar(fit$convergence$iterations),
        final_delta = scalar(fit$convergence$final_delta),
          observed_RMSE = scalar(fit$convergence$observed_RMSE),
        observed_unchanged = unchanged, diagnostic = paste(unique(diagnostic),
          collapse = " | ")))
    }
  }
  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  tuning <- candidates
  summaries <- lapply(seq_len(nrow(ranks)), function(i) {
    z <- results[results$candidate == i, , drop = FALSE]
    good <- z[z$valid, , drop = FALSE]
    data.frame(mean_RMSE = if (nrow(good))
      mean(good$RMSE_masked) else
        NA_real_, sd_RMSE = stats::sd(good$RMSE_masked),
        mean_MAE = if (nrow(good))
      mean(good$MAE_masked) else
        NA_real_, sd_MAE = stats::sd(good$MAE_masked),
        valid = all(z$valid), valid_repetitions = nrow(good))
  })
  tuning <- cbind(tuning, do.call(rbind, summaries))
  good <- which(tuning$valid)
  selected <- if (length(good))
    unname(ranks[good[which.min(tuning$mean_RMSE[good])], ]) else
      rep(NA_integer_, 3)
  names(selected) <- colnames(ranks)
  structure(list(selected_ranks = selected, tuning = tuning, repetition_results = results,
    masks = masks, settings = list(mask_fraction = mask_fraction, repeats = repeats,
      seed = seed, tucker_control = tucker_control, dimensions = d,
        observed_count = length(observed),
      protected_indices = protected, eligible_count = length(eligible), masked_count = nmask),
    diagnostics = list(message = if (length(good))
      NULL else
        "No candidate completed every repetition successfully.",
        invalid_fits = sum(!results$valid)),
    results = tuning), class = "tensorpls_tucker_tuning")
}
.tucker_rank_grid <- function(dimensions) {
  grid <- expand.grid(rank_subject = 1:4, rank_feature = 1:4,
    rank_time = 1:4, KEEP.OUT.ATTRS = FALSE)
  valid <- apply(grid, 1L, function(ranks) {
    all(ranks <= dimensions) && all(vapply(1:3, function(mode)
      ranks[mode] <=
      prod(ranks[-mode]), logical(1)))
  })
  grid <- grid[valid, , drop = FALSE]
  rownames(grid) <- NULL
  grid
}
