#' Stratified subject folds for nested TensorPLS evaluation
#' @param Y Class factor/vector, binary response matrix, one-hot matrix, or
#' a three-way response array whose class coding is identical at all times.
#' @param outer_folds,inner_folds Numbers of stratified folds (at least two).
#' @param seed Integer random seed. The caller's RNG state is preserved.
#' @param subject_ids Unique identifiers, in Y row order.
#' @return A list with outer assignments named by subject, and named inner
#' assignments containing only the corresponding OUTER TRAIN subjects.
#' @export
make_tensorpls_folds <- function(
  Y, outer_folds = 5L, inner_folds = 5L, seed = 123L, subject_ids = NULL
) {
  encoded <- .nested_response(Y)
  n <- nrow(encoded$Y)
  if (is.null(subject_ids))
    subject_ids <- as.character(seq_len(n))
  .nested_ids(subject_ids, n)
  .nested_integer(outer_folds, "outer_folds", 2L)
  .nested_integer(inner_folds, "inner_folds", 2L)
  .nested_integer(seed, "seed", 0L)
  .vip_with_seed(seed, function() {
    outer <- .nested_stratify(encoded$class, outer_folds)
    names(outer) <- subject_ids
    inner <- lapply(seq_len(outer_folds), function(f) {
      train <- which(outer != f)
      stats::setNames(.nested_stratify(encoded$class[train], inner_folds), subject_ids[train])
    })
    names(inner) <- as.character(seq_len(outer_folds))
    list(outer = outer, inner = inner, seed = seed)
  })
}
.nested_integer <- function(x, name, min = 1L) {
  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    !is.finite(x) ||
    x < min ||
    x > .Machine$integer.max ||
    x != as.integer(x)
  )
    stop(name, " must be one integer >= ", min, ".", call. = FALSE)
}
.nested_ids <- function(ids, n) {
  if (
    !is.character(ids) ||
    length(ids) != n ||
    anyNA(ids) ||
    any(!nzchar(ids)) ||
    anyDuplicated(ids)
  )
    stop("Subject identifiers must be unique, nonempty strings.", call. = FALSE)
}
.nested_stratify <- function(classes, k) {
  if (any(table(classes) < k))
    stop("Each class needs at least as many subjects as folds.")
  folds <- integer(length(classes))
  for (level in unique(classes)) {
    idx <- which(classes == level)
    # sample.int avoids sample(x) interpreting a singleton numeric x as 1:x.
    idx <- idx[sample.int(length(idx))]
    folds[idx] <- rep_len(sample.int(k), length(idx))
  }
  folds
}
.nested_response <- function(Y) {
  if (is.array(Y) && length(dim(Y)) == 3L) {
    first <- matrix(Y[, , 1], dim(Y)[1], dim(Y)[2], dimnames = dimnames(Y)[1:2])
    for (t in seq_len(dim(Y)[3]))
      if (!isTRUE(all.equal(
        unname(matrix(Y[, , t], dim(Y)[1], dim(Y)[2])), unname(first)
      )))
        stop("Nested PLS-DA requires time-invariant class coding; class-2 response tensors are not supported.")
    Y <- first
  }
  if (is.matrix(Y)) {
    if (
      !is.numeric(Y) ||
      any(!is.finite(Y)) ||
      !all(Y %in% c(0, 1))
    )
      stop("Matrix Y must be finite binary or one-hot class coding.")
    if (ncol(Y) == 1L) {
      classes <- as.character(Y[, 1])
      levels <- c("0", "1")
      binary <- TRUE
      if (!all(levels %in% classes))
        stop("Both classes must be present.")
      if (is.null(colnames(Y)))
        colnames(Y) <- "Class1"
    } else {
      if (any(rowSums(Y) != 1))
        stop("Multiresponse Y must be one-hot class coding.")
      levels <- colnames(Y)
      if (is.null(levels))
        levels <- paste0("Class", seq_len(ncol(Y)))
      if (
        anyNA(levels) ||
        any(!nzchar(levels)) ||
        anyDuplicated(levels)
      )
        stop("Class names must be unique.")
      colnames(Y) <- levels
      classes <- levels[max.col(Y, ties.method = "first")]
      if (any(colSums(Y) == 0))
        stop("Every declared class must be observed.")
      binary <- FALSE
    }
  } else {
    if (
      !is.atomic(Y) ||
      !is.null(dim(Y)) ||
      anyNA(Y)
    )
      stop("Y must contain class labels without missing values.")
    f <- if (is.factor(Y))
      droplevels(Y) else factor(Y)
    levels <- levels(f)
    classes <- as.character(f)
    if (length(levels) < 2L)
      stop("At least two classes are required.")
    binary <- length(levels) == 2L
    Y <- if (binary)
      matrix(as.numeric(f == levels[2]), ncol = 1L, dimnames = list(NULL, levels[2])) else
      vapply(levels, function(z) as.numeric(classes == z), numeric(length(classes)))
  }
  list(Y = Y, class = classes, levels = levels, binary = binary)
}
.nested_validate_folds <- function(folds, ids, classes) {
  a <- folds$outer
  if (
    !is.numeric(a) ||
    !identical(names(a), ids) ||
    any(!is.finite(a)) ||
    any(a < 1 | a > .Machine$integer.max) ||
    any(a != as.integer(a)) ||
    !identical(sort(unique(as.integer(a))), seq_len(max(a))) ||
    max(a) < 2L
  )
    stop("Invalid outer fold assignments or subject order.")
  if (
    length(folds$inner) != max(a) ||
    !identical(names(folds$inner), as.character(seq_len(max(a))))
  )
    stop("An inner fold assignment is required for every outer fold.")
  for (f in seq_len(max(a))) {
    tr <- which(a != f)
    b <- folds$inner[[as.character(f)]]
    if (
      !is.numeric(b) ||
      !identical(names(b), ids[tr]) ||
      any(!is.finite(b)) ||
      any(b < 1 | b > .Machine$integer.max) ||
      any(b != as.integer(b)) ||
      !identical(sort(unique(as.integer(b))), seq_len(max(b))) ||
      max(b) < 2L
    )
      stop("Inner folds must contain exactly the corresponding OUTER TRAIN subjects, in order.")
    for (j in seq_len(max(b)))
      if (!setequal(classes[tr[b != j]], unique(classes)))
        stop("Every INNER TRAIN split must contain all classes.")
  }
  invisible(TRUE)
}
.nested_imputer <- function(X, method) {
  if (method == "none" && anyNA(X))
    stop("Missing X needs imputation='training_mean' or an externally justified complete input.")
  means <- colMeans(X, na.rm = TRUE)
  all_missing <- !is.finite(means)
  means[all_missing] <- 0
  list(
    method = method,
    means = means,
    all_missing = all_missing,
    convergence = list(
      converged = TRUE,
      iterations = 0L,
      final_delta = 0,
      tolerance = NA_real_,
      method = method,
      iterative = FALSE
    )
  )
}
.nested_impute <- function(object, X) {
  if (
    ncol(X) != length(object$means) ||
    any(is.infinite(X))
  )
    stop("Invalid imputation prediction dimensions/values.")
  if (object$method == "none" && anyNA(X))
    stop("Missing held-out X with imputation='none'.")
  if (anyNA(X)) {
    idx <- which(is.na(X), arr.ind = TRUE)
    X[idx] <- object$means[idx[, 2]]
  }
  # A column entirely missing during training is unusable, including at prediction.
  X[, object$all_missing] <- 0
  X
}
.nested_vips <- function(fit, features, times) {
  p <- length(features)
  k <- length(times)
  h <- fit$ncomp
  vip <- fit$VIP
  rownames(vip) <- as.vector(outer(features, times, paste, sep = "_"))
  colnames(vip) <- paste0("Comp", seq_len(h))
  a <- array(vip, c(p, k, h))
  m1 <- matrix(
    apply(a, c(1, 3), mean), p, h,
    dimnames = list(features, paste0("Comp", seq_len(h), "_meanTime"))
  )
  m2_sum <- matrix(apply(a, c(1, 2), sum), p, k, dimnames = list(features, times))
  m2_mean <- matrix(apply(a, c(1, 2), mean), p, k, dimnames = list(features, times))
  list(VIP2D = vip, VIP3Dmodel1 = m1, VIP3Dmodel2 = m2_sum, M2_mean = m2_mean)
}
.nested_rows <- function(X, rows) {
  if (length(dim(X)) == 3L)
    X[rows, , , drop = FALSE] else X[rows, , drop = FALSE]
}
.nested_fit_pipeline <- function(X, Y, features, times, config, imputation) {
  if (is.list(imputation) && identical(imputation$method, "tucker")) {
    if (length(dim(X)) != 3L)
      stop("Tucker training input must remain three-dimensional.")
    key <- paste(dimnames(X)[[1]], collapse = "\r")
    cache <- imputation$cache
    if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE))
      imputer <- get(key, envir = cache) else {
      imputer <- do.call(fit_tucker_imputer, c(list(X_train = X), imputation$control))
      if (!is.null(cache))
        assign(key, imputer, envir = cache)
    }
    if (!isTRUE(imputer$convergence$converged))
      stop("Tucker training imputer did not converge.")
    # Unfold only AFTER the three-way training imputation is complete.
    Xi <- matrix(imputer$X_train_imputed, dim(X)[1])
    colnames(Xi) <- as.vector(outer(features, times, paste, sep = "_"))
  } else {
    imputer <- .nested_imputer(X, imputation)
    Xi <- .nested_impute(imputer, X)
  }
  if (inherits(imputer, "tucker_imputer"))
    imputer$X_train_imputed <- NULL
  vip <- NULL
  if (config$strategy == "Original")
    selected <- features else {
    base <- .tensorpls_fit(Xi, Y, config$ncomp)
    if (!all(base$convergence$converged))
      stop("VIP fit did not converge.")
    vip <- .nested_vips(base, features, times)
    sets <- vip_feature_sets(vip, thr = config$thr, sources = .vip_strategy_sources(config$strategy),
      combinations = c("individual", "intersection"))
    selected <- features[features %in% sets$sets[[config$strategy]]]
  }
  if (length(selected) < config$ncomp)
    stop("Selected genes fewer than requested components.")
  columns <- which(rep(features, length(times)) %in% selected)
  fit <- .tensorpls_fit(Xi[, columns, drop = FALSE], Y, config$ncomp)
  if (!all(fit$convergence$converged))
    stop("Candidate PLS fit did not converge.")
  list(
    fit = fit,
    imputer = imputer,
    columns = columns,
    selected_features = selected,
    config = config,
    vip = vip,
    preprocessing = list(imputer = imputer, scaling = fit$preprocessing),
    convergence = list(imputation = imputer$convergence, pls = fit$convergence)
  )
}
.nested_predict_pipeline <- function(object, X) {
  if (inherits(object$imputer, "tucker_imputer")) {
    # Frozen training core/factors; subjects are projected independently in 3D.
    completed <- predict_tucker_imputer(object$imputer, X)
    Xi <- matrix(completed, dim(completed)[1])
  } else Xi <- .nested_impute(object$imputer, X)
  .tensorpls_predict(object$fit, Xi[, object$columns, drop = FALSE])
}
.nested_classify <- function(pred, encoding) {
  if (encoding$binary)
    encoding$levels[1L + as.integer(pred[, 1] >= 0.5)] else
    encoding$levels[max.col(pred, ties.method = "first")]
}
.nested_metrics <- function(Y, pred, baseline, truth, encoding) {
  if (any(!is.finite(pred)))
    return(list(
      Q2Y = NA_real_,
      Q2Y_per_response = rep(NA_real_, ncol(Y)),
      PRESS = rep(NA_real_, ncol(Y)),
      denominator = rep(NA_real_, ncol(Y)),
      classification = c(
        accuracy = NA_real_,
        balanced_accuracy = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        BER = NA_real_,
        AUC = NA_real_
      )
    ))
  press <- colSums((Y - pred)^2)
  denominator <- colSums((Y - baseline)^2)
  q2 <- ifelse(denominator > 0, 1 - press / denominator, NA_real_)
  assigned <- .nested_classify(pred, encoding)
  lev <- encoding$levels
  recall <- specificity <- numeric(length(lev))
  for (i in seq_along(lev)) {
    pos <- truth == lev[i]
    recall[i] <- if (any(pos))
      mean(assigned[pos] == lev[i]) else NA_real_
    specificity[i] <- if (any(!pos))
      mean(assigned[!pos] != lev[i]) else NA_real_
  }
  binary <- length(lev) == 2L
  auc <- NA_real_
  if (binary) {
    pos <- truth == lev[2]
    np <- sum(pos)
    nn <- sum(!pos)
    score <- if (encoding$binary)
      pred[, 1] else pred[, 2] - pred[, 1]
    if (np && nn)
      auc <- (sum(rank(score, ties.method = "average")[pos]) - np * (np + 1) / 2) / (np * nn)
  }
  list(
    Q2Y = mean(q2),
    Q2Y_per_response = stats::setNames(q2, colnames(Y)),
    PRESS = stats::setNames(press, colnames(Y)),
    denominator = stats::setNames(denominator, colnames(Y)),
    classification = c(
      accuracy = mean(assigned == truth),
      balanced_accuracy = mean(recall),
      sensitivity = if (binary) recall[2] else mean(recall),
      specificity = if (binary) specificity[2] else mean(specificity),
      BER = 1 - mean(recall),
      AUC = auc
    )
  )
}
# Develop exclusively on OUTER TRAIN. Neither outer-test X nor Y is an argument.
.nested_develop <- function(X, Y, classes, encoding, features, times, inner, configurations, imputation,
  selection_metric) {
  if (is.list(imputation) && identical(imputation$method, "tucker"))
    imputation$cache <- new.env(parent = emptyenv())  # this OUTER TRAIN only; no Y-dependent reuse
  rows <- list()
  scores <- rep(NA_real_, nrow(configurations))
  row <- 0L
  for (c in seq_len(nrow(configurations))) {
    config <- configurations[c, , drop = FALSE]
    pred <- baseline <- matrix(NA_real_, nrow(Y), ncol(Y))
    valid <- TRUE
    for (j in seq_len(max(inner))) {
      tr <- which(inner != j)
      validation <- which(inner == j)
      diagnostics <- character()
      fit <- tryCatch(
        withCallingHandlers(
          .nested_fit_pipeline(
            .nested_rows(X, tr), Y[tr, , drop = FALSE],
            features, times, config, imputation
          ),
          warning = function(w) {
            diagnostics <<- c(diagnostics, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          diagnostics <<- c(diagnostics, conditionMessage(e))
          NULL
        }
      )
      met <- NULL
      if (!is.null(fit)) {
        pr <- tryCatch(
          .nested_predict_pipeline(fit, .nested_rows(X, validation)),
          error = function(e) {
            diagnostics <<- c(diagnostics, conditionMessage(e))
            NULL
          }
        )
        if (!is.null(pr)) {
          pred[validation, ] <- pr
          baseline[validation, ] <- matrix(
            rep(colMeans(Y[tr, , drop = FALSE]), each = length(validation)),
            length(validation), ncol(Y)
          )
          met <- .nested_metrics(
            Y[validation, , drop = FALSE], pr,
            baseline[validation, , drop = FALSE], classes[validation], encoding
          )
        }
      }
      if (
        is.null(met) ||
        length(diagnostics)
      )
        valid <- FALSE
      row <- row + 1L
      rows[[row]] <- data.frame(
        candidate = c,
        inner_fold = j,
        strategy = config$strategy,
        thr = config$thr,
        ncomp = config$ncomp,
        n_inner_train = length(tr),
        n_inner_validation = length(validation),
        n_features = if (is.null(fit)) NA_integer_ else length(fit$selected_features),
        Q2Y = if (is.null(met)) NA_real_ else met$Q2Y,
        BER = if (is.null(met)) NA_real_ else unname(met$classification["BER"]),
        status = if (is.null(met)) "failed" else
          if (length(diagnostics)) "warning" else "ok",
        diagnostic = paste(diagnostics, collapse = " | ")
      )
    }
    if (valid) {
      met <- .nested_metrics(Y, pred, baseline, classes, encoding)
      scores[c] <- if (selection_metric == "Q2Y")
        met$Q2Y else -unname(met$classification["BER"])
    }
  }
  imputation_diagnostics <- if (is.list(imputation) && identical(imputation$method, "tucker"))
    lapply(
      as.list(imputation$cache),
      function(z) list(training_ids = z$training_ids, convergence = z$convergence)
    ) else NULL
  results <- do.call(rbind, rows)
  tuning <- cbind(configurations, selection_score = scores)
  valid <- which(is.finite(scores))
  if (!length(valid))
    return(list(
      model = NULL,
      inner_results = results,
      tuning = tuning,
      error = "No candidate completed every INNER VALIDATION fold successfully.",
      imputation_diagnostics = imputation_diagnostics
    ))
  # Deterministic ties: fewer components, then declared candidate order.
  chosen <- valid[order(-scores[valid], configurations$ncomp[valid], valid)[1]]
  final <- tryCatch(
    .nested_fit_pipeline(
      X, Y, features, times, configurations[chosen, , drop = FALSE], imputation
    ),
    error = function(e) e
  )
  if (inherits(final, "error"))
    return(list(
      model = NULL,
      inner_results = results,
      tuning = tuning,
      error = conditionMessage(final),
      imputation_diagnostics = imputation_diagnostics
    ))
  list(
    model = final, inner_results = results, tuning = tuning, error = NULL,
    imputation_diagnostics = imputation_diagnostics
  )
}
#' Nested cross-validation for the full PLS-DA feature-selection procedure
#'
#' @param X Raw numeric subjects x features x time tensor, with unique dimnames.
#' @param Y Class coding accepted by [make_tensorpls_folds()].
#' @param outer_folds,inner_folds Stratified fold counts, default five each.
#' @param seed Reproducible fold seed, default 123; preserves caller RNG state.
#' @param feature_strategy One strategy for fixed mode, or candidate strategies
#' for tuned mode. Original skips VIP selection. Names follow [vip_feature_sets()],
#' including M1 (all-time mean), M2 (cumulative-path mean), and ordered
#' intersections. VIP3Dmodel2 uses cumulative-path sum selection.
#' @param vip_thr Fixed percentile or candidate percentiles for tuned mode.
#' @param ncomp Candidate component counts, selected only in inner CV.
#' @param strategy_mode Fixed strategy/threshold, or tune both inside inner CV.
#' @param folds Optional reusable object from [make_tensorpls_folds()].
#' @param tucker_control Named controls for fit_tucker_imputer; ranks are fixed, not tuned.
#' @param imputation tucker fits a training-only 3D imputer and projects unseen subjects.
#' none requires complete X; training_mean estimates a separate
#' feature-time mean from each training split, never from held-out subjects.
#' @param input_provenance raw asserts no previous across-subject learned
#' preprocessing. unknown/preimputed explicitly prevent an inductive claim.
#' @param selection_metric Maximize pooled inner Q2Y or minimize pooled inner BER.
#' @param keep_models Store outer-training models and their VIPs for auditing.
#' @param verbose Print progress for each outer fold.
#' @return A nested_tensorpls object with held-out predictions, shared fold
#' assignments, inner tuning, selected configurations/features, preprocessing,
#' convergence, diagnostics, classification metrics and predictive Q2Y.
#' @details OUTER TRAIN is used for model development; INNER TRAIN fits candidate
#' pipelines; INNER VALIDATION selects configurations; OUTER TEST is used only
#' for final prediction/scoring. Stratification uses labels to define folds,
#' before fitting, and is the only pre-split use of labels.
#'
#' Q2Y=1-PRESS/denominator per response, then averaged. The denominator sums
#' squared differences from each subject's OUTER TRAIN response mean. Repeated
#' time-invariant Y is collapsed to one response block; its identical copies
#' would give the same averaged Q2Y. Multiclass Y is one-hot and equally weighted
#' by response. Binary classification uses threshold 0.5 (positive=second class);
#' one-hot classification uses argmax. AUC is binary only. No outer R2Y is claimed.
#'
#' Feature identities are relearned in every training split, even for a fixed
#' strategy. Failed candidates are ineligible unless all inner folds succeed.
#' Failed outer folds remain in the result and make aggregate metrics NA.
#' Tucker uses frozen training factors and observed-entry least-squares projection.
#' Non-converged Tucker fits are rejected. Training fits may be cached across
#' component candidates within the identical training split; no validation/test
#' data or response labels enter this cache.
#' Mean imputation is a different pipeline. Any prior imputation, scaling,
#' normalization or filtering learned using held-out subjects invalidates the
#' inductive interpretation, regardless of this function's internal safeguards.
#' @export
nested_cv_tensorpls <- function(
  X, Y, outer_folds = 5L, inner_folds = 5L, seed = 123L,
  feature_strategy = "VIP2D_intersect_VIP3Dmodel2", vip_thr = 99, ncomp = 2:5,
  strategy_mode = c("fixed", "tune"), folds = NULL,
  imputation = c("none", "training_mean", "tucker"),
  input_provenance = c("unknown", "raw", "preimputed"),
  selection_metric = c("Q2Y", "BER"), keep_models = TRUE,
  verbose = TRUE, tucker_control = list(ranks = c(4L, 4L, 3L))
) {
  strategy_mode <- match.arg(strategy_mode)
  imputation <- match.arg(imputation)
  input_provenance <- match.arg(input_provenance)
  selection_metric <- match.arg(selection_metric)
  .nested_integer(seed, "seed", 0L)
  imputation_method <- imputation
  if (imputation == "tucker") {
    if (
      !is.list(tucker_control) ||
      is.null(names(tucker_control)) ||
      anyDuplicated(names(tucker_control)) ||
      !all(names(tucker_control) %in% c(
        "ranks", "tolerance", "max_iter", "seed", "projection_tolerance",
        "initialization", "ridge"
      ))
    )
      stop("Invalid tucker_control; use named fitting controls, with fixed ranks.")
    if (is.null(tucker_control$seed))
      tucker_control$seed <- seed
    imputation <- list(method = "tucker", control = tucker_control)
  }
  if (
    !is.numeric(X) ||
    length(dim(X)) != 3L ||
    any(is.infinite(X))
  )
    stop("X must be a numeric three-way tensor without infinite values.")
  ids <- dimnames(X)[[1]]
  features <- dimnames(X)[[2]]
  times <- dimnames(X)[[3]]
  .nested_ids(ids, dim(X)[1])
  .nested_ids(features, dim(X)[2])
  .nested_ids(times, dim(X)[3])
  encoded <- .nested_response(Y)
  Ym <- encoded$Y
  classes <- encoded$class
  encoding <- encoded[c("levels", "binary")]
  if (nrow(Ym) != dim(X)[1])
    stop("X and Y subject counts differ.")
  yids <- if (is.null(dim(Y)))
    names(Y) else dimnames(Y)[[1]]
  if (!is.null(yids) && !identical(yids, ids))
    stop("Y subjects must match X in identical order.")
  if (
    length(dim(Y)) == 3L && (
      dim(Y)[3] != length(times) || (
        !is.null(dimnames(Y)[[3]]) && !identical(dimnames(Y)[[3]], times)
      )
    )
  )
    stop("Y time axis must match X.")
  strategies <- .vip_strategies()
  if (
    !is.character(feature_strategy) ||
    !length(feature_strategy) ||
    anyNA(feature_strategy) ||
    anyDuplicated(feature_strategy) ||
    !all(feature_strategy %in% strategies)
  )
    stop("Unknown or duplicate feature_strategy.")
  if (
    !is.numeric(vip_thr) ||
    !length(vip_thr) ||
    any(!is.finite(vip_thr)) ||
    any(vip_thr < 0 | vip_thr > 100)
  )
    stop("vip_thr must contain percentiles in [0,100].")
  if (
    strategy_mode == "fixed" && (length(feature_strategy) != 1L ||
    length(vip_thr) != 1L)
  )
    stop("Fixed mode requires one strategy and one threshold.")
  if (
    !is.numeric(ncomp) ||
    !length(ncomp) ||
    any(!is.finite(ncomp)) ||
    any(ncomp < 1 | ncomp > 10 | ncomp != as.integer(ncomp))
  )
    stop("ncomp must contain integers from 1 to 10.")
  .vip_flag(keep_models, "keep_models")
  .vip_flag(verbose, "verbose")
  if (is.null(folds))
    folds <- make_tensorpls_folds(Y, outer_folds, inner_folds, seed, ids)
  .nested_validate_folds(folds, ids, classes)
  configurations <- expand.grid(
    strategy = feature_strategy, thr = unique(vip_thr),
    ncomp = sort(unique(ncomp)), stringsAsFactors = FALSE
  )
  configurations$thr[configurations$strategy == "Original"] <- NA_real_
  configurations <- unique(configurations)
  rownames(configurations) <- NULL
  if (imputation_method == "tucker")
    Xm <- X else {
    Xm <- matrix(X, dim(X)[1])
    colnames(Xm) <- as.vector(outer(features, times, paste, sep = "_"))
  }
  if (imputation_method == "none" && anyNA(Xm))
    stop("Raw missing X requires imputation='training_mean'; global Tucker imputation is not leakage-safe.")
  diagnostics <- character()
  if (input_provenance != "raw") {
    diagnostics <- "Input provenance is not raw: prior preprocessing may leak across subjects; no strictly inductive claim."
    warning(diagnostics, call. = FALSE)
  }
  k <- max(folds$outer)
  models <- developments <- preprocess <- convergence <- selected <- vector("list", k)
  names(models) <- names(developments) <- names(preprocess) <-
    names(convergence) <- names(selected) <- as.character(seq_len(k))
  predictions <- baseline <- matrix(NA_real_, nrow(Ym), ncol(Ym), dimnames = dimnames(Ym))
  fold_rows <- list()
  pred_rows <- list()
  errors <- list()
  for (f in seq_len(k)) {
    if (verbose)
      message("OUTER fold ", f, "/", k, ": development uses OUTER TRAIN only.")
    train <- which(folds$outer != f)
    outer_test <- which(folds$outer == f)
    development <- .nested_develop(
      .nested_rows(Xm, train), Ym[train, , drop = FALSE], classes[train],
      encoding, features, times, folds$inner[[as.character(f)]],
      configurations, imputation, selection_metric
    )
    model <- development$model
    developments[[f]] <- development[c("inner_results", "tuning", "error")]
    if (imputation_method == "tucker")
      developments[[f]]$imputation_diagnostics <- development$imputation_diagnostics
    if (is.null(model))
      errors[[as.character(f)]] <- development$error else {
      # OUTER TEST predictors are first passed to a fitted transform here.
      pr <- tryCatch(
        .nested_predict_pipeline(model, .nested_rows(Xm, outer_test)),
        error = function(e) e
      )
      if (inherits(pr, "error"))
        errors[[as.character(f)]] <- conditionMessage(pr) else
        predictions[outer_test, ] <- pr
      selected[[f]] <- model$selected_features
      preprocess[[f]] <- model$preprocessing
      convergence[[f]] <- model$convergence
      if (keep_models)
        models[[f]] <- model
    }
    baseline[outer_test, ] <- matrix(
      rep(colMeans(Ym[train, , drop = FALSE]), each = length(outer_test)),
      length(outer_test), ncol(Ym)
    )
    # Only scoring below accesses OUTER TEST labels, after prediction.
    met <- .nested_metrics(Ym[outer_test, , drop = FALSE], predictions[outer_test, , drop = FALSE],
      baseline[outer_test, , drop = FALSE], classes[outer_test], encoding)
    config <- if (is.null(model))
      list(ncomp = NA_integer_, strategy = NA_character_, thr = NA_real_) else model$config
    candidate_diagnostics <- development$inner_results$diagnostic[
      development$inner_results$status != "ok"
    ]
    fold_rows[[f]] <- data.frame(
      fold = f,
      n_train = length(train),
      n_test = length(outer_test),
      selected_ncomp = config$ncomp,
      strategy = config$strategy,
      vip_thr = config$thr,
      n_selected_features = length(selected[[f]]),
      Q2Y = met$Q2Y,
      PRESS = sum(met$PRESS),
      denominator = sum(met$denominator),
      as.list(met$classification),
      status = if (is.null(errors[[as.character(f)]])) "ok" else "failed",
      diagnostic = paste(
        unique(c(errors[[as.character(f)]], candidate_diagnostics)),
        collapse = " | "
      )
    )
    assigned <- if (any(!is.finite(predictions[outer_test, , drop = FALSE])))
      rep(NA_character_, length(outer_test)) else
      .nested_classify(predictions[outer_test, , drop = FALSE], encoding)
    rows <- data.frame(
      subject_id = rep(ids[outer_test], ncol(Ym)),
      outer_fold = f,
      response = rep(colnames(Ym), each = length(outer_test)),
      Y_true = as.vector(Ym[outer_test, , drop = FALSE]),
      Y_pred = as.vector(predictions[outer_test, , drop = FALSE]),
      training_mean = as.vector(baseline[outer_test, , drop = FALSE]),
      true_class = rep(classes[outer_test], ncol(Ym)),
      predicted_class = rep(assigned, ncol(Ym)),
      selected_ncomp = config$ncomp,
      feature_strategy = config$strategy,
      vip_thr = config$thr,
      n_selected_features = length(selected[[f]])
    )
    rows$selected_features <- rep(list(selected[[f]]), nrow(rows))
    pred_rows[[f]] <- rows
  }
  metrics <- .nested_metrics(Ym, predictions, baseline, classes, encoding)
  freq <- data.frame(
    feature = features,
    count = vapply(
      features,
      function(g) sum(vapply(selected, function(s) g %in% s, logical(1))),
      integer(1)
    )
  )
  n_selection_folds <- sum(lengths(selected) > 0L)
  freq$frequency <- if (n_selection_folds)
    freq$count / n_selection_folds else NA_real_
  jaccard <- matrix(NA_real_, k, k, dimnames = list(names(selected), names(selected)))
  for (i in seq_len(k))
    for (j in seq_len(k)) {
      all_genes <- union(selected[[i]], selected[[j]])
      if (length(selected[[i]]) && length(selected[[j]]))
        jaccard[i, j] <- length(intersect(selected[[i]], selected[[j]])) / length(all_genes)
    }
  candidate_diagnostics <- do.call(rbind, lapply(seq_len(k), function(f) {
    z <- developments[[f]]$inner_results
    z <- z[z$status != "ok", , drop = FALSE]
    z$outer_fold <- rep(f, nrow(z))
    z
  }))
  diagnostics <- unique(c(
    diagnostics,
    candidate_diagnostics$diagnostic[candidate_diagnostics$status == "warning"]
  ))
  result <- structure(
    list(
      outer_predictions = do.call(rbind, pred_rows),
      outer_fold_assignments = folds$outer,
      inner_fold_assignments = folds$inner,
      folds = folds,
      inner_tuning_results = developments,
      outer_fold_metrics = do.call(rbind, fold_rows),
      selected_ncomp_per_outer_fold = stats::setNames(
        vapply(fold_rows, function(z) z$selected_ncomp, numeric(1)),
        names(selected)
      ),
      selected_features_per_outer_fold = selected,
      number_of_selected_features_per_fold = lengths(selected),
      preprocessing = preprocess,
      convergence = convergence,
      models = models,
      Q2Y = metrics$Q2Y,
      Q2Y_per_response = metrics$Q2Y_per_response,
      PRESS = metrics$PRESS,
      Q2_denominator = metrics$denominator,
      classification_metrics = metrics$classification,
      errors = errors,
      warnings = diagnostics,
      candidate_diagnostics = candidate_diagnostics,
      stability = list(
        frequency = freq,
        jaccard = jaccard,
        evaluated_folds = names(selected)[lengths(selected) > 0L]
      ),
      settings = list(
        seed = seed,
        outer_folds = k,
        inner_folds = vapply(folds$inner, max, numeric(1)),
        strategy_mode = strategy_mode,
        configurations = configurations,
        imputation = imputation_method,
        input_provenance = input_provenance,
        selection_metric = selection_metric,
        class_levels = encoding$levels,
        response_coding = if (encoding$binary) "binary 0/1" else "one-hot",
        vip_definition = "global unfolded cumulative VIP; model1 all-time mean; model2 cumulative-component sum; M2_mean cumulative-component mean",
        inductive_conditional_on_raw_input = input_provenance == "raw",
        Q2_definition = "1 - pooled outer PRESS / pooled squared error against fold-specific OUTER TRAIN means, per response then averaged",
        scope = "One subject-level nested-CV repetition; no final full-data model is fitted or scored."
      )
    ),
    class = "nested_tensorpls"
  )
  if (imputation_method == "tucker")
    result$settings$tucker_control <- tucker_control
  result
}
