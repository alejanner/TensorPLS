# Train/apply Tucker imputation. Tensor mode products below are algebraic
# contractions, not a flattened matrix imputation model.
.tk_validate <- function(X) {
  if (!is.numeric(X) || length(dim(X)) != 3L || any(dim(X) < 1L) || any(is.infinite(X)))
    stop('X must be a numeric subjects x features x time tensor without infinite values.')
}
.tk_mode <- function(X, M, mode) {
  d <- dim(X)
  order <- c(mode, setdiff(1:3, mode))
  z <- M %*% matrix(aperm(X, order), d[mode])
  aperm(array(z, c(nrow(M), d[-mode])), order(order))
}
.tk_core <- function(X, A, B, C)
  .tk_mode(.tk_mode(.tk_mode(X, t(A), 1), t(B), 2),
  t(C), 3)
.tk_reconstruct <- function(G, A, B, C)
  .tk_mode(.tk_mode(.tk_mode(G, A, 1), B,
  2), C, 3)
# HOSVD start from the training-only, zero-filled centered tensor.
.tk_initial <- function(X, mode, rank) {
  m <- matrix(aperm(X, c(mode, setdiff(1:3, mode))), dim(X)[mode])
  if (nrow(m) <= ncol(m))
    return(eigen(tcrossprod(m), symmetric = TRUE)$vectors[, seq_len(rank), drop = FALSE])
  e <- eigen(crossprod(m), symmetric = TRUE)
  good <- which(e$values[seq_len(rank)] > max(e$values) * 1e-12 & e$values[seq_len(rank)] >
    0)
  u <- if (length(good))
    sweep(m %*% e$vectors[, good, drop = FALSE], 2, sqrt(e$values[good]),
      "/") else
        matrix(numeric(), nrow(m), 0)
  if (ncol(u) < rank) {
    extra <- matrix(0, nrow(m), rank)
    extra[cbind(seq_len(rank), seq_len(rank))] <- 1
    u <- cbind(u, extra)
  }
  qr.Q(qr(u))[, seq_len(rank), drop = FALSE]
}
# Exact observed-entry least-squares coordinate update with frozen G/B/C.
.tk_solve <- function(D, z, tolerance) {
  s <- svd(D, nu = min(dim(D)), nv = min(dim(D)))
  keep <- which(s$d > max(s$d) * tolerance & s$d > 0)
  if (!length(keep))
    return(matrix(0, ncol(D), ncol(z)))
  s$v[, keep, drop = FALSE] %*% sweep(crossprod(s$u[, keep, drop = FALSE], z), 1,
    s$d[keep], "/")
}
.tk_observed_coordinates <- function(Z, basis, patterns, groups, tolerance, penalty = NULL) {
  a <- matrix(0, dim(Z)[1], ncol(basis))
  for (g in seq_along(patterns)) {
    ids <- which(groups == g)
    obs <- patterns[[g]]
    z <- vapply(ids, function(i)
      as.vector(Z[i, , , drop = FALSE])[obs], numeric(length(obs)))
    z <- matrix(z, length(obs), length(ids))
    D <- basis[obs, , drop = FALSE]
    if (!is.null(penalty)) {
      D <- rbind(D, penalty)
      z <- rbind(z, matrix(0, nrow(penalty), length(ids)))
    }
    a[ids, ] <- t(.tk_solve(D, z, tolerance))
  }
  a
}
# Group identical observation masks for exact batched least-squares updates.
.tk_patterns <- function(observed, mode) {
  z <- aperm(observed, c(mode, setdiff(1:3, mode)))
  patterns <- list()
  groups <- integer(dim(z)[1])
  for (i in seq_len(dim(z)[1])) {
    obs <- which(as.vector(z[i, , , drop = FALSE]))
    g <- which(vapply(patterns, identical, logical(1), obs))
    if (!length(g)) {
      patterns[[length(patterns) + 1L]] <- obs
      g <- length(patterns)
    }
    groups[i] <- g[1]
  }
  list(patterns = patterns, groups = groups)
}
# Equivalent QR compression for time-factor updates (avoid a huge dense SVD).
.tk_time_coordinates <- function(Z, observed, factors, G, ridge, tolerance) {
  ranks <- dim(G)
  out <- matrix(0, dim(Z)[3], ranks[3])
  masks <- qrs <- list()
  for (ti in seq_len(dim(Z)[3])) {
    blocks <- targets <- list()
    row <- 0L
    for (i in seq_len(dim(Z)[1])) {
      obs <- which(observed[i, , ti])
      if (!length(obs))
        next
      group <- which(vapply(masks, identical, logical(1), obs))
      if (!length(group)) {
        masks[[length(masks) + 1L]] <- obs
        group <- length(masks)
        q <- qr(factors[[2]][obs, , drop = FALSE])
        k <- min(length(obs), ranks[2])
        qrs[[group]] <- list(Q = qr.Q(q)[, seq_len(k), drop = FALSE], R = qr.R(q)[seq_len(k),
          order(q$pivot), drop = FALSE])
      }
      q <- qrs[[group[1]]]
      row <- row + 1L
      h <- matrix(.tk_mode(G, factors[[1]][i, , drop = FALSE], 1),
        ranks[2], ranks[3])
      blocks[[row]] <- q$R %*% h
      targets[[row]] <- crossprod(q$Q, Z[i, obs, ti])
    }
    D <- do.call(rbind, blocks)
    target <- matrix(unlist(targets), ncol = 1)
    if (ridge > 0) {
      penalty <- sqrt(ridge) * t(matrix(aperm(G, c(3, 1, 2)),
        ranks[3]))
      D <- rbind(D, penalty)
      target <- rbind(target, matrix(0, nrow(penalty), 1))
    }
    out[ti, ] <- as.vector(.tk_solve(D, target, tolerance))
  }
  out
}
# Exact weighted core fit with a small QR-compressed design. For each i,t,
# B_obs=Q R. ||z-B_obs H a||^2 splits into ||Q' z-R H a||^2 plus a
# constant orthogonal residual. Thus stacking these small blocks preserves
# the full observed-entry least-squares objective without normal equations.
.tk_weighted_core <- function(Z, observed, factors, ranks, ridge = 0) {
  blocks <- targets <- list()
  masks <- qrs <- list()
  row <- 0L
  for (ti in seq_len(dim(Z)[3]))
    for (i in seq_len(dim(Z)[1])) {
    obs <- which(observed[i, , ti])
    if (!length(obs))
      next
    group <- which(vapply(masks, identical, logical(1), obs))
    if (!length(group)) {
      masks[[length(masks) + 1L]] <- obs
      group <- length(masks)
      q <- qr(factors[[2]][obs, , drop = FALSE])
      k <- min(length(obs), ranks[2])
      qrs[[group]] <- list(Q = qr.Q(q)[, seq_len(k), drop = FALSE], R = qr.R(q)[seq_len(k),
        order(q$pivot), drop = FALSE])
    }
    q <- qrs[[group[1]]]
    row <- row + 1L
    blocks[[row]] <- kronecker(factors[[3]][ti, , drop = FALSE], kronecker(q$R,
      factors[[1]][i, , drop = FALSE]))
    targets[[row]] <- crossprod(q$Q, Z[i, obs, ti])
  }
  D <- do.call(rbind, blocks)
  target <- matrix(unlist(targets), ncol = 1)
  if (ridge > 0) {
    D <- rbind(D, sqrt(ridge) * diag(prod(ranks)))
    target <- rbind(target, matrix(0, prod(ranks), 1))
  }
  array(.tk_solve(D, target, 1e-12), ranks)
}
#' Fit a training-only three-way Tucker imputer
#' @param X_train Numeric subjects x features x time tensor; NA denotes missing.
#' @param ranks Fixed multilinear ranks. No rank selection is performed.
#' @param tolerance Relative reconstruction-change convergence tolerance.
#' @param max_iter Maximum observed-entry alternating least-squares iterations.
#' @param seed Seed for random orthonormal starts; caller RNG state is preserved.
#' @param ridge Nonnegative fixed penalty on the squared Frobenius norm of the
#' reconstructed training tensor. Default 0.001; zero gives unregularized fitting.
#' This parameter is not selected using held-out data.
#' @param initialization HOSVD (default, deterministic) or seeded random starts.
#' @param projection_tolerance Relative singular-value cutoff for unseen projection.
#' @return A tucker_imputer containing G, B, C, training coordinates A_train,
#' training-only centering, ranks, completed training tensor, and diagnostics.
#' @details The training tensor is centered using each feature/time cell's
#' observed training mean. Factors and core minimize observed-entry least squares
#' by alternating SVD solves. This is a nonconvex local fit; reconstruction
#' replaces missing values only after fitting.
#' Non-convergence emits a warning; prediction refuses a non-converged fit.
#' @export
fit_tucker_imputer <- function(X_train, ranks = c(4L, 4L, 3L), tolerance = 1e-6,
  max_iter = 1000L, seed = 123L, projection_tolerance = 1e-10, initialization = c("hosvd",
    "random"), ridge = 1e-3) {
  initialization <- match.arg(initialization)
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge <
    0)
    stop("ridge must be one finite nonnegative value.")
  .tk_validate(X_train)
  d <- dim(X_train)
  if (!is.numeric(ranks) || length(ranks) != 3L || any(!is.finite(ranks)) || any(ranks <
    1 | ranks != as.integer(ranks)) || any(ranks > d))
    stop('ranks must be three positive integers no larger than training dimensions.')
  if (any(vapply(1:3, function(j)
    ranks[j] > prod(ranks[-j]), logical(1))))
    stop('Each Tucker rank must be no larger than the product of the other two.')
  .nested_integer(max_iter, 'max_iter')
  .nested_integer(seed, 'seed', 0L)
  for (v in list(tolerance, projection_tolerance))
    if (!is.numeric(v) || length(v) !=
    1L || !is.finite(v) || v <= 0 || v >= 1)
    stop('Tolerances must be finite scalars strictly between zero and one.')
  observed <- !is.na(X_train)
  if (!any(observed))
    stop('Training tensor has no observed entries.')
  if (any(apply(observed, 1, sum) == 0))
    stop('A training subject has no observed entries.')
  if (any(apply(observed, 2, sum) == 0) || any(apply(observed, 3,
    sum) == 0))
    stop('A training feature or time has no observed entries; its Tucker factor is unidentified.')
  counts <- apply(observed, c(2, 3), sum)
  center <- apply(X_train, c(2, 3), mean, na.rm = TRUE)
  absent <- !is.finite(center)
  # A feature-time cell absent in training uses a documented training-feature mean.
  feature_mean <- apply(X_train, 2, mean, na.rm = TRUE)
  center[absent] <- matrix(rep(feature_mean, d[3]), d[2], d[3])[absent]
  Z <- sweep(X_train, c(2, 3), center, '-')
  Z[!observed] <- 0
  scale_ref <- max(sqrt(sum(Z[observed]^2)), sqrt(.Machine$double.eps))
  factors <- if (initialization == "hosvd")
    lapply(1:3, function(j)
      .tk_initial(Z, j, ranks[j])) else
        .vip_with_seed(seed, function()
      lapply(1:3, function(j)
      qr.Q(qr(matrix(stats::rnorm(d[j] *
    ranks[j]), d[j], ranks[j])))[, seq_len(ranks[j]), drop = FALSE]))
  mode_patterns <- lapply(1:3, function(j)
    .tk_patterns(observed, j))
  # Missing values remain zero in this training objective; only observed
  # entries enter any least-squares factor or core update.
  G <- .tk_core(Z, factors[[1]], factors[[2]], factors[[3]])
  previous <- NULL
  history <- vector('list', max_iter)
  delta <- Inf
  for (it in seq_len(max_iter)) {
    for (j in 1:3) {
      other <- setdiff(1:3, j)
      Gu <- matrix(aperm(G, c(j, other)), ranks[j])
      if (j == 3L)
        a <- .tk_time_coordinates(Z, observed, factors, G, ridge, projection_tolerance) else {
        D <- kronecker(factors[[other[2]]], factors[[other[1]]]) %*%
          t(Gu)
        permuted <- aperm(Z, c(j, other))
        pattern <- mode_patterns[[j]]
        a <- .tk_observed_coordinates(permuted, D, pattern$patterns, pattern$groups,
          projection_tolerance, if (ridge > 0)
          sqrt(ridge) * t(Gu) else
            NULL)
      }
      qa <- qr(a)
      factors[[j]] <- qr.Q(qa)[, seq_len(ranks[j]), drop = FALSE]
      G <- .tk_mode(G, qr.R(qa)[seq_len(ranks[j]), order(qa$pivot), drop = FALSE],
        j)
    }
    G <- .tk_weighted_core(Z, observed, factors, ranks, ridge)
    reconstructed <- .tk_reconstruct(G, factors[[1]], factors[[2]], factors[[3]])
    delta <- if (is.null(previous))
      Inf else
        sqrt(sum((reconstructed - previous)^2)) / max(scale_ref, sqrt(sum(previous^2)))
    residual <- Z[observed] - reconstructed[observed]
    history[[it]] <- data.frame(iteration = it, delta = delta,
      observed_RMSE = sqrt(mean(residual^2)),
      objective = sum(residual^2) + ridge * sum(G^2))
    if (is.finite(delta) && delta < tolerance)
      break
    previous <- reconstructed
  }
  converged <- is.finite(delta) && delta < tolerance
  if (!converged)
    warning('Tucker imputer did not converge: iterations=', it, ', delta=',
      signif(delta, 6), ', tolerance=', tolerance, call. = FALSE)
  complete <- X_train
  reconstructed_raw <- sweep(reconstructed, c(2, 3), center, '+')
  complete[!observed] <- reconstructed_raw[!observed]
  for (j in 1:3)
    rownames(factors[[j]]) <- dimnames(X_train)[[j]]
  # R ordering: feature index varies fastest within each subject's feature/time slice.
  # vec(X_i - center) = (C %x% B) t(G_(1)) a_i = D a_i.
  basis <- kronecker(factors[[3]], factors[[2]]) %*% t(matrix(G, ranks[1]))
  structure(list(G = G, B = factors[[2]], C = factors[[3]], A_train = factors[[1]],
    basis = basis, preprocessing = list(center = center, counts = counts, absent_cells = absent,
      scaling = 'none'), ranks = as.integer(ranks), feature_names = dimnames(X_train)[[2]],
    time_names = dimnames(X_train)[[3]], dimensions = d, training_ids = dimnames(X_train)[[1]],
    X_train_imputed = complete, projection_tolerance = projection_tolerance, seed = seed,
    initialization = initialization, ridge = ridge, convergence = list(converged = converged,
      iterations = it, final_delta = delta, tolerance = tolerance, ranks = as.integer(ranks),
      method = 'observed-entry weighted Tucker ALS', ridge = ridge,
        objective = sum(residual^2) + ridge *
        sum(G^2), observed_RMSE = sqrt(mean(residual^2)), observed_SSE = sum(residual^2),
      n_unobserved_training_cells = sum(absent)), history = do.call(rbind, history[seq_len(it)])),
    class = 'tucker_imputer')
}
#' Project unseen subjects into a frozen Tucker imputer
#' @param fit A converged fit_tucker_imputer result.
#' @param X_test Three-way tensor with identical ordered feature/time axes.
#' @param diagnostics If TRUE return X, coordinates and per-subject diagnostics.
#' @return A three-way completed tensor, or a diagnostics list. Observed entries
#' are copied unchanged. No factor, core or preprocessing parameter is refitted.
#' @details For each subject independently, solve D_obs a = x_obs - center_obs
#' by truncated SVD, D=(C kronecker B) t(G_(1)). Singular values below the relative
#' cutoff use a minimum-norm least-squares solution and report rank deficiency.
#' A wholly unobserved subject has no identifiable coordinates and is rejected.
#' @export
predict_tucker_imputer <- function(fit, X_test, diagnostics = FALSE) {
  .tk_validate(X_test)
  .vip_flag(diagnostics, 'diagnostics')
  if (!inherits(fit, 'tucker_imputer'))
    stop('fit must be a fitted Tucker imputer.')
  if (!isTRUE(fit$convergence$converged))
    stop('Cannot apply a non-converged Tucker imputer.')
  if (!identical(dim(X_test)[2:3], fit$dimensions[2:3]) || !identical(dimnames(X_test)[[2]],
    fit$feature_names) || !identical(dimnames(X_test)[[3]], fit$time_names))
    stop('New tensor must have identical ordered feature/time dimensions and names.')
  completed <- X_test
  n <- dim(X_test)[1]
  rows <- vector('list', n)
  coordinates <- matrix(NA_real_, n, fit$ranks[1], dimnames = list(dimnames(X_test)[[1]],
    NULL))
  for (i in seq_len(n)) {
    x <- as.vector(X_test[i, , , drop = FALSE])
    obs <- which(!is.na(x))
    if (!length(obs))
      stop('Unseen subject ', i, ' has no observed entries; projection is unidentified.')
    D <- fit$basis[obs, , drop = FALSE]
    z <- x[obs] - as.vector(fit$preprocessing$center)[obs]
    s <- svd(D, nu = min(dim(D)), nv = min(dim(D)))
    keep <- which(s$d > max(s$d) * fit$projection_tolerance & s$d > 0)
    a <- if (length(keep))
      as.vector(s$v[, keep, drop = FALSE] %*% (crossprod(s$u[, keep, drop = FALSE],
        z) / s$d[keep])) else
          rep(0, fit$ranks[1])
    coordinates[i, ] <- a
    reconstructed <- as.vector(fit$basis %*% a) + as.vector(fit$preprocessing$center)
    missing <- is.na(x)
    x[missing] <- reconstructed[missing]
    completed[i, , ] <- matrix(x, dim(X_test)[2], dim(X_test)[3])
    rows[[i]] <- data.frame(subject = if (is.null(dimnames(X_test)[[1]]))
      as.character(i) else
        dimnames(X_test)[[1]][i], n_observed = length(obs),
        n_missing = sum(missing),
      effective_rank = length(keep), rank_deficient = length(keep) < fit$ranks[1],
      observed_RMSE = sqrt(mean((z - as.vector(D %*% a))^2)))
  }
  report <- do.call(rbind, rows)
  if (any(report$rank_deficient))
    warning(sum(report$rank_deficient),
      ' subject projections are rank deficient; used minimum-norm observed-entry least squares.', call. = FALSE)
  if (diagnostics)
    list(X = completed, coordinates = coordinates, diagnostics = report) else
      completed
}
