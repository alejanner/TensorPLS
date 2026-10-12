# Training-only NIPALS core shared by ordinary and nested CV. No held-out Y input.
.tensorpls_fit <- function(X, Y, ncomp, tolerance = 1e-6, max_iter = 100L) {
  if (!is.matrix(X) || !is.matrix(Y) || !is.numeric(X) || !is.numeric(Y) || any(!is.finite(X)) ||
    any(!is.finite(Y)))
    stop("PLS needs finite numeric matrices.")
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  if (nrow(Y) != n || ncomp > min(n - 1L, p) || ncomp < 1L)
    stop("Insufficient training dimensions for requested ncomp.")
  xm <- colMeans(X)
  ym <- colMeans(Y)
  xs <- apply(X, 2, stats::sd)
  ys <- apply(Y, 2, stats::sd)
  xs[xs == 0] <- 1
  ys[ys == 0] <- 1
  E <- sweep(sweep(X, 2, xm, "-"), 2, xs, "/")
  F <- sweep(sweep(Y, 2, ym, "-"), 2, ys, "/")
  W <- P <- matrix(0, p, ncomp)
  Tt <- U <- matrix(0, n, ncomp)
  C <- matrix(0, q, ncomp)
  conv <- data.frame(component = seq_len(ncomp), converged = FALSE, iterations = 0L,
    final_delta = NA_real_, tolerance = tolerance)
  for (h in seq_len(ncomp)) {
    u <- F[, 1L]
    if (sum(u^2) < 1e-24)
      u <- F[, which.max(colSums(F^2))]
    w <- rep(1, p)
    for (it in seq_len(max_iter)) {
      wn <- as.vector(crossprod(E, u))
      wn_norm <- sqrt(sum(wn^2))
      if (!is.finite(wn_norm) || wn_norm < 1e-12)
        stop("PLS exhausted numerical rank at component ", h)
      wn <- wn / wn_norm
      tn <- as.vector(E %*% wn)
      tt <- sum(tn^2)
      if (!is.finite(tt) || tt < 1e-12)
        stop("PLS score has zero numerical variance.")
      cn <- as.vector(crossprod(F, tn)) / tt
      cc <- sum(cn^2)
      if (!is.finite(cc) || cc < 1e-24)
        stop("PLS response covariance exhausted.")
      un <- as.vector(F %*% cn) / cc
      delta <- sum((wn - w)^2)
      if (delta < tolerance)
        break
      w <- wn
      u <- un
    }
    conv[h, c("converged", "iterations", "final_delta")] <- list(delta < tolerance, it, delta)
    pn <- as.vector(crossprod(E, tn)) / tt
    W[, h] <- wn
    P[, h] <- pn
    Tt[, h] <- tn
    U[, h] <- un
    C[, h] <- cn
    E <- E - tcrossprod(tn, pn)
    F <- F - tcrossprod(tn, cn)
  }
  # VIP uses squared response-score correlations.
  cor2 <- suppressWarnings(stats::cor(Y, Tt)^2)
  cor2[!is.finite(cor2)] <- 0
  r2 <- colMeans(cor2)
  Rd <- matrix(0, ncomp, ncomp)
  for (j in seq_len(ncomp))
    Rd[seq_len(j), j] <- r2[seq_len(j)]
  denom <- cumsum(r2)
  if (any(denom <= 0))
    stop("VIP has no explained response signal.")
  VIP <- sqrt((W^2) %*% Rd %*% diag(p / denom, ncomp))
  Ws <- W %*% solve(crossprod(P, W), tol = 1e-20)
  Bs <- Ws %*% t(C)
  Br <- sweep(sweep(Bs, 1, xs, "/"), 2, ys, "*")
  intercept <- as.vector(ym - xm %*% Br)
  pred <- sweep(X %*% Br, 2, intercept, "+")
  structure(list(x.scores = Tt, x.loads = P, y.scores = U, y.loads = C, raw.wgs = W,
    mod.wgs = Ws, std.coefs = Bs, reg.coefs = rbind(Br, INTERCEPT = intercept), y.pred = pred,
    resid = Y - pred, VIP = VIP, preprocessing = list(x_mean = xm, x_sd = xs, y_mean = ym,
      y_sd = ys), convergence = conv, ncomp = ncomp, Q2 = NULL, Q2cum = NULL), class = "pls_reg")
}
.tensorpls_predict <- function(object, newdata, ncomp = object$ncomp) {
  if (!is.matrix(newdata) || ncol(newdata) != nrow(object$raw.wgs) || any(!is.finite(newdata)))
    stop("Prediction inputs must match training predictors and be finite.")
  idx <- seq_len(ncomp)
  w <- object$raw.wgs[, idx, drop = FALSE]
  p <- object$x.loads[, idx, drop = FALSE]
  coef <- w %*% solve(crossprod(p, w), tol = 1e-20) %*% t(object$y.loads[, idx, drop = FALSE])
  z <- sweep(sweep(newdata, 2, object$preprocessing$x_mean, "-"), 2,
    object$preprocessing$x_sd, "/")
  sweep(sweep(z %*% coef, 2, object$preprocessing$y_sd, "*"), 2, object$preprocessing$y_mean,
    "+")
}
