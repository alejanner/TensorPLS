# Run from the package root; optionally supply a different output file.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else file.path("inst", "extdata", "vip_example.rds")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261012L)

n <- 24L
p <- 16L
k <- 3L
subjects <- sprintf("SyntheticSubject%02d", seq_len(n))
features <- sprintf("Feature%02d", seq_len(p))
times <- paste0("T", seq_len(k))
classes <- rep(c(0, 1), length.out = n)

subject_effect <- matrix(rnorm(n * p, sd = 0.5), n, p)
feature_effect <- seq(-0.75, 0.75, length.out = p)
time_effect <- c(0.8, 1, 1.2)
X <- array(rnorm(n * p * k), c(n, p, k),
  dimnames = list(Subject = subjects, Feature = features, Time = times))
for (t in seq_len(k)) {
  X[, , t] <- X[, , t] + subject_effect +
    outer(2 * classes - 1, feature_effect) * time_effect[t]
}
Y <- matrix(classes, ncol = 1L, dimnames = list(subjects, "Case"))
fixture <- list(
  X = X,
  Y = Y,
  provenance = paste(
    "Fully synthetic Gaussian longitudinal fixture; seed 20261012;",
    "24 synthetic subjects, 16 synthetic features, 3 time points;",
    "balanced constructed binary classes, subject effects and a fixed class-time signal."
  )
)
stopifnot(all(is.finite(X)), sum(Y == 0) == 12L, sum(Y == 1) == 12L)
saveRDS(fixture, output, version = 2)
