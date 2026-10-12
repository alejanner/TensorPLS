.rank_reference <- function(name) {
  readRDS(test_path("fixtures", "imputation-api", paste0(name, ".rds")))
}
test_that("explicit candidate tuning preserves every pre-existing output exactly", {
  X <- .rank_reference("X_original")
  controls <- .rank_reference("controls")$controls
  candidates <- .rank_reference("custom_candidates")
  old <- .rank_reference("tuning_before")
  current <- tune_tucker_ranks(X, rank_candidates = candidates, repeats = 2L,
    tucker_control = controls,
    verbose = FALSE)
  for (field in names(old))
    expect_identical(current[[field]], old[[field]], info = field)
  expect_identical(class(current), class(old))
  expect_identical(current$results, current$tuning)
  alias <- tune_tucker_ranks(X, candidates = candidates, repeats = 2L, tucker_control = controls,
    verbose = FALSE)
  expect_identical(alias, current)
  expect_identical(as.matrix(current$tuning[, 1:3]), matrix(as.integer(candidates),
    nrow(candidates), 3, dimnames = list(NULL, names(current$tuning)[1:3])))
  good <- which(current$tuning$valid)
  expect_gt(length(good), 0)
  winner <- good[which.min(current$tuning$mean_RMSE[good])]
  expect_identical(unname(current$selected_ranks), as.integer(current$tuning[winner,
    1:3]))
  expect_type(current$selected_ranks, "integer")
  expect_length(current$selected_ranks, 3L)
})
test_that("default ranks reproduce the established grid and ordering", {
  expected <- .rank_reference("default_grid_from_workflow")
  expect_identical(.tucker_rank_grid(c(12, 8, 4)), expected)
  expect_identical(nrow(expected), 37L)
  X <- array(2, c(4, 4, 4))
  tuning <- tune_tucker_ranks(X, verbose = FALSE)
  expect_identical(tuning$tuning[, 1:3], expected)
  expect_identical(tuning$results, tuning$tuning)
  expect_true(all(tuning$tuning$valid))
  expect_identical(unname(tuning$selected_ranks), c(1L, 1L, 1L))
  expect_identical(tuning$settings$repeats, 3L)
  expect_identical(nrow(tuning$repetition_results), 111L)
  expect_identical(tuning, tune_tucker_ranks(X, candidates = NULL, verbose = FALSE))
})
test_that("small dimensions filter grid rows rather than clipping them", {
  X <- array(2, c(3, 2, 2))
  full <- .rank_reference("default_grid_from_workflow")
  expected <- full[apply(full, 1, function(r)
    all(r <= dim(X))), , drop = FALSE]
  rownames(expected) <- NULL
  tuning <- tune_tucker_ranks(X, repeats = 1L, verbose = FALSE)
  expect_identical(tuning$tuning[, 1:3], expected)
  expect_identical(.tucker_rank_grid(dim(X)), expected)
  expect_false(anyDuplicated(tuning$tuning[, 1:3]) > 0)
  for (i in seq_len(nrow(expected))) {
    r <- as.integer(expected[i, ])
    expect_true(all(r <= dim(X)))
    expect_true(all(vapply(1:3, function(j)
      r[j] <= prod(r[-j]), logical(1))))
  }
  expect_error(tune_tucker_ranks(X, candidates = rbind(c(1, 1, 1),
    c(4, 4, 3))), "dimensions")
  expect_error(tune_tucker_ranks(X, candidates = rbind(c(1, 1, 1),
    c(3, 1, 1))), "product")
  expect_error(tune_tucker_ranks(X, candidates = matrix(c(1.5, 1, 1),
    1)), "integers")
  expect_error(tune_tucker_ranks(X, candidates = matrix(c(1, 1, 1),
    1), rank_candidates = matrix(c(1, 1, 1), 1)), "only one")
})
