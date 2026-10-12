native_fixture <- function()
  readRDS(system.file('extdata', 'vip_example.rds', package = 'TensorPLS'))
test_that('expanded ncomp grid is selected by inner scores and reproduced at seed 123', {
  d <- native_fixture()
  run <- function(Y = d$Y, folds = NULL, grid = 1:8)
    nested_cv_tensorpls(d$X,
    Y, outer_folds = 3, inner_folds = 3, seed = 123, ncomp = grid,
    folds = folds, feature_strategy = 'VIP2D_intersect_VIP3Dmodel2', vip_thr = 50,
      imputation = 'training_mean',
    input_provenance = 'raw', verbose = FALSE)
  a <- run()
  b <- run()
  expect_identical(a, b)
  for (f in 1:3) {
    t <- a$inner_tuning_results[[f]]$tuning
    expect_identical(t$ncomp, 1:8)
    good <- which(is.finite(t$selection_score))
    chosen <- good[order(-t$selection_score[good], t$ncomp[good], good)[1]]
    expect_identical(unname(a$selected_ncomp_per_outer_fold[f]), as.numeric(t$ncomp[chosen]))
  }
  changed <- d$Y
  changed[a$folds$outer == 1, 1] <- 1 - changed[a$folds$outer == 1,
    1]
  c <- run(changed, a$folds)
  expect_identical(a$selected_ncomp_per_outer_fold[1], c$selected_ncomp_per_outer_fold[1])
  expect_identical(a$inner_tuning_results[[1]], c$inner_tuning_results[[1]])
  expect_identical(a$models[[1]], c$models[[1]])
  expect_identical(subset(a$outer_predictions, outer_fold == 1)$Y_pred, subset(c$outer_predictions,
    outer_fold == 1)$Y_pred)
  expect_false(identical(a$outer_fold_metrics$PRESS[1], c$outer_fold_metrics$PRESS[1]))
  ten <- run(grid = 1:10)
  expect_true(all(vapply(ten$inner_tuning_results, function(r)
    identical(r$tuning$ncomp,
    1:10), logical(1))))
  expect_true(all(ten$outer_fold_metrics$status == 'ok'))
})
