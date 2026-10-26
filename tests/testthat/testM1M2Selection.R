.m1_m2_vips <- function(p = 6L, k = 3L, H = 2L) {
  features <- paste0('gene_', seq_len(p))
  times <- paste0('time_', seq_len(k))
  scores <- matrix((seq_len(p * k * H) * 37L %% 97) / 10, p * k, H)
  v <- .nested_vips(list(VIP = scores, ncomp = H), features, times)
  components <- paste0('Comp', seq_len(H))
  v$V <- array(v$VIP2D, c(p, k, H), list(features, times, components))
  v$feature_names <- features
  v$time_names <- times
  v$component_names <- components
  v$ncomp <- H
  class(v) <- 'tensorpls_vips'
  v
}

test_that('retained selections match the frozen original outputs', {
  ref <- readRDS(test_path('fixtures', 'm1_m2_selection_reference.rds'))
  for (s in c('Original', 'VIP2D', 'VIP3Dmodel2', 'VIP2D_intersect_VIP3Dmodel2')) {
    expect_identical(vip_feature_sets(ref$vips, strategy = s, thr = 35), ref$strategies[[s]])
  }
  m1 <- vip_feature_sets(ref$vips, strategy = 'M1', thr = 35)
  expect_identical(m1$selected_features, ref$strategies$VIP3Dmodel1$selected_features)
  expect_identical(m1$thresholds$M1, ref$strategies$VIP3Dmodel1$thresholds$VIP3Dmodel1)
  expect_identical(ref$vips$VIP3Dmodel2, apply(ref$vips$V, c(1, 2), sum))
})

test_that('M1 and M2 resolve backing matrices without changing public names', {
  for (dims in list(c(6L, 5L, 1L), c(6L, 2L, 3L), c(1L, 1L, 1L), c(6L, 1L, 2L))) {
    v <- .m1_m2_vips(dims[1], dims[2], dims[3])
    before <- serialize(v, NULL)
    for (thr in c(0, 50, 99, 100)) {
      m1 <- vip_feature_sets(v, strategy = 'M1', thr = thr)
      expected_m1 <- feature_selection(v$VIP3Dmodel1, thr = thr, strip_time = FALSE)
      expect_identical(m1$selected_features, v$feature_names[v$feature_names %in% expected_m1])
      expect_identical(m1$selected_indices, match(expected_m1, v$feature_names))
      expect_identical(m1$thresholds$M1, apply(v$VIP3Dmodel1, 2, quantile, probs = thr / 100, names = FALSE))
      sources <- vip_feature_sets(v, sources = c('VIP2D', 'M1', 'M2'), thr = thr)
      for (s in c('M1', 'M2', 'VIP2D_intersect_M1', 'VIP2D_intersect_M2',
                  'M1_intersect_M2', 'VIP2D_intersect_M1_intersect_M2')) {
        result <- vip_feature_sets(v, strategy = s, thr = thr)
        expect_identical(result$strategy, s)
        expect_identical(result$selected_features,
          v$feature_names[v$feature_names %in% sources$sets[[s]]])
      }
      m2 <- vip_feature_sets(v, strategy = 'M2', thr = thr)
      expected <- feature_selection(v$M2_mean, thr = thr, strip_time = FALSE)
      expect_identical(m2$selected_features, v$feature_names[v$feature_names %in% expected])
      expect_identical(m2$thresholds$M2, apply(v$M2_mean, 2, quantile,
        probs = thr / 100, type = 7, names = FALSE))

    }
    expect_identical(serialize(v, NULL), before)
  }
})

test_that('M2 selection uses mean values even when sum selects different features', {
  v <- .m1_m2_vips()
  v$M2_mean[,] <- 0
  v$M2_mean[2, ] <- 5
  v$VIP3Dmodel2[,] <- 0
  v$VIP3Dmodel2[5, ] <- 9
  v$VIP2D[,] <- 1
  mean_result <- vip_feature_sets(v, strategy = 'VIP2D_intersect_M2', thr = 99)
  sum_result <- vip_feature_sets(v, strategy = 'VIP2D_intersect_VIP3Dmodel2', thr = 99)
  expect_identical(mean_result$selected_features, 'gene_2')
  expect_identical(sum_result$selected_features, 'gene_5')
  expect_identical(mean_result$thresholds$M2,
    apply(v$M2_mean, 2, quantile, probs = .99, names = FALSE))
  fs <- vip_feature_sets(v, sources = c('VIP2D', 'M1', 'M2'), thr = 99)
  expect_identical(fs$settings$sources, c('VIP2D', 'M1', 'M2'))
  expect_identical(fs$selection$source, c('VIP2D', 'M1', 'M2'))
  expect_identical(fs$summary$sources[[5]], c('VIP2D', 'M2'))
})

test_that('ties, column aggregation, intersections and validation remain explicit', {
  v <- .m1_m2_vips(6L, 2L, 2L)
  v$M2_mean <- matrix(c(1, 1, 2, 2, 3, 3, 3, 3, 2, 2, 1, 1), 6, 2,
    dimnames = dimnames(v$M2_mean))
  for (aggregator in c('any', 'all', 'mean')) {
    z <- vip_feature_sets(v, sources = 'M2', thr = 50, aggregator = aggregator)
    expect_identical(z$sets$M2, feature_selection(v$M2_mean, thr = 50,
      aggregator = aggregator, strip_time = FALSE))
  }
  z <- vip_feature_sets(v, strategy = 'VIP2D_intersect_M2', thr = 50)
  unfolded <- v$VIP2D
  rownames(unfolded) <- rep(v$feature_names, length(v$time_names))
  expected <- intersect(feature_selection(unfolded, thr = 50, strip_time = FALSE),
    feature_selection(v$M2_mean, thr = 50, strip_time = FALSE))
  expect_identical(z$selected_features, v$feature_names[v$feature_names %in% expected])
  expect_false(anyDuplicated(z$selected_features) > 0)
  for (sources in list(c('M1', 'VIP3Dmodel1'), c('M2', 'M2_mean'))) {
    expect_error(vip_feature_sets(v, sources = sources), 'sources must')
  }
  for (s in c('M2_intersect_VIP2D', 'M1_intersect_M1', 'VIP2D_intersect_unknown', NA_character_)) {
    expect_error(vip_feature_sets(v, strategy = s), 'Unsupported')
  }
  expect_error(vip_feature_sets(v, sources = 'bad'), 'sources must')
  missing <- v; missing$M2_mean <- NULL
  expect_error(vip_feature_sets(missing, sources = 'M2'), 'missing.*M2_mean')
  bad <- v; colnames(bad$M2_mean)[1] <- 'wrong'
  expect_error(vip_feature_sets(bad, strategy = 'M2'), 'M2_mean.*labels')
  expect_error(vip_feature_sets(v, strategy = 'M2', thr = 101), 'percentile')
  expect_error(vip_feature_sets(v, strategy = 'M2', sources = 'M2'), 'default sources')
})

test_that('quantile boundary fixtures preserve separate sum and mean decisions', {
  for (H in c(9L, 10L)) {
    v <- .m1_m2_vips(200L, 4L, H)
    scores <- round(matrix(1 + ((seq_len(800 * H) * 37L) %% 997) / 997, 800, H), 1)
    backing <- .nested_vips(list(VIP = scores, ncomp = H), v$feature_names, v$time_names)
    for (nm in names(backing)) v[[nm]] <- backing[[nm]]
    v$V[] <- scores
    for (pair in list(c('M2', 'M2_mean'), c('VIP3Dmodel2', 'VIP3Dmodel2'))) {
      z <- vip_feature_sets(v, strategy = pair[1], thr = 99)
      expected <- feature_selection(v[[pair[2]]], thr = 99, strip_time = FALSE)
      expect_identical(z$selected_features, v$feature_names[v$feature_names %in% expected])
    }
    mean_threshold <- apply(v$M2_mean, 2, quantile, probs = .99, names = FALSE)
    sum_threshold <- apply(v$VIP3Dmodel2, 2, quantile, probs = .99, names = FALSE)
    mean_pass <- sweep(v$M2_mean, 2, mean_threshold, '>=')
    sum_pass <- sweep(v$VIP3Dmodel2, 2, sum_threshold, '>=')
    expect_true(any(mean_pass != sum_pass))
  }
})

.m1_m2_data <- function() {
  set.seed(72)
  X <- array(rnorm(24 * 12 * 3), c(24, 12, 3),
    dimnames = list(paste0('s', 1:24), paste0('gene_', 1:12), paste0('time_', 1:3)))
  Y <- matrix(rep(0:1, 12), ncol = 1, dimnames = list(dimnames(X)[[1]], 'Case'))
  X[, 1:3, ] <- X[, 1:3, ] + 2 * Y[, 1]
  X[1, 5, 2] <- NA_real_
  list(X = X, Y = Y, folds = make_tensorpls_folds(Y, 2L, 2L, 123L, dimnames(X)[[1]]))
}

test_that('new strategies dispatch through inner fits and outer all-time refits', {
  d <- .m1_m2_data()
  original_vips <- .nested_vips
  original_pipeline <- .nested_fit_pipeline
  stages <- list()
  local_mocked_bindings(.nested_vips = function(fit, features, times) {
    v <- original_vips(fit, features, times)
    v$M2_mean[,] <- 0
    v$M2_mean[1:6, ] <- 5
    v$VIP3Dmodel2[,] <- 0
    v$VIP3Dmodel2[7:12, ] <- 5
    v
  }, .nested_fit_pipeline = function(X, Y, features, times, config, imputation) {
    z <- original_pipeline(X, Y, features, times, config, imputation)
    used <- strsplit(config$strategy, '_intersect_', fixed = TRUE)[[1]]
    mapping <- c(VIP2D = 'VIP2D', M1 = 'VIP3Dmodel1', M2 = 'M2_mean',
      VIP3Dmodel2 = 'VIP3Dmodel2')
    sets <- lapply(used, function(src) {
      a <- z$vip[[mapping[[src]]]]
      if (src == 'VIP2D') rownames(a) <- rep(features, length(times))
      feature_selection(a, thr = config$thr, strip_time = FALSE)
    })
    expected <- features[features %in% Reduce(intersect, sets)]
    expect_identical(z$selected_features, expected)
    expect_identical(z$columns, which(rep(features, length(times)) %in% expected))
    stages[[length(stages) + 1L]] <<- config$strategy
    z
  }, .package = 'TensorPLS')
  strategies <- c('M1', 'M2', 'VIP2D_intersect_M1', 'VIP2D_intersect_M2',
    'M1_intersect_M2', 'VIP2D_intersect_M1_intersect_M2',
    'VIP2D_intersect_VIP3Dmodel2')
  for (strategy in strategies) {
    stages <- list()
    cv <- nested_cv_tensorpls(d$X, d$Y, folds = d$folds, outer_folds = 2L, inner_folds = 2L,
      ncomp = 1:2, feature_strategy = strategy, vip_thr = 35,
      imputation = 'training_mean', input_provenance = 'raw', verbose = FALSE)
    expect_length(cv$errors, 0)
    expect_true(all(cv$outer_fold_metrics$status == 'ok'))
    expect_length(stages, 10)
    expect_true(all(unlist(stages) == strategy))
    expect_identical(cv$folds, d$folds)
    expect_true(all(cv$settings$configurations$strategy == strategy))
    expect_true(all(cv$outer_predictions$feature_strategy == strategy))
    stages <- list()
    tuning <- tune_tensorpls_components(d$X, d$Y, folds = d$folds$outer,
      cv_folds = 2L, ncomp = 1:2, feature_strategy = strategy, vip_thr = 35,
      imputation = 'training_mean', verbose = FALSE)
    expect_null(tuning$diagnostics$error)
    expect_length(stages, 5)
    expect_true(all(unlist(stages) == strategy))
    expect_identical(tuning$settings$feature_strategy, strategy)
    expect_identical(tuning$refit$config$strategy, strategy)
  }
})

test_that('new strategy candidates retain tuning metadata and RNG behavior', {
  d <- .m1_m2_data()
  set.seed(19)
  rng <- .Random.seed
  cv <- nested_cv_tensorpls(d$X, d$Y, folds = d$folds, outer_folds = 2L, inner_folds = 2L,
    ncomp = 1:2, feature_strategy = c('M1', 'M2'), strategy_mode = 'tune', vip_thr = 35,
    imputation = 'training_mean', input_provenance = 'raw', verbose = FALSE)
  expect_identical(.Random.seed, rng)
  expect_length(cv$errors, 0)
  expect_setequal(cv$settings$configurations$strategy, c('M1', 'M2'))
  expect_true(all(cv$outer_predictions$feature_strategy %in% c('M1', 'M2')))
  expect_error(nested_cv_tensorpls(d$X, d$Y, feature_strategy = 'M2_intersect_VIP2D',
    imputation = 'training_mean', input_provenance = 'raw', verbose = FALSE), 'Unknown')
  expect_error(tune_tensorpls_components(d$X, d$Y, feature_strategy = 'M2_intersect_VIP2D',
    imputation = 'training_mean', verbose = FALSE), 'supported')
})

test_that('canonical defaults and complete accepted and rejected lists are explicit', {
  accepted <- c('Original', 'VIP2D', 'M1', 'M2', 'VIP3Dmodel2',
    'VIP2D_intersect_M1', 'VIP2D_intersect_M2', 'M1_intersect_M2',
    'VIP2D_intersect_M1_intersect_M2', 'VIP2D_intersect_VIP3Dmodel2')
  removed <- c('VIP3Dmodel1', 'M2_mean', 'VIP2D_intersect_VIP3Dmodel1',
    'VIP2D_intersect_M2_mean', 'VIP3Dmodel1_intersect_VIP3Dmodel2',
    'VIP2D_intersect_VIP3Dmodel1_intersect_VIP3Dmodel2')
  expect_identical(.vip_strategies(), accepted)
  expect_identical(names(.vip_source_fields()), c('VIP2D', 'M1', 'M2', 'VIP3Dmodel2'))
  v <- .m1_m2_vips()
  default <- vip_feature_sets(v)
  expect_identical(default$settings$sources, c('VIP2D', 'M1', 'M2'))
  expect_identical(default$selection$source, c('VIP2D', 'M1', 'M2'))
  expect_identical(names(default$sets), c('VIP2D', 'M1', 'M2',
    'VIP2D_intersect_M1', 'VIP2D_intersect_M2', 'M1_intersect_M2',
    'VIP2D_intersect_M1_intersect_M2'))
  expect_identical(default, vip_feature_sets(v, sources = c('VIP2D', 'M1', 'M2')))
  d <- .m1_m2_data()
  for (s in accepted) expect_identical(vip_feature_sets(v, strategy = s)$strategy, s)
  for (s in removed) {
    expect_error(vip_feature_sets(v, strategy = s), 'Unsupported strategy')
    expect_error(nested_cv_tensorpls(d$X, d$Y, feature_strategy = s), 'Unknown')
    expect_error(tune_tensorpls_components(d$X, d$Y, feature_strategy = s), 'supported')
    expect_error(vip_feature_sets(v, sources = s), 'sources must')
  }
})
