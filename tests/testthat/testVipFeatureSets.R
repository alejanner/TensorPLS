vip_test_input <- function() {
  make <- function(selected, rn = LETTERS[1:5], cn = "Comp1") {
    matrix(as.numeric(LETTERS[1:5] %in% selected), 5, 1, dimnames = list(rn,
      cn))
  }
  list(VIP2D = make(c("A", "B", "C"), paste0(LETTERS[1:5], "_T1")),
    VIP3Dmodel1 = make(c("B", "C", "D")), VIP3Dmodel2 = make(c("C",
      "D", "E"), cn = "T1"))
}
test_that("individuals, every intersection and unions are exact", {
  v <- vip_test_input()
  before <- v
  fs <- vip_feature_sets(v, thr = 60, combinations = c("individual", "intersection", "union"))
  expect_s3_class(fs, "vip_feature_sets")
  expect_equal(fs$sets$VIP2D, c("A", "B", "C"))
  expect_equal(fs$sets$VIP3Dmodel1, c("B", "C", "D"))
  expect_equal(fs$sets$VIP3Dmodel2, c("C", "D", "E"))
  expect_equal(fs$sets$VIP2D_intersect_VIP3Dmodel1, c("B", "C"))
  expect_equal(fs$sets$VIP2D_intersect_VIP3Dmodel2, "C")
  expect_equal(fs$sets$VIP3Dmodel1_intersect_VIP3Dmodel2, c("C", "D"))
  expect_equal(fs$sets$VIP2D_intersect_VIP3Dmodel1_intersect_VIP3Dmodel2, "C")
  expect_equal(fs$sets$VIP2D_union_VIP3Dmodel1, c("A", "B", "C", "D"))
  expect_equal(fs$sets$VIP2D_union_VIP3Dmodel2, LETTERS[1:5])
  expect_equal(fs$sets$VIP3Dmodel1_union_VIP3Dmodel2, LETTERS[2:5])
  expect_equal(fs$sets$VIP2D_union_VIP3Dmodel1_union_VIP3Dmodel2, LETTERS[1:5])
  expect_equal(fs$summary$n_features, unname(lengths(fs$sets)))
  expect_equal(fs$summary$sources[[4]], c("VIP2D", "VIP3Dmodel1"))
  expect_equal(fs$summary$operation[4], "intersection")
  expect_identical(v, before)
})
test_that("threshold, cols and aggregation use feature_selection semantics", {
  v <- vip_test_input()
  v$VIP3Dmodel1 <- cbind(a = 1:5, b = 5:1)
  rownames(v$VIP3Dmodel1) <- LETTERS[1:5]
  for (a in c("any", "all", "mean"))
    for (t in c(0, 50, 99, 100)) {
    fs <- vip_feature_sets(v, thr = t, sources = "VIP3Dmodel1", combinations = "individual",
      cols = 1:2, aggregator = a)
    expect_equal(fs$sets[[1]], feature_selection(v$VIP3Dmodel1, thr = t, cols = 1:2,
      aggregator = a))
    expect_equal(fs$settings$thr, t)
  }
  expect_equal(vip_feature_sets(v, thr = 75, sources = "VIP3Dmodel1", combinations = "individual",
    cols = "a")$sets[[1]], c("D", "E"))
  expect_error(vip_feature_sets(v, cols = "absent"), "unknown column")
})
test_that("suffix mapping and duplicate identifiers are handled without renaming", {
  v <- vip_test_input()
  v$VIP2D <- matrix(1:6, 6, 1, dimnames = list(c("A|t1", "A|t2",
    "B|t1", "B|t2", "C|t1", "C|t2"), "Comp1"))
  fs <- vip_feature_sets(v, thr = 0, pattern = "\\|t\\d+$")
  expect_equal(fs$sets$VIP2D, c("A", "B", "C"))
  expect_equal(fs$selection$n_collapsed_rows[1], 3)
  expect_equal(vip_feature_sets(v, thr = 0, strip_time = FALSE)$sets$VIP2D,
    rownames(v$VIP2D))
  expect_error(vip_feature_sets(v), "inconsistent.*pattern")
  v$VIP3Dmodel1 <- matrix(c(1, 2, 3), 3, 1, dimnames = list(c("A",
    "A", "B"), "Comp1"))
  expect_equal(vip_feature_sets(v, thr = 0, sources = "VIP3Dmodel1",
    combinations = "individual")$sets[[1]],
    c("A", "B"))
  expect_equal(feature_selection(v$VIP3Dmodel1, thr = 0), c("A", "B"))
  one <- matrix(3, 1, 1, dimnames = list("A", "c1"))
  expect_equal(feature_selection(one), "A")
  # Gene names with a time-looking suffix must not be truncated automatically.
  v <- vip_test_input()
  rownames(v$VIP3Dmodel1) <- rownames(v$VIP3Dmodel2) <- paste0(LETTERS[1:5],
    "_T1")
  rownames(v$VIP2D) <- paste0(LETTERS[1:5], "_T1_T1")
  expect_equal(vip_feature_sets(v, thr = 60)$sets$VIP2D, c("A_T1", "B_T1", "C_T1"))
})
test_that("empty intersections are represented and malformed scores rejected", {
  v <- vip_test_input()
  v$VIP2D[, 1] <- c(5, 0, 0, 0, 0)
  v$VIP3Dmodel1[, 1] <- c(0, 5, 0, 0, 0)
  fs <- vip_feature_sets(v, thr = 100)
  expect_equal(fs$sets$VIP2D_intersect_VIP3Dmodel1, character())
  expect_equal(fs$summary$n_features[4], 0)
  for (bad in list(-1, 101, NA_real_, Inf, "95", c(20, 30)))
    expect_error(vip_feature_sets(v,
    thr = bad), "percentile")
  expect_error(vip_feature_sets(v[-1]), "missing.*VIP2D")
  expect_error(vip_feature_sets(NULL), "missing")
  bad <- v
  bad$VIP2D <- c("A", "B")
  expect_error(vip_feature_sets(bad), "numeric matrix")
  bad <- v
  rownames(bad$VIP2D) <- NULL
  expect_error(vip_feature_sets(bad), "row names")
  bad <- v
  bad$VIP2D[1, 1] <- Inf
  expect_error(vip_feature_sets(bad), "infinite")
  bad <- v
  bad$VIP3Dmodel1[, ] <- NA_real_
  expect_error(vip_feature_sets(bad), "missing/empty")
  expect_error(vip_feature_sets(bad, na.rm = FALSE), "missing values")
  expect_error(vip_feature_sets(v, combinations = "bad"), "combinations")
  expect_error(vip_feature_sets(v, sources = c("VIP2D", "VIP2D")), "distinct")
  expect_error(vip_feature_sets(v, strip_time = NA), "strip_time")
  expect_error(vip_feature_sets(v, pattern = c("a", "b")), "pattern")
})
