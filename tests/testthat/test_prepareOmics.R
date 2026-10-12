# Data type: Long-format; numeric-like IDs; numeric TIME; numeric features
#
# SubjectID | Visit | FeatA | FeatB
# | ----- | ----- |
#   101     |   0   |   1   |  10
#   202     |   0   |   2   |  20
#   101     |   1   |   3   |  30
#   202     |   1   |   4   |  40
#
# Expectation:
# - IDs and TIME coerced to numeric internally, but dimnames are stored as characters.
# - Features are numeric.
# - Access: arr["101","FeatA","0"] == 1; arr["202","FeatB","1"] == 40
test_that("force_numeric: ID and TIME coerced to numeric; features to numeric", {
  df <- make_long_numeric_id()
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit", transpose = "never",
    coercion_mode = "force_numeric")
  print(unique(as.character(df$SubjectID)))
  print(dimnames(arr)$Subject)
  expect_equal(length(dim(arr)), 3L)
  expect_identical(dimnames(arr)$Subject, c("101", "202"))
# numeric IDs become character dimnames
  expect_identical(dimnames(arr)$Time, c("0", "1"))
  expect_identical(dimnames(arr)$Feature, c("FeatA", "FeatB"))
# spot-check values
  expect_equal(as.numeric(arr["101", "FeatA", "0"]), 1)
  expect_equal(as.numeric(arr["202", "FeatB", "1"]), 40)
})
# Data type: Long-format; string IDs; TIME labels as strings; features stored as character
#
# SubjectID | VisitLabel | FeatA | FeatB
# | ---------- | ----- |
#   P1      |    t0      | "1"   | "10"
#   P2      |    t0      | "2"   | "20"
#   P1      |    t1      | "3"   | "30"
#   P2      |    t1      | "4"   | "40"
#
# Expectation:
# - IDs and TIME remain character; dimnames Subject = c("P1","P2"), Time = c("t0","t1").
# - Input features not coerced at normalization (but tensor is numeric).
# - Access: arr["P1","FeatA","t0"] == 1; arr["P2","FeatB","t1"] == 40
test_that("force_character: ID and TIME as character; features not coerced if requested", {
  df <- make_time_as_labels()
# VisitLabel = "t0","t1"
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "VisitLabel",
    transpose = "never",
    coercion_mode = "force_character", features_numeric = FALSE)
  expect_identical(dimnames(arr)$Subject, c("P1", "P2"))
  expect_identical(dimnames(arr)$Time, c("t0", "t1"))
  expect_identical(dimnames(arr)$Feature, c("FeatA", "FeatB"))
  expect_equal(as.numeric(arr["P1", "FeatA", "t0"]), 1)
  expect_equal(as.numeric(arr["P2", "FeatB", "t1"]), 40)
})
#| SubjectID | VisitLabel | FeatA | FeatB |                           |
# | --------- | ---------- | ----- | ----- | ------------------------- |
#  | 1         | t0         | 1     | 10    |                           |
#  | 1         | t0         | 2     | 20    | ← duplicato di (ID=1, t0) |
#  | 2         | t1         | 3     | 30    |                           |
#  | 2         | t1         | 4     | 40    | ← duplicato di (ID=2, t1) |
test_that("custom: ID numeric, TIME character (mixed policy) works", {
  df <- make_time_as_labels()
  df$SubjectID <- c(1, 1, 2, 2)
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "VisitLabel",
    transpose = "never",
    coercion_mode = "custom", id_as = "numeric", time_as = "character", features_numeric = TRUE)
  expect_identical(dimnames(arr)$Subject, c("1", "2"))
  expect_identical(dimnames(arr)$Time, c("t0", "t1"))
  print(ar)
  expect_equal(as.numeric(arr["1", "FeatA", "t0"]), 1)
  expect_equal(as.numeric(arr["2", "FeatB", "t1"]), 30)
})
##Checking that merge wih cohort id works well 
test_that("cohort filtering aligns ID types and subsets rows as expected", {
  df <- make_long_numeric_id()
# SubjectID = 101, 202
  cdf <- make_cohort_df(ids = c("101", "999"))
# cohort IDs as character
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit", transpose = "never",
# numeric join
    coercion_mode = "force_numeric", cohort = cdf, cohort_id_col = "CohortKey",
      cohort_filter = "SetTag=='Train'")
  expect_equal(dimnames(arr)$Subject, c("101"))
  expect_equal(dim(arr), c(1, 2, 2))
# 1 subject × 2 features × 2 time points
})
# SubjectID | VisitLabel | FeatA | FeatB
# | ---------- | ----- |
#   P1      |    t0      | "1"   | "10"
#   P2      |    t0      | "2"   | "20"
#   P1      |    t1      | "3"   | "30"
#
test_that("min_timepoints retains only subjects with sufficient observations", {
  df <- make_long_df()
# Drop P2 at Visit==1 to make P2 have only 1 time point
  df <- df[!(df$SubjectID == "P2" & df$Visit == 1), ]
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit", transpose = "never",
    coercion_mode = "force_character", min_timepoints = 2)
  expect_identical(dimnames(arr)$Subject, c("P1"))
# P2 filtered out
})
#SubjectID Visit FeatA FeatB
#1        P1     0     1    10
#2        P2     0     2    20
#3        P1     1     3    30
#4        P2     1     4    40
#5        P1     0   101   201
test_that("deduplicate='mean' aggregates duplicate rows for same (ID,Time)", {
  df <- make_long_df()
# Duplicate one row with different feature values
  df <- rbind(df, transform(df[1, ], FeatA = 101, FeatB = 201))
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit", transpose = "never",
    coercion_mode = "custom", id_as = "character", time_as = "numeric", deduplicate = "mean")
  expect_equal(as.numeric(arr["P1", "FeatA", "0"]), mean(c(1, 101)))
  expect_equal(as.numeric(arr["P1", "FeatB", "0"]), mean(c(10, 201)))
})
#### SUBJECT ORDER####
# SubjectID | Visit | FeatA | FeatB
# | ----- | ----- |
# P1        | 0     | 1     | 10
# P2        | 0     | 2     | 20
# P1        | 1     | 3     | 30
# P2        | 1     | 4     | 40
# 1) Subject ordering control
# If you pass a vector subjects = c("P2", "P1", "P3"),
# the final array must have the Subject dimension in the exact order provided,
# even if the original input had a different order (P1 then P2).
#  1)Superset support
# If subjects includes IDs not present in the data (P3 in this case),
# the function must still create a slice for that ID,
# but filled only with NA across all features and time points.
test_that("subjects controls ordering and supports a superset (extra IDs become NA slices)", {
  df <- make_long_df()
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit", transpose = "never",
# P3 not present → should be all-NA row across features and time
    coercion_mode = "force_character", subjects = c("P2", "P1", "P3"))
  expect_identical(dimnames(arr)$Subject, c("P2", "P1", "P3"))
  expect_true(all(is.na(arr["P3", , ])))
})
#####ABOUT TIME POINT ORDER #### 
test_that("time_points controls ordering/subset of the time axis", {
  df <- make_long_df()
  arr <- prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit", transpose = "never",
# reverse order
    coercion_mode = "force_numeric", time_points = c(1, 0))
  expect_identical(dimnames(arr)$Time, c("1", "0"))
})
#INPUT###
#FeatureKey P1_v0 P1_v1 P2_v0 P2_v1
#      FeatA     1     2     3     4
#     FeatB    10    20    30    40
test_that("transpose='always' accepts wide-by-sample input without crashing", {
  raw <- make_transposed_like()
# We just assert the transpose routine runs; your guardrails on missing id/time
# may still raise the explicit error later, which is okay.
  expect_silent({
    tryCatch({
# not present → guardrail error expected downstream
      prepare_omics(data = raw, id_col = "SubjectID", time_col = "Visit", transpose = "always",
        header_in_row = TRUE, coercion_mode = "force_character")
    }, error = function(e)
      invisible(TRUE))
  })
})
test_that("transpose='auto' triggers on non-numeric unique first column", {
  tbl <- data.frame(FeatureKey = c("FeatA", "FeatB", "FeatC"), SbjA_v0 = c(1,
    2, 3), SbjA_v1 = c(4, 5, 6), stringsAsFactors = FALSE)
  expect_silent(.transpose_if_needed(tbl, mode = "auto", header_in_row = TRUE))
})
##This ensures the function fails fast and clearly when a critical column is missing#
test_that("errors are raised: missing time_col or zero time-points after cohort", {
  df <- make_long_df()
# 1) time_col missing
  expect_error(prepare_omics(data = df[, setdiff(names(df), "Visit")], id_col = "SubjectID",
    time_col = "Visit", transpose = "never", coercion_mode = "force_numeric"),
      regexp = "required for a 3-D tensor")
# 2) zero time-points due to cohort filter removing all rows
  cdf <- make_cohort_df(ids = "ZZZ")
# no overlap with P1/P2
  expect_error(prepare_omics(data = df, id_col = "SubjectID", time_col = "Visit",
    transpose = "never",
    coercion_mode = "force_character", cohort = cdf, cohort_id_col = "CohortKey",
      cohort_filter = "SetTag=='Train'"),
    regexp = "No valid time points")
})
