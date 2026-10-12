# Test unitari per .create_array3d()
test_that(".create_array3d creates correct dimensions", {
  test_data <- data.frame(id = c(1, 1, 2, 2), time = c(-6,
    -3, -6, -3), feat1 = c(1.1, 1.2, 2.1, 2.2),
    feat2 = c(3.1, 3.2, 4.1, 4.2))
  result <- .create_array3d(test_data, id_col = "id", time_col = "time")
# Dimensioni: 2 soggetti, 2 feature, 2 tempi
  expect_equal(dim(result), c(2, 2, 2))
  expect_equal(length(dim(result)), 3)
})
test_that(".create_array3d has correct dimension names", {
  test_data <- data.frame(subject_id = c("A", "A", "B", "B"), timepoint = c(-6,
    -3, -6, -3), gene1 = c(1.1, 1.2, 2.1, 2.2),
    gene2 = c(3.1, 3.2, 4.1, 4.2))
  result <- .create_array3d(test_data, id_col = "subject_id", time_col = "timepoint")
  expect_equal(names(dimnames(result)), c("Subject", "Feature", "Time"))
  expect_equal(dimnames(result)$Subject, c("A", "B"))
  expect_equal(dimnames(result)$Feature, c("gene1", "gene2"))
  expect_equal(dimnames(result)$Time, c("-6", "-3"))
})
test_that(".create_array3d handles explicit feature columns", {
  test_data <- data.frame(id = c(1, 1, 2, 2), time = c(-6,
    -3, -6, -3), feat1 = c(1.1, 1.2, 2.1, 2.2),
    feat2 = c(3.1, 3.2, 4.1, 4.2), feat3 = c(5.1, 5.2, 6.1,
# Extra feature
      6.2))
# Specifica solo feat1 e feat3
  result <- .create_array3d(test_data, id_col = "id", time_col = "time", feature_cols = c("feat1",
    "feat3"))
  expect_equal(dim(result)[2], 2)
# Solo 2 feature
  expect_equal(dimnames(result)$Feature, c("feat1", "feat3"))
})
test_that(".create_array3d handles custom subjects and time_points", {
  test_data <- data.frame(id = c(1, 1, 2, 2, 3, 3),
    time = c(-6, -3, -6, -3, -6, -3), feat = c(1.1,
      1.2, 2.1, 2.2, 3.1, 3.2))
# Specifica solo soggetti 1 e 3, e solo tempo -6
  result <- .create_array3d(test_data, id_col = "id", time_col = "time", subjects = c(1,
    3), time_points = c(-6))
  expect_equal(dim(result), c(2, 1, 1))
# 2 soggetti, 1 feature, 1 tempo
  expect_equal(dimnames(result)$Subject, c("1", "3"))
  expect_equal(dimnames(result)$Time, "-6")
})
test_that(".create_array3d handles duplicates with different strategies", {
# Dati con duplicati per stesso soggetto/tempo
# Stesso soggetto, stesso tempo
  test_data <- data.frame(id = c(1, 1, 1), time = c(-6, -6,
    -6), feat = c(1.0, 1.1, 1.2))
# Test "first"
  result_first <- .create_array3d(test_data, "id", "time", deduplicate = "first")
  expect_equal(result_first[1, 1, 1], 1.0)
# Test "last"
  result_last <- .create_array3d(test_data, "id", "time", deduplicate = "last")
  expect_equal(result_last[1, 1, 1], 1.2)
# Test "mean"
  result_mean <- .create_array3d(test_data, "id", "time", deduplicate = "mean")
  expect_equal(result_mean[1, 1, 1], 1.1)
# Media di 1.0, 1.1, 1.2
})
test_that(".create_array3d handles missing data correctly", {
# Soggetto 2 ha solo un timepoint
  test_data <- data.frame(id = c(1, 1, 2), time = c(-6, -3,
    -6), feat = c(1.1, 1.2, 2.1))
  result <- .create_array3d(test_data, "id", "time")
# Soggetto 2 dovrebbe avere NA al tempo -3
  expect_true(is.na(result[2, 1, 2]))
# Soggetto 2, feature 1, tempo -3
  expect_equal(result[2, 1, 1], 2.1)
# Soggetto 2, feature 1, tempo -6
})
test_that(".create_array3d produces numeric array", {
  test_data <- data.frame(id = c(1, 2), time = c(-6, -6), feat = c("1.1",
# Character numbers
    "2.2"))
  result <- .create_array3d(test_data, "id", "time")
  expect_true(is.numeric(result))
  expect_equal(storage.mode(result), "double")
})
