# Test unitari per .read_or_pass()
test_that(".read_or_pass returns object unchanged when not a file path", {
  # Test con data.frame
  df <- data.frame(x = 1:3, y = 4:6)
  expect_identical(.read_or_pass(df), df)
  # Test con matrice
  mat <- matrix(1:6, nrow = 2)
  expect_identical(.read_or_pass(mat), mat)
  # Test con lista
  lst <- list(a = 1, b = 2)
  expect_identical(.read_or_pass(lst), lst)
})
test_that(".read_or_pass returns string unchanged when file doesn't exist", {
  nonexistent_file <- "this_file_does_not_exist.csv"
  expect_identical(.read_or_pass(nonexistent_file), nonexistent_file)
})
test_that(".read_or_pass handles multiple strings", {
  # Deve restituire il vettore identico se lunghezza > 1
  multiple_strings <- c("file1.csv", "file2.csv")
  expect_identical(.read_or_pass(multiple_strings), multiple_strings)
})
test_that(".read_or_pass reads CSV files correctly", {
  # Crea file CSV temporaneo
  temp_csv <- tempfile(fileext = ".csv")
  test_data <- data.frame(id = 1:3, value = c(1.1, 2.2, 3.3))
  write.csv(test_data, temp_csv, row.names = FALSE)
  result <- .read_or_pass(temp_csv)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  # Cleanup
  unlink(temp_csv)
})
test_that(".read_or_pass reads RDS files correctly", {
  # Crea file RDS temporaneo
  temp_rds <- tempfile(fileext = ".rds")
  test_data <- list(x = 1:5, y = letters[1:5])
  saveRDS(test_data, temp_rds)
  result <- .read_or_pass(temp_rds)
  expect_equal(result, test_data)
  # Cleanup
  unlink(temp_rds)
})
test_that(".read_or_pass fails gracefully with unsupported extensions", {
  # Crea file con estensione non supportata
  temp_unknown <- tempfile(fileext = ".xyz")
  writeLines("test content", temp_unknown)
  expect_error(.read_or_pass(temp_unknown), "Unsupported file extension")
  # Cleanup
  unlink(temp_unknown)
})
test_that(".read_or_pass handles arrow dependencies correctly", {
  # Test che requireNamespace funzioni
  # Questo test potrebbe essere skippato se arrow non è disponibile
  skip_if_not_installed("arrow")
  # Se arrow è disponibile, testa feather
  temp_feather <- tempfile(fileext = ".feather")
  test_data <- data.frame(x = 1:3, y = 4:6)
  # Solo se arrow è davvero disponibile
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_feather(test_data, temp_feather)
    result <- .read_or_pass(temp_feather)
    expect_s3_class(result, "data.frame")
    unlink(temp_feather)
  }
})
