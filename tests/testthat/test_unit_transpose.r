# Test unitari per .transpose_if_needed()
test_that(".transpose_if_needed mode='never' returns unchanged dataframe", {
  df <- data.frame(label = c("A", "B", "C"), x = 1:3, y = 4:6)
  result <- .transpose_if_needed(df, mode = "never")
  expect_identical(result, df)
})
test_that(".transpose_if_needed mode='always' transposes dataframe", {
  df <- data.frame(label = c("A", "B", "C"), x = 1:3, y = 4:6)
  result <- .transpose_if_needed(df, mode = "always", header_in_row = TRUE)
# Dopo trasposizione: 3 colonne (A, B, C), 2 righe (x, y)
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), 2)
  expect_equal(colnames(result), c("A", "B", "C"))
})
test_that(".transpose_if_needed mode='always' without header_in_row", {
  df <- data.frame(label = c("A", "B"), x = 1:2)
  result <- .transpose_if_needed(df, mode = "always", header_in_row = FALSE)
# Senza header_in_row, mantiene tutte le righe
  expect_equal(nrow(result), 2)
# label + x
  expect_equal(ncol(result), 2)
# 2 osservazioni originali
})
test_that(".transpose_if_needed mode='auto' detects when to transpose", {
# Caso che DOVREBBE essere trasposto: prima colonna character con valori unici
# Character, tutti unici
  df_should_transpose <- data.frame(id = c("Subject1", "Subject2", "Subject3"), x = 1:3,
    y = 4:6)
  result_auto <- .transpose_if_needed(df_should_transpose, mode = "auto")
  expect_equal(ncol(result_auto), 3)
# Trasposto
# Caso che NON dovrebbe essere trasposto: prima colonna numerica
# Numerica
  df_no_transpose <- data.frame(id = 1:3, x = 4:6, y = 7:9)
  result_no_auto <- .transpose_if_needed(df_no_transpose, mode = "auto")
  expect_identical(result_no_auto, df_no_transpose)
# Non trasposto
})
test_that(".transpose_if_needed mode='auto' handles duplicates correctly", {
# Prima colonna character ma con duplicati → NON trasporre
# Duplicato!
  df_duplicates <- data.frame(id = c("A", "A", "B"), x = 1:3,
    y = 4:6)
  result <- .transpose_if_needed(df_duplicates, mode = "auto")
  expect_identical(result, df_duplicates)
# Non dovrebbe trasporre
})
test_that(".transpose_if_needed handles non-dataframe input", {
# Test con matrice
  mat <- matrix(1:6, nrow = 2)
  result <- .transpose_if_needed(mat, mode = "always")
  expect_identical(result, mat)
# Dovrebbe restituire identico
# Test con NULL
  result_null <- .transpose_if_needed(NULL, mode = "always")
  expect_null(result_null)
})
test_that(".transpose_if_needed handles edge cases", {
# Data.frame vuoto
  df_empty <- data.frame()
  result_empty <- .transpose_if_needed(df_empty, mode = "auto")
  expect_identical(result_empty, df_empty)
# Data.frame con una sola colonna
  df_single <- data.frame(x = 1:3)
  result_single <- .transpose_if_needed(df_single, mode = "auto")
# Prima colonna è numerica, non dovrebbe trasporre
  expect_identical(result_single, df_single)
# Data.frame con una sola riga
  df_one_row <- data.frame(a = 1, b = 2, c = 3)
  result_one_row <- .transpose_if_needed(df_one_row, mode = "always")
  expect_equal(nrow(result_one_row), 2)
# 3 colonne → 2 righe (dopo rimozione header)
})
