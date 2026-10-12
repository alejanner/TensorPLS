# test_plot_nplsda_scores.R
# Unit tests for plot_nplsda_scores function
create_mock_data <- function(n_obs = 50, n_comp = 3, n_classes = 2) {
  list(scores = matrix(rnorm(n_obs * n_comp), nrow = n_obs, ncol = n_comp),
    expl_var = matrix(c(0.127,
    0.127, 0.180, 0.180, 0.077, 0.204, 0.220, 0.401, 0.039, 0.243,
    0.315, 0.716), nrow = 3, ncol = 4, byrow = TRUE, dimnames = list(c("t1",
    "t2", "t3"), c("R2X", "R2Xcum", "R2Y", "R2Ycum"))), class_vec = factor(sample(0:(n_classes -
    1), n_obs, replace = TRUE), levels = 0:(n_classes - 1), labels = paste0("Class",
    0:(n_classes - 1))))
}
test_that("Function works with basic inputs", {
  data <- create_mock_data()
  # Should not throw error
  expect_silent({
    pdf(NULL)  # Null device to avoid creating actual plot
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec)
    dev.off()
  })
})
test_that("Function handles different component selections", {
  data <- create_mock_data(n_comp = 5)
  # Test different component combinations
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec,
      pc1 = 1, pc2 = 2)
    dev.off()
  })
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec,
      pc1 = 2, pc2 = 3)
    dev.off()
  })
})
test_that("Function handles multiple classes correctly", {
  # Test with 3 classes
  data <- create_mock_data(n_classes = 3)
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec)
    dev.off()
  })
  # Test with 5 classes
  data <- create_mock_data(n_classes = 5)
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec)
    dev.off()
  })
})
test_that("Variance parameter works correctly", {
  data <- create_mock_data()
  # Test with variance = TRUE
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec,
      variance = TRUE)
    dev.off()
  })
  # Test with variance = FALSE
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec,
      variance = FALSE)
    dev.off()
  })
})
test_that("Function fails with invalid inputs", {
  data <- create_mock_data()
  # Test with out of bounds component index
  expect_error({
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec,
      pc1 = 10  # Out of bounds
)
  })
  # Test with mismatched dimensions
  wrong_class_vec <- factor(c("A", "B"))  # Wrong length
  expect_error({
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = wrong_class_vec)
  })
})
test_that("Custom title is applied correctly", {
  data <- create_mock_data()
  custom_title <- "My Custom NPLS-DA Plot"
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
      class_vec = data$class_vec,
      main_title = custom_title)
    dev.off()
  })
})
test_that("Function handles edge cases", {
  # Minimum case: 2 observations, 2 components, 2 classes
  minimal_data <- create_mock_data(n_obs = 2, n_comp = 2, n_classes = 2)
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = minimal_data$scores, expl_var = minimal_data$expl_var,
      class_vec = minimal_data$class_vec)
    dev.off()
  })
  # Single class (edge case - should still work but not very meaningful)
  single_class_data <- create_mock_data(n_classes = 1)
  expect_silent({
    pdf(NULL)
    plot_nplsda_scores(scores_matrix = single_class_data$scores,
      expl_var = single_class_data$expl_var,
      class_vec = single_class_data$class_vec)
    dev.off()
  })
})
test_that("Function preserves graphical parameters", {
  data <- create_mock_data()
  # Save original par
  pdf(NULL)
  original_mar <- par("mar")
  original_xpd <- par("xpd")
  plot_nplsda_scores(scores_matrix = data$scores, expl_var = data$expl_var,
    class_vec = data$class_vec)
  # Check that parameters are restored
  expect_equal(par("mar"), original_mar)
  expect_equal(par("xpd"), original_xpd)
  dev.off()
})
