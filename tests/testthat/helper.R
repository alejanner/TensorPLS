# Helper functions condivise tra tutti i test
create_synthetic_omics_data <- function(n_subjects = 10, n_features = 2, n_times = 3) {
  data.frame(Individual.Id = rep(1:n_subjects, each = n_times), Time.to.IA = rep(c(-6,
    -3, 0), n_subjects), Feature1 = rnorm(n_subjects * n_times), Feature2 = rnorm(n_subjects *
    n_times))
}
create_large_synthetic_data <- function(n_subjects, n_features, n_times) {
  feature_cols <- paste0("Feature", 1:n_features)
  base_data <- data.frame(Individual.Id = rep(1:n_subjects, each = n_times),
    Time.to.IA = rep(seq(-12,
    0, length.out = n_times), n_subjects))
  feature_data <- matrix(rnorm(n_subjects * n_times * n_features), ncol = n_features)
  colnames(feature_data) <- feature_cols
  cbind(base_data, feature_data)
}
validate_omics_output <- function(arr) {
  is.array(arr) && length(dim(arr)) %in% c(2, 3) && all(is.finite(arr) |
    is.na(arr))
}
