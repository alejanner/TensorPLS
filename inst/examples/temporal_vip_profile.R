# Run: Rscript temporal_vip_profile.R output_directory [saved_model_bundle.rds]
# The optional bundle is the saved public before/after full-data analysis.
library(TensorPLS)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args))
  args[1] else
    "temporal-vip-plots"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (length(args) >= 2L) {
  saved <- readRDS(args[2])
  model <- saved$before
  selected_features <- saved$selected$selected_features
  # Names from the saved public vip_feature_sets(..., strategy=...) result.
  # Both profiles describe this same full-data model, not a subset refit.
  stopifnot(inherits(model, "tensorpls_model"), length(selected_features) > 0L)
} else {
  set.seed(123)
  X <- array(rnorm(40 * 12 * 5), c(40, 12, 5), dimnames = list(paste0("subject",
    1:40), paste0("feature", 1:12), c("-12", "-9", "-6",
    "-3", "0")))
  Y <- factor(rep(c("Control", "Case"), 20))
  model <- fit_tensorpls_pls(X, Y, ncomp = 3)
  # An illustrative user-supplied set; no selection strategy is required.
  selected_features <- c("feature2", "feature4", "feature6", "feature8")
}
vips <- tensorpls_vips(model)
# A: all-feature descriptive profile.
all_summary <- temporal_vip_summary(vips)
p_all <- plot_temporal_vip(vips)
# B: any named signature drawn from the same model's features.
selected_summary <- temporal_vip_summary(vips, features = selected_features, dispersion = "IQR")
p_selected <- plot_temporal_vip(vips, features = selected_features)
p_reference <- plot_temporal_vip(vips, features = selected_features, reference = "all")
p_iqr <- plot_temporal_vip(vips, features = selected_features, dispersion = "IQR")
p_map <- plot_temporal_vip_heatmap(vips, features = selected_features, top_n = 30)
plots <- list(temporal_vip_all = p_all, temporal_vip_selected = p_selected,
  temporal_vip_selected_vs_all = p_reference,
  temporal_vip_selected_iqr = p_iqr, temporal_vip_heatmap = p_map)
for (name in names(plots))
  ggsave(file.path(out, paste0(name, ".png")), plots[[name]],
  width = 7.5, height = 5, units = "in", dpi = 300, bg = "white")
write.csv(all_summary, file.path(out, "all_summary.csv"), row.names = FALSE)
write.csv(selected_summary, file.path(out, "selected_summary.csv"), row.names = FALSE)
write.csv(data.frame(feature = selected_features), file.path(out,
  "selected_features.csv"), row.names = FALSE)
stopifnot(identical(all_summary$importance, unname(colMeans(vips$M2_mean))),
  identical(selected_summary$importance,
  unname(colMeans(vips$M2_mean[selected_features, , drop = FALSE]))))
print(selected_summary)
# Nested CV assesses out-of-sample performance separately. These profiles use
# the final full-data VIP representation, not Q2 or held-out predictions.
