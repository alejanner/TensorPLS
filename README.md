# TensorPLS

**PLS-DA for longitudinal omics data in R**

TensorPLS is an R package for analysing longitudinal omics data with PLS-DA. It handles missing measurements, selects molecular features using VIP scores and helps examine their importance across sampling times. Predictive performance is evaluated with subject-level nested cross-validation.

Data are organised as **subjects × features × time**. After missing-data handling, the tensor is unfolded along the subject mode for matrix-based PLS-DA. VIP scores are then mapped back to their original feature–time positions for interpretation.

## Installation

```r
install.packages("remotes")
remotes::install_github("alejanner/TensorPLS")
library(TensorPLS)
```

## Input data

- **`X`**: a numeric 3D array (`subjects × features × time`) with unique, non-empty dimension names. Missing values are represented by `NA`.
- **`Y`**: class labels in the same subject order, supplied as a factor or a `0/1` response matrix. One-hot encoding is supported for multiclass outcomes.

For a binary factor, the second level is the positive class. Longitudinal tables can be converted into an array with `prepare_omics()`.

## Quick start

The following example uses simulated data with a class-associated signal and missing measurements.

```r
library(TensorPLS)
set.seed(123)

n <- 60L
p <- 40L
X <- array(
  rnorm(n * p * 4L),
  dim = c(n, p, 4L),
  dimnames = list(
    paste0("S", seq_len(n)),
    paste0("Gene", seq_len(p)),
    c("T0", "T1", "T2", "T3")
  )
)

Y <- factor(rep(c("Control", "Case"), each = n / 2L),
            levels = c("Control", "Case"))
names(Y) <- dimnames(X)[[1]]

X[Y == "Case", 1:8, 2:4] <- X[Y == "Case", 1:8, 2:4] + 1.2
X[sample.int(length(X), round(0.05 * length(X)))] <- NA_real_
```

### Predictive evaluation

`nested_cv_tensorpls()` evaluates imputation, component tuning, VIP-based feature selection and refitting without using outer-test subjects to develop the model. Here, feature selection retains the intersection of features selected by VIP2D and M1.

```r
selection_strategy <- "VIP2D_intersect_M1"  # VIP2D ∩ M1

cv <- nested_cv_tensorpls(
  X = X,
  Y = Y,
  outer_folds = 5L,
  inner_folds = 5L,
  ncomp = 1:3,
  feature_strategy = selection_strategy,
  vip_thr = 80,
  strategy_mode = "fixed",
  imputation = "training_mean",
  input_provenance = "raw",
  selection_metric = "Q2Y",
  seed = 123L,
  keep_models = FALSE,
  verbose = FALSE
)

stopifnot(length(cv$errors) == 0L)

c(
  Q2Y = cv$Q2Y,
  AUROC = unname(cv$classification_metrics["AUC"]),
  balanced_accuracy = unname(cv$classification_metrics["balanced_accuracy"])
)

head(cv$outer_predictions[
  , c("subject_id", "outer_fold", "Y_true", "Y_pred")
])
```

Each subject contributes one outer held-out prediction. The manuscript applications used VIP2D ∩ M2 with a 99th-percentile threshold; VIP2D ∩ M1 and the 80th percentile here illustrate another available selection option.

### Choose the number of PLS components

For molecular interpretation on the full cohort, first select the number of components from the all-feature model using cross-validation:

```r
development <- tune_tensorpls_components(
  X = X,
  Y = Y,
  ncomp = 1:3,
  cv_folds = 5L,
  feature_strategy = "Original",
  imputation = "training_mean",
  selection_metric = "Q2Y",
  seed = 123L
)

H <- development$selected_ncomp
development$tuning
```

`H` is the selected number of cumulative PLS components. `"Original"` means that this tuning step uses all features, before VIP selection.

### Fit the model and inspect VIP scores

For the full-cohort analysis, fit a model with the selected `H` and calculate VIP scores:

```r
X_complete <- ImputeTensor(X, method = "mean")
model_all <- fit_tensorpls_pls(X_complete, Y, ncomp = H)
vips_all <- tensorpls_vips(model_all)

VIP2D <- vips_all$VIP2D
V <- vips_all$V
M1 <- vips_all$VIP3Dmodel1
M2 <- vips_all$M2_mean

head(VIP2D)
head(M1)
head(M2)
```

### Understanding VIP2D, M1 and M2

**VIP2D** contains the cumulative VIP of each feature–time predictor, with one column for each number of PLS components (`h`). Reshaping these scores gives **`V[f, t, h]`**, the structured VIP representation used in the manuscript. Two summaries are obtained from the same values:

- **`M1[f, h]`** averages `V` over time: how important is a feature across the sampling window for a model with `h` components?
- **`M2[f, t]`** averages `V` over the cumulative component path: how does the feature's importance vary across sampling times?

For example, suppose GeneA has these VIP2D values (illustrative numbers):

| Feature–time predictor | 1 component | 2 components |
| --- | ---: | ---: |
| GeneA at T0 | 1.4 | 1.2 |
| GeneA at T1 | 0.4 | 0.8 |

**M1** averages the two time points for each cumulative model: `M1[GeneA, 1] = (1.4 + 0.4)/2 = 0.9` and `M1[GeneA, 2] = (1.2 + 0.8)/2 = 1.0`. **M2** instead averages the cumulative models at each time: `M2[GeneA, T0] = (1.4 + 1.2)/2 = 1.3` and `M2[GeneA, T1] = (0.4 + 0.8)/2 = 0.6`.

Thus, VIP2D shows individual feature–time scores, M1 summarizes a feature across time, and M2 retains its temporal profile. Here, `h = 2` means the model containing the first **two** components, not the isolated contribution of component 2. In the R API, M1 is returned as `VIP3Dmodel1` and M2 as `M2_mean`.

### Select features and refit

The same selection strategy can be applied to the full-cohort VIP scores:

```r
selected <- vip_feature_sets(
  vips_all,
  strategy = selection_strategy,
  thr = 80
)

selected$selected_features

X_selected <- X_complete[, selected$selected_indices, , drop = FALSE]
model_selected <- fit_tensorpls_pls(X_selected, Y, ncomp = H)
vips_selected <- tensorpls_vips(model_selected)
```

The intersection is taken at the **molecular-feature level**: a feature must pass both VIP2D and M1 selection criteria, then all its time points are retained for model refitting. This is distinct from intersecting individual feature–time positions.

## Plots and selected-feature table

### Score plot

```r
if (H >= 2L) {
  plot_nplsda_scores(model_selected, comp = c(1L, 2L))
} else {
  boxplot(model_selected$x_scores[, 1] ~ Y,
          xlab = "Class", ylab = "PLS component 1 score")
}
```

The score plot shows the fitted subjects in the PLS component space; it is not a plot of held-out predictions. `plot_nplsda_scores()` draws directly, so use `png()`/`dev.off()` if you need to save it.

### Temporal VIP profile

```r
p_time <- plot_temporal_vip(
  vips_selected,
  features = selected$selected_features,
  reference = "none",
  summary = "mean",
  dispersion = "IQR"
)
p_time
ggplot2::ggsave("temporal_vip.png", p_time, width = 7, height = 4)
```

This shows the mean `M2` of the selected features at each sampling time, with the interquartile range across features.

### Final table and VIP heatmap

```r
ranking <- temporal_vip_table(
  vips_selected,
  nested = cv,
  order_by = "mean_M2"
)

head(ranking, 10)
write.csv(ranking, "tensorpls_selected_features.csv", row.names = FALSE)
head(read.csv("tensorpls_selected_features.csv"))

top10 <- head(ranking$feature, 10L)
p_heatmap <- plot_temporal_vip_heatmap(vips_selected, features = top10)
p_heatmap
ggplot2::ggsave("vip_heatmap.png", p_heatmap, width = 7, height = 5)
```

In the exported table, `feature` identifies the selected molecular feature, `mean_M2` is its average M2 across time, `peak_M2` is its highest M2 value and `peak_time` identifies the corresponding sampling time. Columns named `M2_<time>` contain its importance at each time point. With `nested = cv`, `cv_selection_count` and `cv_selection_frequency` report how often the feature was selected across outer training folds. These values describe model-assigned VIP importance, not changes in gene expression or causal effects.

## Missing values

In `nested_cv_tensorpls()` and `tune_tensorpls_components()`, use `imputation = "training_mean"` for training-derived feature–time means, `"tucker"` for train/apply Tucker reconstruction, or `"none"` for complete data. Tucker ranks can be provided through `tucker_control = list(ranks = c(4L, 4L, 3L))`.


## Applications

TensorPLS was evaluated on two public longitudinal transcriptomic datasets with **5 × 5 nested cross-validation**. Metrics were calculated from pooled outer held-out predictions.

| Dataset | Subjects | Features | Time points | Imputation | AUROC | Balanced accuracy | Q²Y |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| [GSE68310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68310) — influenza vs other respiratory viruses | 109 | 47,162 | 4 | Tucker | 0.9715 | 0.9309 | 0.6697 |
| [GSE182616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182616) — burn-injury severity | 66 | 50,599 | 6 | Training mean | 0.8272 | 0.7564 | 0.2939 |


## Documentation

For function arguments and returned objects, see `?tune_tensorpls_components`, `?nested_cv_tensorpls`, `?fit_tensorpls_pls`, `?tensorpls_vips`, `?vip_feature_sets`, `?plot_nplsda_scores`, `?plot_temporal_vip`, `?plot_temporal_vip_heatmap` and `?temporal_vip_table` in R.

## Citation

The TensorPLS manuscript is in preparation. Until it is published, please cite the [TensorPLS repository](https://github.com/alejanner/TensorPLS), specifying the release or commit used.

## License

TensorPLS is distributed under the **GNU General Public License v3 (GPL-3.0)**. See [LICENSE](LICENSE.md).

