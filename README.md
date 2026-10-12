# TensorPLS

An R package for exploring patterns in three-way omics data using PLS-DA.

## Author
Alessandro Giordano  
PhD Student 

## License
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

You are free to use, modify, and redistribute this software under the terms of this license.  
Any derivative work **must also be released under the same license** and must **credit the original author**.

See the `LICENSE` file for details.

## Citation
If you use **TensorPLS** in academic work, please cite:

Alessandro Giordano (2026). *TensorPLS: An R package for multi-way PLS-DA analysis*.  
GitHub repository: https://github.com/alejanner/TensorPLS



---
# Prerequisites

TensorPLS depends on several R packages.  
Most will be installed automatically with `devtools::install_github()`.  

**Note:** `mixOmics` is hosted on **Bioconductor**, so you need to install it manually first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")
```
# Installation

To install TensorPLS from GitHub, use the [`devtools`](https://cran.r-project.org/package=devtools) package:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install TensorPLS from GitHub
devtools::install_github("alejanner/TensorPLS")
```

## Why TensorPLS?

Multi-omics and longitudinal studies often generate data with a natural **multi-way structure**  
(e.g., *subjects × features × time*).  


---

# 1) Core questions addressed by TensorPLS

When applying PLS-DA with TensorPLS, you can answer:

1. **Is the model able to discriminate groups?**  

2. **Which variables/features mostly contribute to group discrimination?**  

3. **Is the model able to predict outside the sample (Q²)?**  

4. **How much variation is captured (R²)?**  

5. **Which time points or blocks contribute most?**  

---

# 2) Algorithms and data structures used

### Tensors
TensorPLS natively handles **3D arrays** (tensors) that encode multiple modes  
(e.g., *samples × variables × time/blocks*).  

- The **3D structure** is preserved during preprocessing and **imputation of missing values**, ensuring that the multi-way information is respected when filling gaps.  
- For the **modeling step (PLS-DA)**, the tensor is **flattened into a 2D matrix** (samples × unfolded features), since standard PLS-DA operates on matrix data.  

### Partial Least Squares discriminant analysis
PLS finds latent components by **maximizing the covariance between `X` (predictors) and `Y` (response)**, not just variance in `X`.  
In **PLS-DA**, `Y` encodes class membership (e.g., one-hot/dummy coding).

# 3) Prepare a tensor

This small example constructs an incomplete subjects × features × time tensor
from a long-format table. Supply your own data and column names for analysis.

```r
library(TensorPLS)
set.seed(123)
data <- expand.grid(subject = seq_len(12), time = seq_len(4))
for (j in seq_len(8)) data[[paste0("feature", j)]] <- rnorm(nrow(data))
data$feature1[1] <- NA_real_
X <- prepare_omics(data, id_col = "subject", time_col = "time",
                   transpose = "never")
Y <- matrix(rep(c(0, 1), 6), ncol = 1,
            dimnames = list(dimnames(X)[[1]], "Class"))
dim(X)
# 12 8 4
```

## Full-data tensor completion

For descriptive analysis, complete the tensor using all supplied subjects:

```r
X_complete <- ImputeTensor(X, method = "mean")
X_complete <- ImputeTensor(X, method = "tucker")
```

Mean completion uses one mean per feature-time predictor and fills entirely
missing predictors with zero. Tucker uses the fitter's reference ranks
`c(4L, 4L, 3L)` and unchanged controls. These ranks are not universally optimal;
if they exceed the tensor dimensions, supply valid ranks or tune them.

```r
rank_tuning <- tune_tucker_ranks(X)
rank_tuning$selected_ranks
X_complete <- ImputeTensor(
  X, method = "tucker", ranks = rank_tuning$selected_ranks
)

rank_tuning$results
X_complete <- ImputeTensor(X, method = "tucker", ranks = c(3L, 3L, 2L))
```

The default search uses ranks `1:4` per mode, filtered by tensor dimensions and
Tucker validity, with subject rank varying fastest. It contains 37 candidates
when all dimensions are at least four. Recommendation uses the lowest mean
masked RMSE among candidates valid in every repetition; exact ties retain
candidate order. Users may choose any other valid tuple. Custom grids can be
passed as `candidates = my_grid`; `rank_candidates = my_grid` remains supported.
Both `results` and `tuning` expose the same candidate summaries.

Advanced Tucker controls can be supplied through `tucker_control`, using the
same names and defaults as `fit_tucker_imputer()`.

For predictive validation, give incomplete X directly to
`nested_cv_tensorpls(X, Y, imputation = "training_mean")` or its Tucker option.
Do not globally complete X before cross-validation: preprocessing must be
learned within the appropriate training partitions. Full-data completion does
not replace the nested pipeline.
