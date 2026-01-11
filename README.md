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
``
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

### Partial Least Squares (PLS) & N-PLS-DA
PLS finds latent components by **maximizing the covariance between `X` (predictors) and `Y` (response)**, not just variance in `X`.  
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/PLSDA.png" alt="PLS-DA" width="400">
</p>
In **PLS-DA**, `Y` encodes class membership (e.g., one-hot/dummy coding).  
**N-PLS-DA** extends PLS-DA to **multi-way (tensor) `X`**, extracting components that respect the tensor modes.


# 3) A workflow diagram
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/overviewTensorPLS.png" alt="Pipeline overview" width="400">
</p>

# 4) Tutorial: How to use TensorPLS:

Here we show a minimal workflow using TensorPLS on example data, this example consider a very large dataset with 136 Subjects, 21285 features (Genes) and 5 time points. So we need to 
expect for those type of dataset a significant computational time in Steps like imputation (if needed). 
In this example we have an external metadata file with information about subjects that user can use to filter the dataset (cohort_filter parameter). I will also include some examples from the GCTOF dataset to highlight specific differences when needed, while the complete pipeline is illustrated using the Gene Expression dataset.

```r
library(TensorPLS)

# Load example data
raw_path    <- system.file("extdata", "GeneExpressionDataProcessed.csv", package = "TensorPLS")
cohort_path <- system.file("extdata", "CohortData.csv", package = "TensorPLS")

# Build tensor (Subjects × Features × Time). 
## Parameters used in the example

X <- prepare_omics(
  data          = raw_path,              # Input CSV with omics data
  id_col        = "Individual.Id",       # Subject identifier column
  time_col      = "Time.to.IA",          # Time-point column
  transpose     = "always",              # Force transposition (raw file is wide-by-feature)
  coercion_mode = "force_numeric",       # Coerce IDs and Time to numeric, features to numeric
  cohort        = cohort_path,           # Metadata file with subject annotations
  cohort_id_col = "Group.Id",            # Column in metadata matching subject IDs
  cohort_filter = "Model.or.Validation=='Model'" # Keep only the "Model" subset
)
#> [1]   136 21285     5
```

## Imputation with Tucker decomposition

Now that we have created the tensor `X`, we need to address **missing data**.  
TensorPLS uses a **Tucker decomposition** combined with **ALS (Alternating Least Squares)** for imputation.

On large datasets this step can take hours (2–3h for the Gene Expression dataset).  
On smaller datasets (e.g. metabolomics, 136 × 514 × 5), imputation may only take a few minutes.

---

## Choosing the number of components

Before imputation, we must decide how many components to use for each tensor mode:

- **Mode 1** = Subjects  
- **Mode 2** = Features  
- **Mode 3** = Time  

This decision is made using **Pareto (elbow) plots** and **heatmaps**.

---

### Example: Gene Expression dataset

Pareto plot:  
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/paretoGE.png" alt="Elbow Gene Expression Dataset" width="400">
</p>

We see a good tradeoff at **11 components**.  
From the heatmap below, we select **(4, 4, 3)**.  

<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/heatmapComponents.png" alt="Heatmap Gene Expression" width="400">
</p>

---

### At this step, I also provide a plot from the GCTOF dataset (used for testing), where a more distinct elbow is observed at 10 components.

Pareto plot:  
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/paretoGCTOFX.png" alt="Elbow GCTOF dataset" width="400">
</p>



## Running the imputation

So now we are ready to perform the **imputation step**!  
We need to pass to the imputation function the number of components to be used for the Tucker decomposition — a vital step in the imputation method.

```r
# Gene Expression dataset
fullarrayGeneExpression <- ImputemethodPackage(
  X        = X,
  fac      = c(4, 4, 3),
  conver   = 1e-07,
  max.iter = 1000
)
```
### PLS-DA after imputation  

Once imputation is completed, the next step is **PLS-DA modeling**.  
Before running PLS-DA, we need to **tune the number of components**.  

We use the **Q² metric** to guide this choice:  
- Q² measures the model’s predictive ability.  
- We look for the “elbow”: the point where adding more components no longer increases explained variance.  

 **Important**: PLS-DA requires an outcome `Y` as input.  

#### Structure of Y  
- `Y` must be a **3D array**, just like `X`.  
- **Dimensions:** `Subjects × 1 × Time`  
  - **Subjects** → same number of individuals as in `X`  
  - **1** → only one feature (the outcome/class)  
  - **Time** → same number of time points as in `X`  

 **Example**:  
If `X` has shape `136 × 21285 × 5`, then `Y` must have shape `136 × 1 × 5`.  

Here `Y` encodes class membership, e.g.:  
- `0 = control`  
- `1 = case`  

---

#### Tuning the number of components

```r
nCompGE <- ncomp_elbow_nplsda(fullarrayGeneExpression, outcomedummyarray136, reps = 10)
```
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/nCompGE.png" alt="Elbow pls-DA" width="400">
</p>

We see a good tradeoff at **3 components**.  

### PLS-DA analysis  

Now that we are ready for **PLS-DA analysis**, we can use several functions in TensorPLS to answer the core questions:  

- Is the model able to discriminate groups?  
- How much variance is explained (**R²**)?  
- How well does the model predict unseen data (**Q²**)?  
- Which features contribute the most (**VIP scores**)?  

---

#### 1) Running N-PLS-DA with VIPs

```r
nplsda_vipsGE <- nplsda_vips(
  X = fullarrayGeneExpression, 
  Y = outcomedummyarray136,
  ncomp = 3,               # number of components chosen in the tuning step
  slice_vip = TRUE
)



```
Function outputs

This function outputs:

Q² → predictive ability of the model

R² → explained variance

VIP scores → feature importance metrics

### Example: cumulative Q² values with 3 components

| Component | Q² (t1) | Q² (t2) | Q² (t3) | Q² (t4) | Q² (t5) |
|-----------|---------|---------|---------|---------|---------|
| 1         | 0.1980  | 0.1980  | 0.1980  | 0.1980  | 0.1980  |
| 2         | 0.5184  | 0.5184  | 0.5184  | 0.5184  | 0.5184  |
| 3         | 0.7721  | 0.7721  | 0.7721  | 0.7721  | 0.7721  |

The model reaches about **77% predictive power (Q²)**.

Explained variance (R²):
| Component | R2X    | R2Xcum | R2Y    | R2Ycum |
| --------- | ------ | ------ | ------ | ------ |
| t1        | 0.0674 | 0.0674 | 0.4631 | 0.4631 |
| t2        | 0.0918 | 0.1592 | 0.2836 | 0.7466 |
| t3        | 0.0580 | 0.2172 | 0.1601 | 0.9068 |

R2Ycum shows that the cumulative explained variance of Y is about 90%.

### Understanding VIPs

The concept of **VIP (Variable Importance in Projection)** is central to PLS and N-PLS-DA.  

TensorPLS provides three complementary views:

| View            | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| **VIP2D**       | For component *h*, how important is feature *f* at time *t*.               |
| **VIP3D Model 1** | On component *h*, how important is feature *f* on average across time.    |
| **VIP3D Model 2** | At time *t*, how important is feature *f* overall across components.      |

---

### Assessing group discrimination

To visualize how well the model separates groups, we compute **variates**:

```r
nplsda_vipsVariatesGE= compute_npls_variates(X = fullarrayGeneExpression, Y = outcomedummyarray136, ncomp =3)
###How to define classe vec:
class_vec <- factor(
  outcomedummyarray136[, 1, 1],
  levels = c(0, 1),
  labels = c("Class0", "Class1")
)
head(class_vec)
229251 235421 249696 254394 259207 265155 
Class0 Class0 Class1 Class0 Class1 Class1 
Levels: Class0 Class1
```

Then we plot the scores to check group separation:
```r
plot_nplsda_scores(scores_matrix = nplsda_vipsVariatesGE$NPLSDAvariates$X,nplsdaVipIntersectionGETotal$explvar, class_vec = class_vec,pc1 = 1,pc2 = 2,variance = TRUE)
```
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/tvariatesGE.png" alt="Mo pls-DA" width="400">
</p>
Each dot represents a sample, and the two colors correspond to the two classes.

We observe a partial separation between groups: even without feature selection, the model is already able to discriminate between classes.

Now, to improve the ability of the model to discriminate groups we have added a Feature Selection Steps. 
### Feature Selection to improve discrimination

To improve the ability of the model to discriminate groups, we add a **Feature Selection step**.  
The logic of feature selection in TensorPLS is based on a combination of two approaches:

1. **Feature selection on the VIPs** objects returned from the `nplsda_vips` function.  
2. Evaluate the **intersection or union** between selected variables.  
3. **Re-run PLS-DA** on these different combinations and compare results in terms of metrics like **Q²** and **R²**.  
4. Choose the **best model**.  

---

> At the moment, only **point 1** is implemented as a function.  
> Points **2–4** must be carried out **manually** by the user.  
> An automatic pipeline for this process is not yet available, so users must compare results across different situations.  

---

### Example: Feature selection workflow

Here we show an example of how to perform feature selection:  

---

#### 1. Feature selection on the VIPs

We apply the `feature_selection` function on the three VIP objects returned from `nplsda_vips`.

```r
selectedVip3DModel2 <- feature_selection(nplsda_vipsGE$VIP3Dmodel2, thr = 99, strip_time = TRUE)
selectedVip3DModel1 <- feature_selection(nplsda_vipsGE$VIP3Dmodel1, thr = 99, strip_time = TRUE)
selectedVip2D       <- feature_selection(
  nplsda_vipsGE$VIP2D,
  thr = 99,
  strip_time = TRUE,
  pattern = "_-?\\d+$"
)
```
> **Note**  
> - The output of **VIP2D** is formatted as `Feature_Time`, e.g. `Feature_-n` where *n* is one of the time points (*-12, -9, -6, -3, 0*).  
> - To remove this suffix and compare features across time, use the **regex pattern** (`pattern`) together with `strip_time = TRUE`.  
> - You may need to **adapt the pattern** depending on your data format.


Now that we have completed the **first step of feature selection**, the next step is to inspect the **relationships between the different VIPs** selected at the chosen percentile threshold.  

You can use a preliminary function to visualize the overlap between selected features:  

```r
feature_lists <- list(
  "VIP 3D – Model 2" = unique(selectedVip3DModel2),
  "VIP 3D – Model 1" = unique(selectedVip3DModel1),
  "VIP 2D – Model"   = unique(selectedVip2D)
)

draw_venn_diagram(feature_lists)

```
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/vennGE.png" alt="Mo pls-DA" width="400">
</p>

From the Venn diagram, we can observe that **VIP 3D Model 2** and **VIP 2D** share the largest overlap.  

However, this information alone is **not sufficient** to decide whether to use their **union** or **intersection**.  
To make this decision, it is necessary to **re-run PLS-DA** on each combination of selected features and compare performance metrics such as **Q²** and **R²**.  

---

### Example: Re-running PLS-DA with different feature sets

In this example, we show how to compare models obtained from different **intersections** or **unions** of selected features.  
The first step is to **recreate the tensor (`fullarray`)** using only the selected variables.

```r
# Taking intersection between VIP 3D Model 2 and VIP 3D Model 1
intersectionGEModel1Model2 <- intersect(selectedVip3DModel2, selectedVip3DModel1)

# Create the full array containing only those features
fullarrayIntersectionGEModel1Model2 <- fullarrayGeneExpression[,
  is.element(colnames(fullarrayGeneExpression), intersectionGEModel1Model2),
]

# Tune the number of components (if needed) and re-run N-PLS-DA
nplsda_vipsGEModel1Model2 <- nplsda_vips(
  fullarrayIntersectionGEModel1Model2,
  outcomedummyarray136,
  ncomp = 3
)
```
You can repeat this process for each combination of intersections and unions between VIP sets.
### Best performing model

In our case, the **best model** was obtained using the **intersection between VIP2D and VIP 3D Model 2**,  
reaching the following performance:

- **Q² = 0.9397**  
- **R² = 0.98**

---

### Variates after feature selection

Re-plotting the variates with this feature set, we can observe an almost **complete separation between groups**:

<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/VIP2DModel32.png" alt="Variates after Feature Selection" width="400">
</p>

### Time-course contribution analysis

Up to now, we have shown how **feature selection** improves discrimination and predictive performance.  
But in **time-course data**, it’s not only *which variables* matter — it’s also **when they matter**.  
In other words: *which time points are driving the model the most?*  

To answer this, we need to **compute the factors**.  

> **Note**: The number of components (`ncomp`) must always be the same as the one used in the PLS-DA analysis.

```r
factorsGE <- compute_npls_factors(
  X = fullarrayIntersectionGEModel1Model2,
  Y = outcomedummyarray136,
  ncomp = 3
)

plot_nplsda_blockX_mode3(factorsGE, edge = c(0.2, 0.3))
```
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/GEfactorsIntersection.png" alt="Tines that drives the components" width="400">
</p>

### Interpreting time-point contributions

What do we see in the plot?  
Each point represents a **time point** projected onto the first two components.

- Points **far from the origin (0,0)** → contribute strongly to shaping the components, i.e. they are **influential for class separation**.  
- Points **close to the origin** → have little impact on the components.  
- **Distance from the origin ∝ importance** of the time point.  
- **Direction** indicates which component is being driven.  
- **Sign (positive/negative)** reflects the *phase*: time points with opposite signs drive the components in opposite directions.  
- **Clustered points** suggest time points with a similar role.  
- **Outliers** indicate key time points strongly influencing the model.  

---

### Plotting VIP2D and Assessing Feature Robustness

Another  functionality offered by TensorPLS is the ability to **plot VIP2D scores** and understand what they represent.  
VIP2D highlights the **importance of each feature at each time point** in driving class separation.  

In addition to simple ranking, TensorPLS provides a way to **evaluate the robustness** of these Top-N VIP2D features through **permutation testing**.  
This allows you to distinguish between features that appear influential in a single model fit and those that remain stable across repeated resampling.

---

#### How it works

1. The function extracts the **Top-N features (per time point)** with the highest VIP2D scores.  
2. Class labels (`Y`) are **randomly permuted** across samples, and the PLS-DA model is re-fitted.  
3. At each permutation, the **Top-N list** is recomputed.  
4. For every feature, the **number of times** it reappears in the Top-N list across permutations is counted.  
5. A **permutation p-value** is computed as:  

where *R* is the total number of permutations.  

---

#### Interpretation

-  A **low permutation p-value** (e.g., < 0.05) → the feature consistently ranks among the Top-N and is unlikely to appear by chance.  
-  A **high permutation p-value** → the feature is unstable and less reliable as a biomarker.  

---

#### In the plots

- **Orange bars + orange dots** = features selected in the Top-N by VIP score.  
- **Green dots** = features also significant according to permutation testing (robust features).  

This visualization helps distinguish between:  
- **Influential but unstable features** (orange only)  
- **Robust features** that remain stable across resampling (green).  

---

#### Example usage

```r
plot_vip2d_with_groups_nogaps(
  X      = fullarrayIntersectionGEModel1Model2,
  Y      = outcomedummyarray136,
  ncomp  = 3,
  topN   = 20,
  perms  = 1000
)
```
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/plotVip2DPer.png" alt="Tines that drives the components" width="400">
</p>
