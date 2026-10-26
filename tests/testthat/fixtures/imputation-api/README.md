# Imputation regression fixtures

This directory contains synthetic data and frozen numerical references
used to test the stability of mean and Tucker3 imputation in TensorPLS.

`X_original.rds` contains a named 12 × 8 × 4 synthetic tensor generated
with seed 903, including seven missing entries.

`X_mean_before.rds` and `X_tucker_before.rds` contain the corresponding
imputed tensors. Additional files store fitting controls, reference
outputs, rank-tuning configurations and function definitions used in
regression tests.