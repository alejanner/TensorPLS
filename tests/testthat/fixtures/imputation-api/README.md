# Exact imputation regression references

These objects were generated before the public imputation API edits, from the
active working tree based on Git commit
`66b558be7d6d931437df0051db1794fa6942f6ae` (with pre-existing uncommitted changes).
The separate verification record preserves that entire working source, its
SHA256 manifest, function definitions, controls and generation scripts:
`/Users/alessandro/Downloads/TensorPLS_imputation_api_verification/before/`.

`X_original.rds` is a named 12 x 8 x 4 Gaussian Tucker tensor generated with
seed 903, with seven missing cells. `controls.rds` records the exact Tucker
controls. The mean and Tucker outputs and missing-cell vectors are frozen
references; do not regenerate them to accommodate a failing test.

`tuning_before.rds` records two explicit candidates and two masks with the same
controls. `default_grid_from_workflow.rds` contains the ordered rank columns
from the existing prepared rank-search workflow's `candidate_grid()`.
`functions.rds` preserves formals and deparsed executable bodies.

The tests require exact identity, without rounding or tolerances. Numeric
references were produced with the R/BLAS environment recorded in the external
baseline's `sessionInfo.txt`; a differing runtime may require investigation,
not replacement of the reference values.
