test_that("the modern namespace has no retired fitting or VIP bindings", {
  retired <- c("nplsda_vips", "pls_reg", "compute_npls_variates",
    "compute_npls_loadings", "compute_npls_factors",
      "plot_vip2d_with_groups_nogaps", ".validate", ".npls_factors_with_core_arrays",
    ".validate_factors_input_unified", ".compute_global_loadings",
      ".compute_slice_loadings", ".validate_loadings_input_unified",
        ".compute_consensus_loadings_unified", ".validate_variates_input",
          ".compute_global_variates", ".compute_slice_variates",
            ".compute_consensus_variates_exact")
  namespace <- asNamespace("TensorPLS")
  for (symbol in retired) {
    expect_false(exists(symbol, envir = namespace, inherits = FALSE), info = symbol)
    expect_false(symbol %in% getNamespaceExports(namespace), info = symbol)
  }
  for (symbol in c("fit_tensorpls_pls", ".tensorpls_fit", ".tensorpls_predict",
    "tensorpls_vips", "nested_cv_tensorpls")) {
    expect_true(is.function(get(symbol, envir = namespace)), info = symbol)
  }
})
