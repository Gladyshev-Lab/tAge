# Test helpers ---------------------------------------------------------------

# Silence the density plots the preprocessing functions emit.
suppressMessages(requireNamespace("Biobase", quietly = TRUE))

# Build the example ExpressionSet shipped with the package.
.tage_example_eset <- function() {
  make_ExpressionSet(
    load_example_expression_data(),
    load_example_metadata(),
    verbose = FALSE
  )
}

# Download (and cache) a clock model from Zenodo for prediction tests.
# Returns the local path, or NULL if the download failed / no network.
.tage_test_model <- function(filename,
                             record = "18763485") {
  cache <- tools::R_user_dir("tAge", "cache")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cache, filename)
  if (!file.exists(dest)) {
    # download_clocks() raises the timeout for the gigabyte Bayesian ridge
    # models and removes partial or non-pickle files itself.
    tryCatch(
      download_clocks(filename, dest_dir = cache, record = record, quiet = TRUE,
                      timeout = 120),   # enough for the ~1 MB elastic net models; BR models are skipped
      error = function(e) NULL, warning = function(w) NULL
    )
  }
  if (file.exists(dest) && file.info(dest)$size > 0) dest else NULL
}

# Skip a test unless reticulate and the required Python modules are available.
skip_if_no_python <- function() {
  testthat::skip_if_not_installed("reticulate")
  ok <- tryCatch(
    all(vapply(c("joblib", "pandas", "sklearn", "numpy"),
               reticulate::py_module_available, logical(1))),
    error = function(e) FALSE
  )
  if (!isTRUE(ok)) {
    testthat::skip("Python modules (joblib/pandas/scikit-learn/numpy) not available")
  }
}
