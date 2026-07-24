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
    url <- sprintf("https://zenodo.org/records/%s/files/%s?download=1",
                   record, filename)
    ok <- tryCatch(
      utils::download.file(url, dest, mode = "wb", quiet = TRUE) == 0,
      error = function(e) FALSE, warning = function(e) FALSE
    )
    if (!isTRUE(ok) && file.exists(dest)) unlink(dest)
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
