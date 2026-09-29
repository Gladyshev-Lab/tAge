# Pseudobulk aggregation: cells are added in order until a sample holds at
# least the coverage threshold, the rule of tage.pp.aggregate in Python.

test_that("cells are grouped until the threshold is reached, then a new group starts", {
  expect_equal(tAge:::.coverage_groups(c(60, 60, 60, 60), 100), c(1L, 1L, 2L, 2L))
  # a cell above the threshold closes its group; no group is left empty
  expect_equal(tAge:::.coverage_groups(c(50, 250, 30), 100), c(1L, 1L, 2L))
  expect_equal(tAge:::.coverage_groups(c(100, 100), 100), c(1L, 2L))
})

.toy_counts <- function(n_cells = 200, seed = 1) {
  set.seed(seed)
  m <- matrix(stats::rpois(50 * n_cells, 3), nrow = 50,
              dimnames = list(paste0("gene", 1:50), paste0("cell", seq_len(n_cells))))
  Matrix::Matrix(m, sparse = TRUE)
}

test_that("every pseudobulk but the last reaches the threshold and counts are kept", {
  counts <- .toy_counts()
  res <- tAge:::.aggregate_matrix(counts, 1000, shuffle = FALSE, seed = NULL)
  k <- length(res$coverages)
  expect_true(all(res$coverages[-k] >= 1000))
  expect_true(all(res$cell_counts > 0))
  expect_equal(sum(res$cell_counts), ncol(counts))
  expect_equal(unname(colSums(res$counts)), res$coverages)
  expect_equal(unname(rowSums(res$counts)), unname(Matrix::rowSums(counts)))

  dropped <- tAge:::.aggregate_matrix(counts, 1000, shuffle = FALSE, seed = NULL,
                                      drop_incomplete = TRUE)
  expect_true(all(dropped$coverages >= 1000))
  expect_equal(length(dropped$coverages), sum(res$coverages >= 1000))
})

test_that("aggregate_pseudobulk and aggregate_on_obs_columns follow the same rule", {
  testthat::skip_if_not_installed("Seurat")
  counts <- .toy_counts()
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  obj$donor <- rep(c("d1", "d2"), each = 100)

  pb <- aggregate_pseudobulk(obj, coverage_threshold = 1000, verbose = FALSE)
  cov <- Biobase::pData(pb)$cumulative_coverage
  expect_true(all(cov[-length(cov)] >= 1000))
  expect_true(all(Biobase::pData(pb)$n_cells > 0))
  expect_equal(sum(Biobase::exprs(pb)), sum(counts))

  by_donor <- aggregate_on_obs_columns(obj, "donor", coverage_threshold = 1000,
                                       drop_incomplete = TRUE, verbose = FALSE)
  expect_true(all(Biobase::pData(by_donor)$cumulative_coverage >= 1000))
  expect_setequal(unique(Biobase::pData(by_donor)$donor), c("d1", "d2"))
})
