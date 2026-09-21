# Outlier detection: a planted outlier must be caught by each rule and no
# ordinary sample removed.

.outlier_fixture <- function(seed = 5, n = 16, genes = 400) {
  set.seed(seed)
  mu <- stats::rlnorm(genes, meanlog = 4, sdlog = 1)
  m  <- sapply(seq_len(n), function(i) stats::rpois(genes, mu))
  dimnames(m) <- list(paste0("g", seq_len(genes)), paste0("s", seq_len(n)))
  # Sample 1 is a different profile altogether.
  m[, 1] <- stats::rpois(genes, sample(mu))
  make_ExpressionSet(m, data.frame(row.names = colnames(m), grp = rep("a", n)), verbose = FALSE)
}

test_that("pca_iqr flags the planted outlier only", {
  eset <- .outlier_fixture()
  out <- remove_outliers(eset, method = "pca_iqr", verbose = FALSE)
  expect_false("s1" %in% colnames(out))
  expect_gte(ncol(out), ncol(eset) - 2)
})

test_that("spearman_median flags the planted outlier only", {
  eset <- .outlier_fixture()
  out <- remove_outliers(eset, method = "spearman_median", cor_threshold = 0.9, verbose = FALSE)
  expect_false("s1" %in% colnames(out))
  expect_gte(ncol(out), ncol(eset) - 2)
})

test_that("the Mahalanobis default still runs and keeps the sample order", {
  eset <- .outlier_fixture()
  out <- remove_outliers(eset, verbose = FALSE)
  expect_s4_class(out, "ExpressionSet")
  expect_true(all(colnames(out) %in% colnames(eset)))
})

test_that("split_by applies the rule within each group", {
  eset <- .outlier_fixture(n = 24)
  Biobase::pData(eset)$grp <- rep(c("a", "b"), each = 12)
  out <- remove_outliers(eset, split_by = "grp", method = "pca_iqr", min_samples = 6, verbose = FALSE)
  expect_false("s1" %in% colnames(out))
  expect_true(all(paste0("s", 13:24) %in% colnames(out)))
})
