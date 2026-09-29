# flag_outliers(): the rule of the Python package's tage.pp.flag_outliers.

.profile <- function(seed, n_genes = 3000) {
  set.seed(seed)
  p <- stats::rlnorm(n_genes, 0, 2)
  p / sum(p)
}

# Negative binomial counts (dispersion 0.05) around one expression profile.
.counts <- function(profile, n, seed, depth = 2e7) {
  set.seed(seed)
  lib <- depth * stats::rlnorm(n, 0, 0.3)
  mu <- outer(profile, lib)
  matrix(stats::rpois(length(mu), stats::rgamma(length(mu), shape = 20, scale = mu / 20)),
         nrow = length(profile))
}

.eset <- function(m, ...) {
  colnames(m) <- paste0("s", seq_len(ncol(m)))
  rownames(m) <- paste0("g", seq_len(nrow(m)))
  make_ExpressionSet(m, data.frame(row.names = colnames(m), ...), verbose = FALSE)
}

test_that("normal samples are not flagged and nothing is removed", {
  e <- flag_outliers(.eset(.counts(.profile(1), 12, 2)))
  expect_false(any(e$tage_outlier))
  expect_equal(unname(ncol(e)), 12L)
  expect_true(all(c("tage_outlier_distance", "tage_outlier_z") %in% colnames(Biobase::pData(e))))
})

test_that("a sample from another tissue is flagged", {
  m <- .counts(.profile(1), 12, 3)
  m[, 5] <- .counts(.profile(9), 1, 4)[, 1]
  e <- flag_outliers(.eset(m))
  expect_equal(which(e$tage_outlier), 5L)
})

test_that("depth does not make an outlier", {
  m <- .counts(.profile(1), 12, 5)
  deeper <- m; deeper[, 1] <- deeper[, 1] * 3
  expect_equal(flag_outliers(.eset(m))$tage_outlier_distance,
               flag_outliers(.eset(deeper))$tage_outlier_distance)
  set.seed(6)
  shallow <- m; shallow[, 1] <- stats::rbinom(nrow(m), m[, 1], 0.2)
  expect_false(flag_outliers(.eset(shallow))$tage_outlier[1])
})

test_that("strata are scored separately; a small stratum is reported, not scored", {
  m <- cbind(.counts(.profile(1), 8, 7), .counts(.profile(9), 8, 8), .counts(.profile(1), 3, 9))
  e <- expect_warning(
    flag_outliers(.eset(m, tissue = rep(c("A", "B", "C"), c(8, 8, 3))), split_by = "tissue"),
    "tissue = C: 3 samples"
  )
  expect_false(any(e$tage_outlier))
  expect_true(all(is.na(e$tage_outlier_z[e$tissue == "C"])))
  expect_false(anyNA(e$tage_outlier_z[e$tissue != "C"]))
})

test_that("the paper's PCA rule is available", {
  m <- .counts(.profile(1), 12, 10)
  m[, 1] <- .counts(.profile(9), 1, 11)[, 1]
  e <- flag_outliers(.eset(m), method = "pca_iqr")
  expect_true(e$tage_outlier[1])
  expect_true(all(c("tage_outlier_pc1", "tage_outlier_pc2") %in% colnames(Biobase::pData(e))))
})

test_that("the Klotho example flags the samples the Python package flags", {
  # Cross-checked with tage.pp.flag_outliers on the same data.
  e <- flag_outliers(.tage_example_eset(), split_by = "Tissue")
  expect_equal(Biobase::sampleNames(e)[e$tage_outlier], c("RNA_93M", "RNA_99M"))
  expect_equal(max(e$tage_outlier_z[e$Tissue == "Kidney"]), 1.3, tolerance = 0.05)
})

test_that("input is checked and the removed function points to the new one", {
  m <- .counts(.profile(1), 8, 12)
  bad <- m; bad[1, 1] <- -1
  expect_error(flag_outliers(.eset(bad)), "raw counts")
  expect_error(flag_outliers(.eset(m), split_by = "nope"), "not a column")
  expect_error(flag_outliers(.eset(m), method = "union"))
  expect_error(remove_outliers(.eset(m)), "flag_outliers")
})
