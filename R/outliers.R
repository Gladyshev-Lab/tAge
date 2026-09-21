#' Remove outlier samples
#'
#' Three detectors are available, all working on \code{log10(counts + 1)}
#' within each level of \code{split_by} (tissue, dataset, cell type):
#' \describe{
#'   \item{\code{"mahalanobis"}}{(default) robust Mahalanobis distance in
#'     PCA space (MCD covariance, chi-squared cutoff at
#'     \code{threshold_quantile}) plus PC1/PC2 beyond 3 IQR.}
#'   \item{\code{"pca_iqr"}}{the paper's rule for bulk data: samples whose
#'     PC1 or PC2 score lies more than \code{iqr_factor} (1.5) interquartile
#'     ranges below the first or above the third quartile.}
#'   \item{\code{"spearman_median"}}{the paper's rule for the integrated
#'     meta-dataset: samples whose Spearman correlation with the median
#'     expression profile of their group is below \code{cor_threshold}
#'     (0.5).}
#' }
#'
#' @param eset An ExpressionSet object.
#' @param n_components Integer. Number of PCA components (Mahalanobis). Default 10.
#' @param threshold_quantile Numeric in (0,1). Chi-squared quantile for the
#'   Mahalanobis cutoff. Default 0.99.
#' @param split_by Character or NULL. Column in pData to split by before
#'   outlier detection. If NULL, all samples are analyzed together. Default NULL.
#' @param min_samples Integer. Minimum samples in a group to run outlier
#'   detection. Groups below this are kept as-is. Default 10.
#' @param verbose Logical. Default TRUE.
#' @param method One of \code{"mahalanobis"}, \code{"pca_iqr"},
#'   \code{"spearman_median"}.
#' @param iqr_factor Numeric. Interquartile-range multiplier for
#'   \code{"pca_iqr"}. Default 1.5.
#' @param cor_threshold Numeric. Minimum Spearman correlation with the group
#'   median profile for \code{"spearman_median"}. Default 0.5.
#' @return ExpressionSet with outlier samples removed.
#' @export
remove_outliers <- function(
  eset,
  n_components = 10,
  threshold_quantile = 0.99,
  split_by = NULL,
  min_samples = 10,
  verbose = TRUE,
  method = c("mahalanobis", "pca_iqr", "spearman_median"),
  iqr_factor = 1.5,
  cor_threshold = 0.5
) {
  method <- match.arg(method)
  single <- function(e) {
    switch(
      method,
      mahalanobis     = .remove_outliers_single(e, n_components, threshold_quantile, verbose),
      pca_iqr         = .remove_outliers_pca_iqr(e, iqr_factor, verbose),
      spearman_median = .remove_outliers_spearman(e, cor_threshold, verbose)
    )
  }

  if (is.null(split_by)) {
    return(single(eset))
  }

  if (!split_by %in% colnames(Biobase::pData(eset))) {
    stop(paste0("'", split_by, "' not found in pData"))
  }

  groups <- unique(Biobase::pData(eset)[[split_by]])

  if (verbose) {
    cat("Removing outliers by group:", split_by, "\n")
    cat("  - Groups:", length(groups), "\n\n")
  }

  keep_samples <- character()

  for (grp in groups) {
    mask <- Biobase::pData(eset)[[split_by]] == grp
    eset_sub <- eset[, mask]
    n_sub <- ncol(eset_sub)

    if (n_sub < min_samples) {
      if (verbose) cat("  [", grp, "] ", n_sub, " samples \u2014 too few, keeping all\n")
      keep_samples <- c(keep_samples, colnames(eset_sub))
      next
    }

    if (verbose) cat("  [", grp, "] ", n_sub, " samples\n")

    eset_filtered <- single(eset_sub)

    n_removed <- n_sub - ncol(eset_filtered)
    if (verbose && n_removed > 0) {
      cat("    - Removed:", n_removed, "outliers\n")
    } else if (verbose) {
      cat("    - No outliers\n")
    }

    keep_samples <- c(keep_samples, colnames(eset_filtered))
  }

  eset_final <- eset[, keep_samples]

  if (verbose) {
    n_total_removed <- ncol(eset) - ncol(eset_final)
    pct <- round(n_total_removed / ncol(eset) * 100, 1)
    cat("\n\u2713 Total: removed", n_total_removed, "of", ncol(eset),
        "samples (", pct, "%)\n")
    cat("  - Samples retained:", ncol(eset_final), "\n")
  }

  return(eset_final)
}


.remove_outliers_single <- function(
  eset,
  n_components = 10,
  threshold_quantile = 0.99,
  verbose = TRUE
) {
  expr_data <- Biobase::exprs(eset)
  n_samples <- ncol(expr_data)
  n_genes <- nrow(expr_data)

  if (n_samples < 6) {
    if (verbose) cat("    Too few samples for outlier detection, skipping.\n")
    return(eset)
  }

  # Clamp n_components: need n_samples > n_components + 1, and keep ratio reasonable
  max_comp <- max(2, floor((n_samples - 1) / 2))
  max_comp <- min(max_comp, n_genes - 1)
  if (n_components > max_comp) {
    if (verbose) cat("    Clamping PCA components from", n_components, "to", max_comp, "\n")
    n_components <- max_comp
  }

  # PCA on log-transformed counts
  log_data <- log10(expr_data + 1)
  gene_vars <- apply(log_data, 1, stats::var)
  log_data <- log_data[gene_vars > 0, , drop = FALSE]

  pca_result <- prcomp(t(log_data), center = TRUE, scale. = TRUE, rank. = n_components)
  pca_coords <- pca_result$x[, seq_len(n_components), drop = FALSE]
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

  # Robust Mahalanobis distance using MCD estimator
  if (requireNamespace("robustbase", quietly = TRUE) && n_samples >= 2 * n_components + 1) {
    mcd <- robustbase::covMcd(pca_coords, alpha = 0.75)
    center <- mcd$center
    cov_mat <- mcd$cov + diag(1e-6, n_components)
  } else {
    center <- colMeans(pca_coords)
    cov_mat <- stats::cov(pca_coords) + diag(1e-6, n_components)
  }

  maha_dist <- mahalanobis(pca_coords, center, cov_mat)

  # Chi-squared threshold
  threshold <- qchisq(threshold_quantile, df = n_components)
  outlier_maha <- maha_dist > threshold

  # Additionally flag samples with extreme PC1/PC2 values (simple IQR method)
  outlier_iqr <- rep(FALSE, n_samples)
  for (pc in 1:min(2, n_components)) {
    vals <- pca_coords[, pc]
    q1 <- quantile(vals, 0.25)
    q3 <- quantile(vals, 0.75)
    iqr <- q3 - q1
    outlier_iqr <- outlier_iqr | (vals < q1 - 3 * iqr) | (vals > q3 + 3 * iqr)
  }

  # Union of both methods
  outlier_mask <- outlier_maha | outlier_iqr

  if (verbose) {
    n_maha <- sum(outlier_maha)
    n_iqr <- sum(outlier_iqr & !outlier_maha)
    if (n_maha > 0 || n_iqr > 0) {
      cat("    - Mahalanobis outliers:", n_maha, "\n")
      if (n_iqr > 0) cat("    - Additional IQR outliers:", n_iqr, "\n")
    }

    .plot_pca_outliers(pca_coords, var_explained, outlier_mask, maha_dist, threshold,
                        title = "PCA outlier detection")
  }

  if (sum(outlier_mask) == 0) return(eset)
  if (all(outlier_mask)) {
    warning("All samples flagged as outliers, returning original data")
    return(eset)
  }

  eset[, !outlier_mask]
}


# Log-scaled, variance-filtered PCA scores shared by the PCA-based detectors.
.tage_log_pca <- function(expr_data, rank) {
  log_data <- log10(expr_data + 1)
  gene_vars <- apply(log_data, 1, stats::var)
  log_data <- log_data[gene_vars > 0, , drop = FALSE]
  stats::prcomp(t(log_data), center = TRUE, scale. = TRUE, rank. = rank)
}

# The paper's bulk-data rule: PC1 or PC2 beyond Q1 - k*IQR / Q3 + k*IQR.
.remove_outliers_pca_iqr <- function(eset, iqr_factor = 1.5, verbose = TRUE) {
  expr_data <- Biobase::exprs(eset)
  n_samples <- ncol(expr_data)
  if (n_samples < 4) {
    if (verbose) cat("    Too few samples for outlier detection, skipping.\n")
    return(eset)
  }
  pca <- .tage_log_pca(expr_data, rank = 2)
  scores <- pca$x[, 1:2, drop = FALSE]
  var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100

  outlier <- rep(FALSE, n_samples)
  for (pc in 1:2) {
    v  <- scores[, pc]
    q  <- stats::quantile(v, c(0.25, 0.75))
    iqr <- q[2] - q[1]
    outlier <- outlier | v < q[1] - iqr_factor * iqr | v > q[2] + iqr_factor * iqr
  }

  if (verbose) {
    cat("    - PC1/PC2 IQR outliers:", sum(outlier), "\n")
    .plot_pca_outliers(scores, var_explained, outlier, rep(1, n_samples),
                       threshold = iqr_factor, title = "PCA outlier detection (IQR rule)")
  }
  if (!any(outlier)) return(eset)
  if (all(outlier)) {
    warning("All samples flagged as outliers, returning original data")
    return(eset)
  }
  eset[, !outlier]
}

# The paper's meta-dataset rule: Spearman correlation of each sample with the
# median profile of its group below cor_threshold.
.remove_outliers_spearman <- function(eset, cor_threshold = 0.5, verbose = TRUE) {
  expr_data <- Biobase::exprs(eset)
  n_samples <- ncol(expr_data)
  if (n_samples < 3) {
    if (verbose) cat("    Too few samples for outlier detection, skipping.\n")
    return(eset)
  }
  log_data <- log10(expr_data + 1)
  ref <- apply(log_data, 1, stats::median, na.rm = TRUE)
  rho <- apply(log_data, 2, function(x) {
    suppressWarnings(stats::cor(x, ref, method = "spearman", use = "pairwise.complete.obs"))
  })
  outlier <- is.na(rho) | rho < cor_threshold

  if (verbose) {
    cat("    - Spearman correlation with the group median profile: min",
        round(min(rho, na.rm = TRUE), 3), "; outliers:", sum(outlier), "\n")
  }
  if (!any(outlier)) return(eset)
  if (all(outlier)) {
    warning("All samples flagged as outliers, returning original data")
    return(eset)
  }
  eset[, !outlier]
}


.plot_pca_outliers <- function(pca_coords, var_explained, outlier_mask, maha_dist, threshold, title) {
  df <- data.frame(
    PC1 = pca_coords[, 1],
    PC2 = pca_coords[, 2],
    outlier = ifelse(outlier_mask, "Outlier", "Inlier"),
    maha_dist = maha_dist
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2, fill = .data$outlier, size = .data$maha_dist)) +
    ggplot2::geom_point(alpha = 0.9, shape = 21, colour = TAGE_SURFACE, stroke = 0.5) +
    ggplot2::scale_fill_manual(values = c("Inlier" = TAGE_REFERENCE_COLOR, "Outlier" = "#e34948")) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      fill = NULL,
      caption = paste0("threshold ", round(threshold, 2),
                        " \u00b7 outliers ", sum(outlier_mask), "/", length(outlier_mask))
    ) +
    theme_tage(base_size = 10, grid = "both")

  print(p)
}
