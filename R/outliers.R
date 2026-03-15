#' Remove outlier pseudobulk samples using PCA-based Mahalanobis distance
#'
#' @param eset An ExpressionSet object.
#' @param n_components Integer. Number of PCA components. Default 10.
#' @param threshold_quantile Numeric in (0,1). Chi-squared quantile for outlier
#'   cutoff. Default 0.99.
#' @param split_by Character or NULL. Column in pData to split by before
#'   outlier detection. If NULL, all samples are analyzed together. Default NULL.
#' @param min_samples Integer. Minimum samples in a group to run outlier
#'   detection. Groups below this are kept as-is. Default 10.
#' @param verbose Logical. Default TRUE.
#' @return ExpressionSet with outlier samples removed.
#' @export
remove_outliers <- function(
  eset,
  n_components = 10,
  threshold_quantile = 0.99,
  split_by = NULL,
  min_samples = 10,
  verbose = TRUE
) {
  if (is.null(split_by)) {
    return(.remove_outliers_single(eset, n_components, threshold_quantile, verbose))
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
      if (verbose) cat("  [", grp, "] ", n_sub, " samples — too few, keeping all\n")
      keep_samples <- c(keep_samples, colnames(eset_sub))
      next
    }

    if (verbose) cat("  [", grp, "] ", n_sub, " samples\n")

    eset_filtered <- .remove_outliers_single(
      eset_sub, n_components, threshold_quantile, verbose = verbose
    )

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
  gene_vars <- apply(log_data, 1, var)
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
    cov_mat <- cov(pca_coords) + diag(1e-6, n_components)
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
                        title = "PCA — Outlier Detection")
  }

  if (sum(outlier_mask) == 0) return(eset)
  if (all(outlier_mask)) {
    warning("All samples flagged as outliers, returning original data")
    return(eset)
  }

  eset[, !outlier_mask]
}


.plot_pca_outliers <- function(pca_coords, var_explained, outlier_mask, maha_dist, threshold, title) {
  df <- data.frame(
    PC1 = pca_coords[, 1],
    PC2 = pca_coords[, 2],
    outlier = ifelse(outlier_mask, "Outlier", "Inlier"),
    maha_dist = maha_dist
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = outlier, size = maha_dist)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("Inlier" = "steelblue", "Outlier" = "firebrick")) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggplot2::labs(
      title = title,
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = NULL,
      caption = paste0("Threshold: ", round(threshold, 2),
                        " | Outliers: ", sum(outlier_mask), "/", length(outlier_mask))
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "top",
                    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  print(p)
}
