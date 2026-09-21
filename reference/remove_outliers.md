# Remove outlier samples

Three detectors are available, all working on `log10(counts + 1)` within
each level of `split_by` (tissue, dataset, cell type):

- `"mahalanobis"`:

  (default) robust Mahalanobis distance in PCA space (MCD covariance,
  chi-squared cutoff at `threshold_quantile`) plus PC1/PC2 beyond 3 IQR.

- `"pca_iqr"`:

  the paper's rule for bulk data: samples whose PC1 or PC2 score lies
  more than `iqr_factor` (1.5) interquartile ranges below the first or
  above the third quartile.

- `"spearman_median"`:

  the paper's rule for the integrated meta-dataset: samples whose
  Spearman correlation with the median expression profile of their group
  is below `cor_threshold` (0.5).

## Usage

``` r
remove_outliers(
  eset,
  n_components = 10,
  threshold_quantile = 0.99,
  split_by = NULL,
  min_samples = 10,
  verbose = TRUE,
  method = c("mahalanobis", "pca_iqr", "spearman_median"),
  iqr_factor = 1.5,
  cor_threshold = 0.5
)
```

## Arguments

- eset:

  An ExpressionSet object.

- n_components:

  Integer. Number of PCA components (Mahalanobis). Default 10.

- threshold_quantile:

  Numeric in (0,1). Chi-squared quantile for the Mahalanobis cutoff.
  Default 0.99.

- split_by:

  Character or NULL. Column in pData to split by before outlier
  detection. If NULL, all samples are analyzed together. Default NULL.

- min_samples:

  Integer. Minimum samples in a group to run outlier detection. Groups
  below this are kept as-is. Default 10.

- verbose:

  Logical. Default TRUE.

- method:

  One of `"mahalanobis"`, `"pca_iqr"`, `"spearman_median"`.

- iqr_factor:

  Numeric. Interquartile-range multiplier for `"pca_iqr"`. Default 1.5.

- cor_threshold:

  Numeric. Minimum Spearman correlation with the group median profile
  for `"spearman_median"`. Default 0.5.

## Value

ExpressionSet with outlier samples removed.
