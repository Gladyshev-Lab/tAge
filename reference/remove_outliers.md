# Remove outlier pseudobulk samples using PCA-based Mahalanobis distance

Remove outlier pseudobulk samples using PCA-based Mahalanobis distance

## Usage

``` r
remove_outliers(
  eset,
  n_components = 10,
  threshold_quantile = 0.99,
  split_by = NULL,
  min_samples = 10,
  verbose = TRUE
)
```

## Arguments

- eset:

  An ExpressionSet object.

- n_components:

  Integer. Number of PCA components. Default 10.

- threshold_quantile:

  Numeric in (0,1). Chi-squared quantile for outlier cutoff. Default
  0.99.

- split_by:

  Character or NULL. Column in pData to split by before outlier
  detection. If NULL, all samples are analyzed together. Default NULL.

- min_samples:

  Integer. Minimum samples in a group to run outlier detection. Groups
  below this are kept as-is. Default 10.

- verbose:

  Logical. Default TRUE.

## Value

ExpressionSet with outlier samples removed.
