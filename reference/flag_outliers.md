# Flag outlier samples

Flags samples within each stratum; nothing is removed. The rule is the
one of the Python package's `tage.pp.flag_outliers`, and both give the
same numbers.

## Usage

``` r
flag_outliers(
  eset,
  split_by = NULL,
  method = c("distance", "pca_iqr"),
  threshold = 5,
  iqr_factor = 1.5,
  min_cpm = 10,
  min_samples = 6
)
```

## Arguments

- eset:

  ExpressionSet of raw counts (bulk or pseudobulk).

- split_by:

  Column of the phenoData defining the strata (tissue, dataset, cell
  type) that are scored separately, as in
  [`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md).
  `NULL`: one stratum.

- method:

  `"distance"` or `"pca_iqr"`.

- threshold:

  Robust z-score above which `"distance"` flags a sample. Default 5.

- iqr_factor:

  Interquartile ranges beyond the quartiles for `"pca_iqr"`. Default
  1.5.

- min_cpm:

  Genes with a lower mean CPM in the stratum are left out. Default 10.

- min_samples:

  Strata with fewer samples are not scored (a warning says which).
  Default 6.

## Value

The ExpressionSet with phenoData columns `tage_outlier` (logical) and,
for `"distance"`, `tage_outlier_distance` and `tage_outlier_z`; for
`"pca_iqr"`, `tage_outlier_pc1` and `tage_outlier_pc2`. Samples of
unscored strata get `FALSE` and `NA`.

## Details

`method = "distance"` (default) scores every sample by `1 - r`, the
Pearson correlation of its log2(CPM + 1) profile with the median profile
of its stratum, and flags it when the robust z-score of that distance
(median and MAD within the stratum) exceeds `threshold`. It is the
paper's "correlation with the median profile of the tissue" rule with a
threshold relative to the stratum instead of a fixed rho \< 0.5, which
passes even a sample from another tissue. Only samples further from the
median than the rest are flagged.

`method = "pca_iqr"` is the paper's rule: PC1 or PC2 more than
`iqr_factor` interquartile ranges beyond the quartiles. The paper
applied it to hundreds of samples; with a dozen, the first components
pick out the noisiest sample and 1.5 IQR flags a normal sample in most
strata.

Both work on log2(CPM + 1) of the genes with a mean CPM of at least
`min_cpm` in the stratum. Depth only moves a sample's log profile by a
constant, which the correlation ignores, but low-count genes add noise
that grows as depth falls, and PCA on unnormalised counts mostly finds
depth.

Remove the flagged samples explicitly, after looking at them:
`eset <- eset[, !eset$tage_outlier]`.

## Examples

``` r
eset <- make_ExpressionSet(load_example_expression_data(), load_example_metadata(),
                           verbose = FALSE)
eset <- flag_outliers(eset, split_by = "Tissue")
Biobase::pData(eset)[eset$tage_outlier, c("Tissue", "tage_outlier_z")]
#>                  Tissue tage_outlier_z
#> RNA_93M Skeletal muscle       9.443985
#> RNA_99M Skeletal muscle       7.983690
eset <- eset[, !eset$tage_outlier]
```
