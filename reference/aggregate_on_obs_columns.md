# Aggregate single-cell data into pseudobulk samples within obs column groups

This function first splits the data by all unique combinations of the
specified observation (metadata) columns, then applies coverage-based
aggregation within each group. This is useful for creating pseudobulk
samples stratified by biological variables such as sample, cell type,
tissue, or condition.

## Usage

``` r
aggregate_on_obs_columns(
  seurat_obj,
  obs_column_names,
  coverage_threshold = 1e+06,
  assay = "RNA",
  layer = "counts",
  shuffle = FALSE,
  seed = NULL,
  new_sample_prefix = "",
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell RNA-seq data.

- obs_column_names:

  Character vector of metadata column names to stratify by. For example,
  c("sample_id", "cell_type") will first split cells into groups defined
  by each unique sample_id x cell_type combination, then aggregate
  within each group.

- coverage_threshold:

  Integer specifying the minimum cumulative read count per pseudobulk
  sample. Default is 1e6.

- assay:

  Character string specifying which assay to use. Default is "RNA".

- layer:

  Character string specifying which layer to use. Default is "counts".

- shuffle:

  Logical indicating whether to randomly shuffle cells before
  aggregation within each group. Default is FALSE.

- seed:

  Integer seed for reproducibility when shuffle = TRUE. Default is NULL.

- new_sample_prefix:

  Character string prefix for pseudobulk sample names. Default is "".

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

## Value

An ExpressionSet object with pseudobulk expression data. Metadata
includes the stratification columns, cumulative_coverage, and n_cells.

## Examples

``` r
if (FALSE) { # \dontrun{
# Aggregate within each sample
eset <- aggregate_on_obs_columns(seurat_obj,
                                 obs_column_names = "sample_id",
                                 coverage_threshold = 1e6)

# Aggregate within each sample x cell_type combination
eset <- aggregate_on_obs_columns(seurat_obj,
                                 obs_column_names = c("sample_id", "cell_type"),
                                 coverage_threshold = 5e5)
} # }
```
