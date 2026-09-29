# Aggregate single-cell data into pseudobulk samples based on read coverage

This function aggregates single-cell RNA-seq data into pseudobulk
samples by sequentially accumulating cells until a cumulative read
coverage threshold is reached; the next cell starts a new sample. Every
pseudobulk sample holds at least `coverage_threshold` reads, except the
last one, made of the cells left over, which is dropped with
`drop_incomplete = TRUE`. The rule is the one of
[`aggregate_on_obs_columns`](https://gladyshev-lab.github.io/tAge/reference/aggregate_on_obs_columns.md)
and of the Python package's `tage.pp.aggregate`.

## Usage

``` r
aggregate_pseudobulk(
  seurat_obj,
  coverage_threshold = 1e+06,
  assay = "RNA",
  layer = "counts",
  shuffle = FALSE,
  seed = NULL,
  new_sample_prefix = "",
  drop_incomplete = FALSE,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell RNA-seq data.

- coverage_threshold:

  Integer specifying the minimum cumulative read count per pseudobulk
  sample. Cells are sequentially added to a group until this threshold
  is met, then a new group begins. Default is 1e6 (1 million reads).

- assay:

  Character string specifying which assay to use. Default is "RNA".

- layer:

  Character string specifying which layer to use. Default is "counts".

- shuffle:

  Logical indicating whether to randomly shuffle cells before
  aggregation. Shuffling breaks any ordering present in the data.
  Default is FALSE.

- seed:

  Integer seed for reproducibility when shuffle = TRUE. Default is NULL.

- new_sample_prefix:

  Character string prefix for pseudobulk sample names. Default is "".

- drop_incomplete:

  Logical. Drop the last pseudobulk sample when its leftover cells do
  not reach `coverage_threshold`. Default FALSE.

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

## Value

An ExpressionSet object with pseudobulk expression data and metadata
including cumulative_coverage (total reads) and n_cells per pseudobulk
sample.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic coverage-based aggregation
eset <- aggregate_pseudobulk(seurat_obj, coverage_threshold = 1e6)

# With shuffling for randomized grouping
eset <- aggregate_pseudobulk(seurat_obj, coverage_threshold = 5e5,
                             shuffle = TRUE, seed = 42)
} # }
```
