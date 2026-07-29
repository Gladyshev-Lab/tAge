# Aggregate single-cell data into pseudobulk samples based on read coverage

This function aggregates single-cell RNA-seq data into pseudobulk
samples by sequentially accumulating cells until a cumulative read
coverage threshold is reached. Each pseudobulk sample is guaranteed to
contain at least the specified number of total reads, ensuring
sufficient sequencing depth for downstream transcriptomic age
prediction.

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
