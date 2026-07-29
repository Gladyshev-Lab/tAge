# Subset Seurat object by metadata criteria

Convenience function to subset Seurat objects with multiple criteria.

## Usage

``` r
subset_seurat_by_metadata(seurat_obj, filters, min_cells = 1, verbose = TRUE)
```

## Arguments

- seurat_obj:

  A Seurat object to subset.

- filters:

  Named list of filtering criteria. Each element should be a metadata
  column name with a vector of values to keep.

- min_cells:

  Minimum number of cells to retain after filtering. If fewer cells
  remain, returns NULL with a warning. Default is 1.

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

## Value

A subsetted Seurat object, or NULL if fewer than min_cells remain.

## Examples

``` r
if (FALSE) { # \dontrun{
# Subset by tissue and age
seurat_subset <- subset_seurat_by_metadata(
  seurat_obj,
  filters = list(
    tissue = c("Kidney", "Liver"),
    age = c("3m", "24m")
  )
)

# Subset by multiple criteria
seurat_subset <- subset_seurat_by_metadata(
  seurat_obj,
  filters = list(
    tissue = "Tongue",
    sex = "male",
    age = "24m"
  ),
  min_cells = 100
)
} # }
```
