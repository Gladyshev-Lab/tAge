# Load AnnData h5ad file and convert to Seurat object

This function loads a single-cell dataset in h5ad format (AnnData) and
converts it to a Seurat object for downstream analysis.

## Usage

``` r
load_h5ad_to_seurat(h5ad_path, verbose = TRUE)
```

## Arguments

- h5ad_path:

  Character string specifying the path to the h5ad file.

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

## Value

A Seurat object containing the single-cell data.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- load_h5ad_to_seurat("data.h5ad")
} # }
```
