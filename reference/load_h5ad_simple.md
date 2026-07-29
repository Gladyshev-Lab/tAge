# Load h5ad file using reticulate and anndata

This function provides an alternative method to load h5ad files using
the reticulate package and Python's anndata library, then converts to
Seurat.

## Usage

``` r
load_h5ad_simple(h5ad_path, python_path = NULL, verbose = TRUE)
```

## Arguments

- h5ad_path:

  Character string specifying the path to the h5ad file.

- python_path:

  Character string specifying the path to Python executable. Default is
  NULL (uses current Python environment).

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

## Value

A Seurat object containing the single-cell data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using default Python
seurat_obj <- load_h5ad_simple("data.h5ad")

# Using specific Python environment
seurat_obj <- load_h5ad_simple("data.h5ad", python_path = "./.venv/bin/python")
} # }
```
