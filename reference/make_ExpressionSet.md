# Create ExpressionSet object from expression data and phenotype data

This function creates an ExpressionSet object from expression data and
phenotype data, which is the standard format used in Bioconductor for
microarray and RNA-seq data.

## Usage

``` r
make_ExpressionSet(exprs_data, phenodata, verbose = TRUE)
```

## Arguments

- exprs_data:

  A data frame or matrix containing expression data where rows are genes
  and columns are samples.

- phenodata:

  A data frame containing phenotype data where rows are samples and
  columns are phenotypic variables.

- verbose:

  Logical indicating whether to print progress messages and create
  density plots. Default is TRUE.

## Value

An ExpressionSet object containing the expression data and phenotype
information.

## Examples

``` r
# Load example data
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()

# Create ExpressionSet
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 
```
