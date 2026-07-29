# Apply log10 transformation to expression data

This function applies a log10(x + 1) transformation to the expression
data in an ExpressionSet. This transformation is commonly used to
stabilize variance and make the data more suitable for downstream
analysis.

## Usage

``` r
log_transform(eset, verbose = TRUE)
```

## Arguments

- eset:

  An ExpressionSet object containing expression data.

- verbose:

  Logical indicating whether to print progress messages and create
  density plots. Default is TRUE.

## Value

An ExpressionSet object with log10-transformed expression data.

## Examples

``` r
# Load example data and create ExpressionSet
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 


# Apply log transformation
log_eset <- log_transform(eset, verbose = TRUE)
#> ✓ Log transformation completed
```
