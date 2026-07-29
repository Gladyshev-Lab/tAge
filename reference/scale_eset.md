# Scale expression data to have zero mean and unit variance

This function applies z-score scaling to the expression data in an
ExpressionSet. It calls [`scale`](https://rdrr.io/r/base/scale.html),
which operates column-wise, so scaling is performed *per sample* (each
sample/column is scaled to zero mean and unit variance across genes).
This is the "Scaling" normalisation strategy from the paper and matches
the TACO reference application. Per-gene standardisation is handled
separately inside the trained clock model (its `StandardScaler` step),
using training-set statistics.

## Usage

``` r
scale_eset(eset, verbose = TRUE)
```

## Arguments

- eset:

  An ExpressionSet object containing expression data.

- verbose:

  Logical indicating whether to print progress messages and create
  density plots. Default is TRUE.

## Value

An ExpressionSet object with scaled expression data.

## Examples

``` r
# Load example data and create ExpressionSet
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 


# Scale expression data
scaled_eset <- scale_eset(eset, verbose = TRUE)
#> ✓ Scaling completed
```
