# Perform YuGene normalization

This function performs YuGene normalization on an ExpressionSet. YuGene
is a non-parametric normalization method that transforms expression data
to a uniform distribution by computing cumulative proportions and
applying a specific transformation.

## Usage

``` r
YuGene(eset, verbose = TRUE)
```

## Arguments

- eset:

  An ExpressionSet object containing expression data.

- verbose:

  Logical indicating whether to print progress messages and create
  density plots. Default is TRUE.

## Value

An ExpressionSet object with YuGene-normalized expression data.

## Examples

``` r
# Load example data and create ExpressionSet
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 


# Perform YuGene normalization
yugene_eset <- YuGene(eset, verbose = TRUE)
#> ✓ YuGene normalization completed
```
