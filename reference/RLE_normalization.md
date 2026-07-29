# Perform RLE (Relative Log Expression) normalization

This function performs RLE normalization on an ExpressionSet using the
edgeR package. RLE normalization is commonly used for RNA-seq data to
correct for library size differences and composition bias.

## Usage

``` r
RLE_normalization(original_dataset, verbose = FALSE)
```

## Arguments

- original_dataset:

  An ExpressionSet object containing raw expression data.

- verbose:

  Logical indicating whether to print progress messages and create
  density plots. Default is FALSE.

## Value

An ExpressionSet object with RLE-normalized expression data.

## Examples

``` r
# Load example data and create ExpressionSet
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 


# Perform RLE normalization
rle_eset <- RLE_normalization(eset, verbose = TRUE)
#> calcNormFactors has been renamed to normLibSizes
#> ✓ RLE normalization completed
```
