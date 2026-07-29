# Filter genes based on expression thresholds

This function filters genes from an ExpressionSet based on count and
percentage thresholds. Genes are retained if they have expression values
above the count threshold in at least the specified percentage of
samples.

## Usage

``` r
filter_genes(
  exprs_set,
  count_threshold = 10,
  percent_threshold = 20,
  verbose = TRUE
)
```

## Arguments

- exprs_set:

  An ExpressionSet object containing expression data.

- count_threshold:

  Numeric threshold for minimum expression count. Default is 10.

- percent_threshold:

  Numeric threshold for minimum percentage of samples that must have
  expression above count_threshold. Default is 20.

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

## Value

A filtered ExpressionSet object with only genes passing the thresholds.

## Examples

``` r
# Load example data and create ExpressionSet
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 


# Filter genes
filtered_eset <- filter_genes(eset, count_threshold = 10, percent_threshold = 20)
#> ✓ Gene filtering completed
#>   - Number of genes before filtering: 57010 
#>   - Number of genes after filtering: 19550 
#>   - Percentage of genes retained: 34.3 %
```
