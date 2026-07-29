# Complete preprocessing pipeline for tAge analysis

This function performs a complete preprocessing pipeline for
transcriptomic age analysis, including gene filtering, normalization,
transformation, scaling, control subtraction, and gene ID conversion. It
returns multiple versions of the processed data suitable for different
analysis approaches.

## Usage

``` r
tAge_preprocessing(
  eset,
  species = "mouse",
  gene_mapping_type = "Gene.Symbol",
  verbose = TRUE,
  control_group_column = NULL,
  control_group_label = NULL,
  count_threshold = 10,
  percent_threshold = 20
)
```

## Arguments

- eset:

  An ExpressionSet object containing raw expression data.

- species:

  Character string specifying the species. Options: "mouse", "rat",
  "human", "monkey", "rhesus". Default is "mouse".

- gene_mapping_type:

  Gene mapping type. Options: "Gene.Symbol", "Ensembl"

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

- control_group_column:

  Character string specifying the column name in phenoData that contains
  control group labels. Default is NULL.

- control_group_label:

  Character string specifying the label for control samples. Default is
  NULL.

- count_threshold:

  Numeric threshold for minimum expression count in gene filtering.
  Default is 10.

- percent_threshold:

  Numeric threshold for minimum percentage of samples that must have
  expression above count_threshold. Default is 20.

## Value

A list containing six processed ExpressionSet objects:

- RLE_normalized:

  RLE-normalized data

- log_transformed:

  Log-transformed data

- scaled:

  Scaled data with gene ID conversion

- scaled_diff:

  Scaled data with control subtraction and gene ID conversion

- yugene:

  YuGene-normalized data with gene ID conversion

- yugene_diff:

  YuGene-normalized data with control subtraction and gene ID conversion

## Examples

``` r
# Load example data
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
eset <- make_ExpressionSet(expr_data, meta_data)
#> ✓ ExpressionSet created successfully
#>   - Number of genes: 57010 
#>   - Number of samples: 24 


# Run complete preprocessing pipeline
processed_data <- tAge_preprocessing(eset, species = "mouse")
#> ✓ Gene filtering completed
#>   - Number of genes before filtering: 57010 
#>   - Number of genes after filtering: 19550 
#>   - Percentage of genes retained: 34.3 %
#> Error in `sampleNames<-`(`*tmp*`, value = sampleNames(phenoData)): 'value' length (24) must equal sample number in AssayData (1)
```
