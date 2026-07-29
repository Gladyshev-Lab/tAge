# Load example metadata

Loads the example metadata included with the tAge package.

## Usage

``` r
load_example_metadata()
```

## Value

A data frame with sample metadata.

## Examples

``` r
meta_data <- load_example_metadata()
head(meta_data)
#>         Mouse.ID Genotype  Sex Tissue
#> Klo93K        93       WT Male Kidney
#> Klo94K        94       WT Male Kidney
#> Klo95K        95       WT Male Kidney
#> Klo96K        96       WT Male Kidney
#> Klo99K        99       WT Male Kidney
#> Klo100K      100       WT Male Kidney
```
