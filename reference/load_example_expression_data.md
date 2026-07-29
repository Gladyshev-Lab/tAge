# Load example expression data

Loads the example expression data included with the tAge package.

## Usage

``` r
load_example_expression_data()
```

## Value

A data frame with expression data where rows are genes and columns are
samples.

## Examples

``` r
expr_data <- load_example_expression_data()
head(expr_data)
#>                    Klo93K Klo94K Klo95K Klo96K Klo99K Klo100K Klo113K Klo114K
#> ENSMUSG00000104478      0      0      0      0      0       0       0       0
#> ENSMUSG00000104385      0      0      0      0      0       0       0       0
#> ENSMUSG00000086053      0      0      0      0      0       0       0       0
#> ENSMUSG00000101231      0      0      0      0      0       0       0       0
#> ENSMUSG00000102135      0      2      6      4      2       0       2       6
#> ENSMUSG00000103282      0      0      0      0      0       0       0       0
#>                    Klo115K Klo119K Klo120K Klo122K RNA_93M RNA_94M RNA_95M
#> ENSMUSG00000104478       0       0       0       0       0       0       0
#> ENSMUSG00000104385       0       0       0       0       0       0       0
#> ENSMUSG00000086053       0       0       2       0       0       0      10
#> ENSMUSG00000101231       0       0       0       0       0       0       0
#> ENSMUSG00000102135       4       4       6       2       2       4       2
#> ENSMUSG00000103282       0       0       0       0       0       0       0
#>                    RNA_96M RNA_99M RNA_100M RNA_113M RNA_114M RNA_115M RNA_119M
#> ENSMUSG00000104478       0       0        0        0        0        2        4
#> ENSMUSG00000104385       0       0        0        0        0        0        0
#> ENSMUSG00000086053       1      18        0        3        2        0        7
#> ENSMUSG00000101231       0       0        0        0        0        0        0
#> ENSMUSG00000102135       4       0        0        2        6        2        2
#> ENSMUSG00000103282       0       0        0        0        0        0        0
#>                    RNA_120M RNA_122M
#> ENSMUSG00000104478        0        0
#> ENSMUSG00000104385        0        0
#> ENSMUSG00000086053        4        0
#> ENSMUSG00000101231        0        0
#> ENSMUSG00000102135        4        2
#> ENSMUSG00000103282        0        0
```
