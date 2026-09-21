# Significance stars for tAge statistics

Maps p-values onto the label set of the published figures: `***` below
0.001, `**` below 0.01, `*` below 0.05 and `^` below 0.1.

## Usage

``` r
tage_significance_stars(p)
```

## Arguments

- p:

  Numeric vector of p-values.

## Value

Character vector of labels; `""` where nothing is significant and `NA`
stays `NA`.

## Examples

``` r
tage_significance_stars(c(1e-4, 0.02, 0.08, 0.5))
#> [1] "***" "*"   "^"   ""   
```
