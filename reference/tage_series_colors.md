# Colour per group level

Assigns the categorical slots in the order the levels are given, with
the reference level in the neutral reference grey, so a colour stays
with its group whatever else is on the figure. This is how every
group-coloured tAge figure picks its colours.

## Usage

``` r
tage_series_colors(levels, reference = NULL, palette = NULL)
```

## Arguments

- levels:

  Group levels, in display order.

- reference:

  Level drawn in the reference grey. Default `NULL` (none).

- palette:

  Override: a named vector `c(level = colour)`, or an unnamed vector of
  colours in level order. Levels missing from a named vector are filled
  in automatically.

## Value

A named character vector of hex colours.

## Examples

``` r
tage_series_colors(c("WT", "KO", "HET"), reference = "WT")
#>        WT        KO       HET 
#> "#8a8987" "#2a78d6" "#eb6834" 
```
