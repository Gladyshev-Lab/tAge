# Size a tAge figure was designed for

The figure functions grow the canvas with the number of clocks, modules
and panels they were given, and record the result. `tage_fig_size()`
reads it back so the same size can be reused for `ggsave()`, a knitr
chunk (`fig.width` / `fig.height`) or a manual tweak.

## Usage

``` r
tage_fig_size(p)
```

## Arguments

- p:

  A plot from
  [`tage_clock_forest`](https://gladyshev-lab.github.io/tAge/reference/tage_clock_forest.md)
  or
  [`tage_module_heatmap`](https://gladyshev-lab.github.io/tAge/reference/tage_module_heatmap.md).

## Value

Named numeric vector with `width` and `height` in inches, or `NULL` when
the plot carries no recorded size.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- tage_clock_forest(results, clocks_meta = clocks,
                       group_column = "Genotype", reference_group = "WT")
tage_fig_size(p)
} # }
```
