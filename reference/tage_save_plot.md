# Save a tAge figure at the size it was designed for

Thin wrapper over
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
that defaults `width` and `height` to the size recorded by the figure
function, so a figure with thirty modules is not squeezed into the same
canvas as one with three.

## Usage

``` r
tage_save_plot(p, filename, width = NULL, height = NULL, dpi = 300, ...)
```

## Arguments

- p:

  A plot from
  [`tage_clock_forest`](https://gladyshev-lab.github.io/tAge/reference/tage_clock_forest.md)
  or
  [`tage_module_heatmap`](https://gladyshev-lab.github.io/tAge/reference/tage_module_heatmap.md).

- filename:

  Output path; the extension picks the device.

- width, height:

  Size in inches. `NULL` (default) uses the recorded size, falling back
  to ggplot2's own default when there is none.

- dpi:

  Resolution for raster devices.

- ...:

  Passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Value

The path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- tage_clock_forest(results, clocks_meta = clocks,
                       group_column = "Genotype", reference_group = "WT")
tage_save_plot(p, "forest.png")                 # recorded size
tage_save_plot(p, "forest.pdf", width = 7)      # fixed width, recorded height
} # }
```
