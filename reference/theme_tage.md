# The tAge ggplot2 theme

Recessive axes (left and bottom only), a light grid on the value axis, a
left-aligned title with the subtitle in secondary ink and the caption in
muted ink, facet strips as plain semibold text, legend on the right.
Every tAge figure uses it; add it to your own plots to match.

## Usage

``` r
theme_tage(base_size = 10, grid = c("y", "x", "both", "none"))
```

## Arguments

- base_size:

  Base font size in points.

- grid:

  `"y"` (default), `"x"`, `"both"` or `"none"`: which axis gets the
  light grid.

## Value

A ggplot2 theme.

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
  ggplot2::geom_point() +
  theme_tage()
} # }
```
