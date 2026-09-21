# Diverging fill scale for signed effects

Blue for negative, neutral grey at zero, red for positive - the scale
every tAge heatmap uses. The limits are made symmetric around zero.

## Usage

``` r
scale_fill_tage_diverging(limit, name = ggplot2::waiver(), ...)
```

## Arguments

- limit:

  Absolute limit of the scale; clamp values to it beforehand.

- name:

  Legend title.

- ...:

  Passed to
  [`ggplot2::scale_fill_gradientn()`](https://ggplot2.tidyverse.org/reference/scale_gradient.html).

## Value

A ggplot2 scale.
