# The tAge colour tokens

The colours every figure in the package is drawn with, for callers who
build their own plots and want them to match: categorical slots for
groups, the reference-group grey, one colour per clock outcome, the
diverging ramp for signed effects and the sequential ramp for
magnitudes, plus the text and grid tokens.

## Usage

``` r
tage_colors()
```

## Value

A named list of hex colours and colour vectors.

## Examples

``` r
tage_colors()$categorical
#> [1] "#2a78d6" "#eb6834" "#1baf7a" "#eda100" "#e87ba4" "#008300" "#4a3aa7"
#> [8] "#e34948"
```
