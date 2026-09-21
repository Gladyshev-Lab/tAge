# Species supported by the clocks

The species accepted by
[`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
and
[`predict_tAge`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md),
with the maximum lifespan used to rescale chronological-age clocks and
the age units each species is reported in by default (months for
rodents, years for primates).

## Usage

``` r
tage_species()
```

## Value

A data frame with columns `species`, `max_lifespan_years` and
`default_units`.

## Examples

``` r
tage_species()
#>   species max_lifespan_years default_units
#> 1   mouse                4.0        months
#> 2     rat                3.8        months
#> 3   human              122.0         years
#> 4  monkey               39.0         years
```
