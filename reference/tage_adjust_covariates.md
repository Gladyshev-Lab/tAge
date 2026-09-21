# Covariate-adjusted tAge values for plotting

Removes the fitted covariate effects from a tAge column while keeping
the group effect, so that box plots show what the model tested: the
covariate model is fitted without the grouping variable, and the mean
covariate effect is added back so the adjusted values stay on the
original scale.

## Usage

``` r
tage_adjust_covariates(
  data,
  value_column,
  covariates,
  split_by = NULL,
  se_column = NULL
)
```

## Arguments

- data:

  Data frame of per-sample predictions.

- value_column:

  Column to adjust.

- covariates:

  Covariate columns to regress out.

- split_by:

  Optional stratifying column, included as a fixed effect for elastic
  net clocks and used to fit one model per stratum for Bayesian ridge
  clocks.

- se_column:

  Per-sample standard deviations for a Bayesian ridge clock. Default
  `NULL` uses [`lm`](https://rdrr.io/r/stats/lm.html).

## Value

Numeric vector of adjusted values, in the row order of `data`; rows with
missing values in the model return `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
results$adjusted <- tage_adjust_covariates(
  results, "yugene_diff_EN_tAge", covariates = "Sex", split_by = "Tissue"
)
} # }
```
