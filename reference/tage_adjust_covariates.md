# Covariate-adjusted tAge values for plotting

Removes the fitted covariate effects from a tAge column while keeping
the group effect, and adds the mean covariate effect back so the
adjusted values stay on the original scale.

## Usage

``` r
tage_adjust_covariates(
  data,
  value_column,
  covariates,
  split_by = NULL,
  se_column = NULL,
  group_column = NULL
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

  Optional stratifying column. With `group_column`, one model per
  stratum, as in
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md);
  without it, a fixed effect for elastic net clocks and one model per
  stratum for Bayesian ridge clocks.

- se_column:

  Per-sample standard deviations for a Bayesian ridge clock. Default
  `NULL` uses [`lm`](https://rdrr.io/r/stats/lm.html).

- group_column:

  Column holding the experimental groups, fitted together with the
  covariates. Default `NULL` keeps the covariate-only model.

## Value

Numeric vector of adjusted values, in the row order of `data`; rows with
missing values in the model return `NA`.

## Details

With `group_column` the covariate effects are those of the model the
statistics fit – `value ~ group + covariates`, one model per stratum,
`lm` or the weighted meta-regression – so the difference between group
means of the adjusted values is the tested estimate (exactly for elastic
net clocks with `variance_strata = "all_data"`). Without it the
covariate model leaves the group out; when the covariates are unevenly
distributed across groups (more females among the treated, say) that
model attributes part of the group effect to the covariates, and the
adjusted values no longer show what was tested.

## Examples

``` r
if (FALSE) { # \dontrun{
results$adjusted <- tage_adjust_covariates(
  results, "yugene_diff_EN_tAge", covariates = "Sex", split_by = "Tissue",
  group_column = "Genotype"
)
} # }
```
