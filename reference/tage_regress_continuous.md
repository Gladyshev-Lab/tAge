# Regress tAge on a continuous predictor

Slope of predicted tAge against a numeric predictor such as
chronological age, dose or time in culture. Elastic net clocks use
[`lm`](https://rdrr.io/r/stats/lm.html); Bayesian ridge clocks use a
REML meta-regression weighted by the per-sample prediction standard
deviation and report z-tests.

## Usage

``` r
tage_regress_continuous(
  data,
  value_columns,
  predictor,
  covariates = NULL,
  split_by = NULL,
  se_columns = NULL,
  p_adjust = TAGE_P_ADJUST_METHODS,
  p_adjust_scope = c("within_column", "across_columns", "global", "none"),
  conf_level = 0.95
)
```

## Arguments

- data:

  Data frame of per-sample predictions.

- value_columns:

  Character vector of tAge columns.

- predictor:

  Numeric column regressed against.

- covariates:

  Additional fixed-effect columns. Default `NULL`.

- split_by:

  Column to stratify on. Default `NULL`.

- se_columns:

  Per-sample standard-deviation columns for Bayesian ridge clocks; see
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md).

- p_adjust:

  Method passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- p_adjust_scope:

  Correction family; see
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md).
  Because there is one test per clock and stratum, `"across_columns"`
  corrects across clocks within a stratum.

- conf_level:

  Two-sided confidence level for `ci_low` / `ci_high`.

## Value

A data frame with `value_column`, `split`, `term`, `n`, `estimate`
(slope per unit of `predictor`), `se`, `ci_low`, `ci_high`, `statistic`,
`df`, `p_value`, `p_adjusted` and `label`.

## Examples

``` r
if (FALSE) { # \dontrun{
tage_regress_continuous(
  results,
  value_columns = "scaled_diff_EN_tAge",
  predictor     = "age_months",
  split_by      = "Tissue"
)
} # }
```
