# Module-clock effect sizes and p-values

Per-module statistics for the module-clock heatmaps. With
`standardize = TRUE` (the default) the two quantities come from
different models: the p-value is the estimated marginal-mean contrast of
the full model fitted on the stratum, while the effect size is the
coefficient of a separate two-group model fitted on values standardised
to unit variance, which makes modules comparable despite their different
native scales. With `standardize = FALSE` the estimate, standard error
and p-value all come from the same contrast.

## Usage

``` r
tage_module_stats(
  data,
  module_columns,
  group_column,
  reference_group,
  compare_groups = NULL,
  covariates = NULL,
  split_by = NULL,
  variance_strata = c("all_data", "subset"),
  standardize = TRUE,
  p_adjust = TAGE_P_ADJUST_METHODS,
  p_adjust_scope = c("across_columns", "within_column", "global", "none"),
  conf_level = 0.95
)
```

## Arguments

- data:

  Data frame with one column per module clock.

- module_columns:

  Character vector of module-clock columns.

- group_column:

  Column defining the experimental groups.

- reference_group:

  Reference level of `group_column`.

- compare_groups:

  Levels compared against the reference. Default `NULL` uses every other
  level.

- covariates:

  Covariate columns. Default `NULL`.

- split_by:

  Column to stratify on. Default `NULL`.

- variance_strata:

  `"all_data"` (default, matching the application) fits the p-value
  model on every group in the stratum; `"subset"` restricts it to the
  compared groups.

- standardize:

  Whether to report standardised effect sizes. Default `TRUE`.

- p_adjust:

  Method passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- p_adjust_scope:

  Correction family. Default `"across_columns"` corrects across modules
  within each comparison and stratum, which is the application's "per
  group" setting; `"global"` matches its "globally" setting.

- conf_level:

  Two-sided confidence level for `ci_low` / `ci_high`.

## Value

A data frame with `module`, `split`, `group1`, `group2`, `n`,
`estimate`, `se`, `ci_low`, `ci_high`, `statistic`, `df`, `p_value`,
`p_adjusted` and `label`.

## Examples

``` r
if (FALSE) { # \dontrun{
tage_module_stats(
  results,
  module_columns  = grep("^module_", names(results), value = TRUE),
  group_column    = "Genotype",
  reference_group = "WT",
  split_by        = "Tissue"
)
} # }
```
