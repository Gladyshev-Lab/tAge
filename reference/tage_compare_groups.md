# Compare tAge between experimental groups

Estimated marginal-mean contrasts of predicted tAge between a reference
group and one or more comparison groups, as used for the clock analyses
in the paper.

## Usage

``` r
tage_compare_groups(
  data,
  value_columns,
  group_column,
  reference_group,
  compare_groups = NULL,
  covariates = NULL,
  split_by = NULL,
  se_columns = NULL,
  method = c("trt.vs.ctrl", "pairwise"),
  variance_strata = c("subset", "all_data"),
  p_adjust = TAGE_P_ADJUST_METHODS,
  p_adjust_scope = c("within_column", "across_columns", "global", "none"),
  conf_level = 0.95
)
```

## Arguments

- data:

  Data frame of per-sample predictions, typically the table returned by
  [`predict_tAge`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md).

- value_columns:

  Character vector of columns holding tAge predictions. Several columns
  are accepted so that `p_adjust_scope = "across_columns"` can correct
  across clocks or modules.

- group_column:

  Column defining the experimental groups.

- reference_group:

  Level of `group_column` used as the reference.

- compare_groups:

  Levels compared against the reference. Default `NULL` uses every other
  level present.

- covariates:

  Character vector of covariate columns entering the model as fixed
  effects. Default `NULL`.

- split_by:

  Column to stratify on, fitting a separate model per level (e.g.
  `"Tissue"`). Default `NULL`.

- se_columns:

  Columns holding per-sample prediction standard deviations from
  Bayesian ridge clocks. `NULL` (default) fits ordinary linear models.
  May be a single name, a vector parallel to `value_columns`, or a
  vector named by value column; `NA` entries fall back to `lm`.

- method:

  `"trt.vs.ctrl"` (default) contrasts every comparison group against the
  reference; `"pairwise"` contrasts all pairs.

- variance_strata:

  `"subset"` (default) estimates the residual variance from the compared
  groups only; `"all_data"` fits the model on every group present in the
  stratum and reports the requested contrasts from it.

- p_adjust:

  Method passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html);
  `"BH"` by default.

- p_adjust_scope:

  Family over which p-values are corrected. `"within_column"` (default)
  corrects across all comparisons and strata of one clock.
  `"across_columns"` corrects across clocks within each comparison and
  stratum (the family used for module clocks). `"global"` corrects
  everything together, `"none"` disables it.

- conf_level:

  Two-sided confidence level for `ci_low` / `ci_high`. The critical
  value follows the test: normal for the meta-regression, t otherwise.

## Value

A data frame with one row per contrast and the columns `value_column`,
`split`, `group1` (reference), `group2`, `n`, `estimate` (always
`group2 - group1`), `se`, `statistic`, `df`, `p_value`, `p_adjusted` and
`label`.

## Details

For elastic net clocks the model is `value ~ group + covariates` fitted
with [`lm`](https://rdrr.io/r/stats/lm.html). For Bayesian ridge clocks
– signalled by supplying `se_columns` – the model is a REML
meta-regression
([`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html))
that weights each sample by its own prediction uncertainty, and the
contrasts are z-tests rather than t-tests. In both cases contrasts come
from
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
with no built-in multiplicity adjustment; correction is applied
afterwards over the scope given by `p_adjust_scope`.

## See also

[`tage_regress_continuous`](https://gladyshev-lab.github.io/tAge/reference/tage_regress_continuous.md)
for numeric predictors,
[`tage_module_stats`](https://gladyshev-lab.github.io/tAge/reference/tage_module_stats.md)
for module clocks, and
[`tage_adjust_covariates`](https://gladyshev-lab.github.io/tAge/reference/tage_adjust_covariates.md)
for the matching plotting values.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "EN")
tage_compare_groups(
  results,
  value_columns   = "yugene_diff_EN_tAge",
  group_column    = "Genotype",
  reference_group = "WT",
  split_by        = "Tissue"
)
} # }
```
