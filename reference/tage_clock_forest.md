# Forest plot of clock effects

One row per clock: the effect between two groups, whiskers spanning the
confidence interval, filled when the adjusted p-value clears
`sig_threshold` and hollow otherwise. Panels are laid out as outcome
(rows) by stratum (columns). This is the figure to reach for when
several clocks are applied to one comparison – a heatmap hides the
uncertainty and a box plot per clock does not fit on a page.

## Usage

``` r
tage_clock_forest(
  data,
  value_columns = NULL,
  group_column,
  reference_group = NULL,
  compare_groups = NULL,
  clocks_meta = NULL,
  covariates = NULL,
  split_by = NULL,
  se_columns = NULL,
  variance_strata = c("subset", "all_data"),
  p_adjust = "BH",
  p_adjust_scope = "across_columns",
  conf_level = 0.95,
  sig_threshold = 0.05,
  sort_by_effect = TRUE,
  label_column = "name",
  outcome_column = "outcome",
  units = NULL,
  title = "Clock effects",
  subtitle = NA,
  caption = NULL,
  base_size = 11,
  row_height = 0.3,
  panel_width = 4.9,
  width = NULL,
  height = NULL,
  stats = NULL
)
```

## Arguments

- data:

  Data frame of per-sample predictions, e.g. from
  [`predict_tAge`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md).

- value_columns:

  Prediction columns to show. May be omitted when `clocks_meta` is
  given, in which case its `filename` column is used.

- group_column:

  Column holding the experimental groups.

- reference_group:

  Level everything is compared against.

- compare_groups:

  Levels compared with the reference. `NULL` uses every other level;
  more than one is drawn as separate coloured series.

- clocks_meta:

  Clock table from
  [`list_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md).
  Supplies the row labels and the outcome each clock predicts, which
  sets the panel rows, the colours and the x-axis units.

- covariates, split_by, se_columns, variance_strata, p_adjust,
  p_adjust_scope, conf_level:

  Passed to
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md).
  `p_adjust_scope` defaults to `"across_columns"`, i.e. correction
  across the clocks in one panel, which is what "filled = significant"
  implies.

- sig_threshold:

  Adjusted p-value below which a marker is filled.

- sort_by_effect:

  Order clocks by effect size within each panel.

- label_column, outcome_column:

  Columns of `clocks_meta` holding the display label and the outcome.

- units:

  Named character vector overriding the x-axis unit per outcome, e.g.
  `c(Chronological = "years")` for human data.

- title, subtitle, caption:

  Figure text. `subtitle` defaults to a description of the markers; pass
  `NULL` to drop it.

- base_size:

  Base font size.

- row_height, panel_width:

  Inches per clock row and per panel column, used to size the canvas so
  that rows stay legible as the clock set grows.

- width, height:

  Absolute figure size in inches, overriding `row_height` and
  `panel_width`.

- stats:

  Pre-computed
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md)
  output to plot instead of recomputing.

## Value

A `ggplot` object. The statistics are attached as the `"tage_stats"`
attribute and the intended size as `"tage_size"`; see
[`tage_save_plot`](https://gladyshev-lab.github.io/tAge/reference/tage_save_plot.md).

## Details

The statistics come from
[`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md),
so the figure and the numbers behind it cannot drift apart.

## See also

[`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md),
[`tage_module_heatmap`](https://gladyshev-lab.github.io/tAge/reference/tage_module_heatmap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
clocks <- list_clocks(type = "EN", tissue = "Multi-Tissue")
tage_clock_forest(results, clocks_meta = clocks,
                  group_column = "Genotype", reference_group = "WT",
                  split_by = "Tissue")
} # }
```
