# Heatmap of module-clock effects

Rows are co-expression modules, columns are strata (or the comparison
when there is no stratification), and each cell carries the effect with
a star when the adjusted p-value clears `sig_threshold`. Modules live on
different native scales, so the fill is the effect divided by a robust
scale of its column while the printed number stays in the original
units.

## Usage

``` r
tage_module_heatmap(
  data,
  module_columns,
  group_column,
  reference_group,
  compare_groups = NULL,
  covariates = NULL,
  split_by = NULL,
  variance_strata = c("all_data", "subset"),
  standardize = TRUE,
  p_adjust = "BH",
  p_adjust_scope = "across_columns",
  sig_threshold = 0.05,
  modules_version = c("4.6", "5.4", "human"),
  module_functions = NULL,
  color_scale = c("robust", "absolute"),
  robust_pct = 95,
  limit = NULL,
  annotate = TRUE,
  digits = 2,
  title = "Module clocks",
  subtitle = NA,
  caption = NULL,
  base_size = 11,
  cell_height = 0.3,
  cell_width = 1.15,
  width = NULL,
  height = NULL,
  stats = NULL
)
```

## Arguments

- data:

  Data frame with one column per module clock.

- module_columns:

  Module-clock columns.

- group_column:

  Column holding the experimental groups.

- reference_group:

  Reference level.

- compare_groups:

  Levels compared against the reference; each becomes a facet. `NULL`
  uses every other level.

- covariates, split_by, variance_strata, standardize, p_adjust,
  p_adjust_scope:

  Passed to
  [`tage_module_stats`](https://gladyshev-lab.github.io/tAge/reference/tage_module_stats.md).

- sig_threshold:

  Adjusted p-value below which a cell is starred.

- modules_version:

  Module set used for the row labels, one of `"4.6"`, `"5.4"` or
  `"human"`.

- module_functions:

  Named character vector overriding the bundled module-to-function
  annotation.

- color_scale:

  `"robust"` divides each column by its `robust_pct` percentile of
  `|effect|` so differently scaled columns share one colour bar;
  `"absolute"` uses the raw effect.

- robust_pct, limit:

  Percentile used by the robust scale, and the fill limit.

- annotate:

  Whether to print the effect in each cell.

- digits:

  Digits used for the printed effect.

- title, subtitle, caption:

  Figure text.

- base_size:

  Base font size.

- cell_height, cell_width:

  Inches per module row and per column, used to size the canvas so that
  row labels stay legible as the module set grows.

- width, height:

  Absolute figure size in inches, overriding `cell_height` and
  `cell_width`.

- stats:

  Pre-computed
  [`tage_module_stats`](https://gladyshev-lab.github.io/tAge/reference/tage_module_stats.md)
  output.

## Value

A `ggplot` object, with the statistics attached as the `"tage_stats"`
attribute and the intended size as `"tage_size"`; see
[`tage_save_plot`](https://gladyshev-lab.github.io/tAge/reference/tage_save_plot.md).

## Details

The statistics come from
[`tage_module_stats`](https://gladyshev-lab.github.io/tAge/reference/tage_module_stats.md).

## See also

[`tage_module_stats`](https://gladyshev-lab.github.io/tAge/reference/tage_module_stats.md),
[`tage_clock_forest`](https://gladyshev-lab.github.io/tAge/reference/tage_clock_forest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tage_module_heatmap(results, module_columns = modules,
                    group_column = "Genotype", reference_group = "WT",
                    split_by = "Tissue")
} # }
```
