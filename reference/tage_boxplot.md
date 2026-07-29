# Box plot of tAge predictions with pairwise significance annotation

Draws a box plot with jittered points for one prediction column, split
by a grouping variable and optionally faceted by a subgroup. Pairwise
comparisons are annotated with brackets; comparisons involving groups
with too few observations, and — when `p_threshold` is set —
non-significant ones, are dropped before plotting so the panel stays
readable.

## Usage

``` r
tage_boxplot(
  data,
  x_var,
  y_var,
  subgroup_var = NULL,
  colors = NULL,
  point_size = 2,
  point_alpha = 0.7,
  box_width = 0.5,
  stat_method = "t.test",
  comparisons = NULL,
  p_label = "p.signif",
  p_threshold = NULL,
  min_group_n = 2,
  font_size = 12,
  theme_type = "classic",
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  legend_position = "right",
  y_center = NULL,
  y_min = NULL,
  y_max = NULL,
  facet_scales = "free_y",
  x_order = NULL,
  width = 10,
  height = 6
)
```

## Arguments

- data:

  Data frame of per-sample values, typically the prediction table
  returned by
  [`predict_tAge`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md).

- x_var:

  Column name defining the groups on the x axis.

- y_var:

  Column name holding the values to plot, e.g. a `*_tAge` column.

- subgroup_var:

  Optional column name used to facet the plot. Default `NULL`.

- colors:

  Optional named vector of fill colours, keyed by the levels of `x_var`.

- point_size, point_alpha:

  Size and opacity of the jittered points.

- box_width:

  Width of the boxes.

- stat_method:

  Test used for the pairwise comparisons, passed to
  `ggpubr::stat_compare_means`. Default `"t.test"`.

- comparisons:

  List of length-2 character vectors giving the pairs to test. Default
  `NULL` tests all pairs of `x_var` levels.

- p_label:

  Label style for the annotations, e.g. `"p.signif"` for stars or
  `"p.format"` for numeric p-values.

- p_threshold:

  If supplied, comparisons whose p-value is not below this threshold are
  removed from the plot. When faceting, a comparison is kept if it is
  significant in at least one facet.

- min_group_n:

  Minimum number of non-missing observations required in both groups for
  a comparison to be shown. When faceting, every facet must meet it.

- font_size:

  Base font size; also scales the annotation text.

- theme_type:

  ggplot2 theme to apply, e.g. `"classic"` or `"minimal"`.

- title, xlab, ylab:

  Plot title and axis labels.

- legend_position:

  Legend placement passed to
  [`ggplot2::theme`](https://ggplot2.tidyverse.org/reference/theme.html).

- y_center:

  If given, the y axis is made symmetric around this value — useful for
  relative predictions centred on zero.

- y_min, y_max:

  Explicit y-axis limits, overriding the automatic range.

- facet_scales:

  Scale behaviour across facets, passed to
  [`ggplot2::facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
  Default `"free_y"`.

- x_order:

  Character vector giving the order of the `x_var` levels.

- width, height:

  Plot size in inches, used to set the inline display size.

## Value

A `ggplot` object.

## See also

[`predict_tAge`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md)
for producing the input table.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "EN")
tage_boxplot(
  results,
  x_var = "Genotype",
  y_var = "yugene_diff_EN_tAge",
  subgroup_var = "Tissue",
  x_order = c("WT", "Klotho KO"),
  ylab = "Relative tAge, months"
)
} # }
```
