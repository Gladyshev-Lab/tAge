# Run tAge pipeline separately per group factor (e.g., "tissue")

Run tAge pipeline separately per group factor (e.g., "tissue")

## Usage

``` r
tAge_by_group(
  eset,
  split_by,
  model_paths,
  species = "mouse",
  mode = "EN",
  control_group_column = NULL,
  control_group_label = NULL,
  count_threshold = 10,
  percent_threshold = 20,
  min_samples = 5,
  verbose = TRUE
)
```

## Arguments

- eset:

  An ExpressionSet from pseudobulk aggregation.

- split_by:

  Character. Column in pData to split by (e.g., "tissue").

- model_paths:

  Named list of model paths.

- species:

  Character. Species name.

- mode:

  Character. "EN" or "BR".

- control_group_column:

  Character or NULL. Column for control subtraction.

- control_group_label:

  Character or NULL. Label for control group.

- count_threshold:

  Integer. Gene filtering threshold.

- percent_threshold:

  Numeric. Gene filtering percent.

- min_samples:

  Integer. Minimum samples per group to run tAge. Default 5.

- verbose:

  Logical.

## Value

Data frame with predictions from all groups combined.
