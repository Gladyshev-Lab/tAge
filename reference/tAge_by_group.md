# Run tAge pipeline separately per group factor (e.g., "tissue")

Preprocesses every level of `split_by` on its own – gene filtering,
normalisation and reference centring all happen within the stratum, as
the clocks were trained and applied in the paper – and predicts on the
combined data. This is
[`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
with `split_by` followed by
[`predict_tAge`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md).

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
  verbose = TRUE,
  gene_mapping_type = "auto",
  return_std = identical(mode, "BR"),
  age_units = c("auto", "months", "years"),
  normalized_age = c("fraction", "percent")
)
```

## Arguments

- eset:

  An ExpressionSet, e.g. from pseudobulk aggregation.

- split_by:

  Character. Column in pData to split by (e.g., "tissue").

- model_paths:

  Named list of model paths.

- species:

  Species of the samples: "mouse", "rat", "human" or "monkey" (see
  [`tage_species`](https://gladyshev-lab.github.io/tAge/reference/tage_species.md)).
  Default is "mouse".

- mode:

  Character string specifying the model type. Must be either "EN" for
  Elastic Net or "BR" for Bayesian Ridge.

- control_group_column:

  Character string specifying the column name in phenoData that contains
  control group labels. Default is NULL.

- control_group_label:

  Character string specifying the label for control samples. Default is
  NULL. With `split_by`, the controls of each stratum are its reference;
  a stratum without controls is centred on all of its samples, with a
  warning.

- count_threshold:

  Numeric threshold for minimum expression count in gene filtering.
  Default is 10.

- percent_threshold:

  Numeric threshold for minimum percentage of samples that must have
  expression above count_threshold. Default is 20.

- min_samples:

  Integer. Strata with fewer samples are left out, with a warning.
  Default 5.

- verbose:

  Logical indicating whether to print progress messages. Default is
  TRUE.

- gene_mapping_type:

  Identifier type of the row names: "Ensembl", "Gene.Symbol", "Entrez"
  or "auto" (default, detected from the gene table; see
  [`map_genes`](https://gladyshev-lab.github.io/tAge/reference/map_genes.md)).

- return_std:

  Logical. Whether to also return the per-sample predictive standard
  deviation, which only Bayesian Ridge models provide. Defaults to
  `TRUE` for `mode = "BR"`. The standard deviations are what
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md)
  weights samples by, so keep them if you intend to run statistics on BR
  predictions.

- age_units:

  Units of chronological-age predictions: `"auto"` (default; months for
  rodents, years for primates), `"months"` or `"years"`.

- normalized_age:

  Scale of normalized-age clocks: `"fraction"` of the expected maximum
  lifespan (default) or `"percent"` (x100, as in the paper).

## Value

Data frame with predictions for all strata combined (the `split_by`
column is part of the sample metadata).
