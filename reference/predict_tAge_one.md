# Predict transcriptomic age with one pre-trained model

Applies one pre-trained Elastic Net (EN) or Bayesian Ridge (BR) clock to
a preprocessed ExpressionSet through the Python model (via reticulate).

## Usage

``` r
predict_tAge_one(
  eset,
  model_path,
  species = NULL,
  mode,
  return_std = identical(mode, "BR"),
  age_units = c("auto", "months", "years"),
  normalized_age = c("fraction", "percent")
)
```

## Arguments

- eset:

  An ExpressionSet object containing processed expression data (one
  element of the list returned by
  [`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)).

- model_path:

  Character string specifying the path to the pre-trained model file.

- species:

  Species of the *samples*: `"mouse"`, `"rat"`, `"human"` or `"monkey"`
  (see
  [`tage_species`](https://gladyshev-lab.github.io/tAge/reference/tage_species.md)).
  Its only effect is the rescaling of chronological-age clocks to age
  units by the species maximum lifespan; mortality and normalized-age
  clocks ignore it. It is *not* the species group the model was trained
  on ("Mouse" / "Rodents" / "Multispecies" in
  [`list_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md)).
  Default `NULL` takes the species recorded by
  [`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md).

- mode:

  Character string specifying the model type. Must be either "EN" for
  Elastic Net or "BR" for Bayesian Ridge.

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

A data frame containing the predicted transcriptomic age results with
sample information and predicted ages, plus a `BR_tAge_std` column when
`return_std` is `TRUE`. The attribute `"tage_units"` names the unit of
the prediction column.
