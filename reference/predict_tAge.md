# Predict transcriptomic age for multiple processed ExpressionSet objects

Applies one model per normalisation (`scaled_diff`, `yugene_diff`, ...)
and returns the sample metadata with one prediction column per
normalisation, named `<normalisation>_<mode>_tAge`.

## Usage

``` r
predict_tAge(
  tAge_eset,
  model_paths,
  species = NULL,
  mode,
  return_std = identical(mode, "BR"),
  age_units = c("auto", "months", "years"),
  normalized_age = c("fraction", "percent")
)
```

## Arguments

- tAge_eset:

  A named list of ExpressionSet objects, each representing a different
  normalization method (e.g., "scaled", "scaled_diff", "yugene",
  "yugene_diff"), as returned by
  [`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md).

- model_paths:

  A named list of model paths corresponding to each normalization
  method.

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

  Logical. Whether to keep the per-sample predictive standard deviation
  of Bayesian Ridge clocks. Defaults to `TRUE` for `mode = "BR"`, adding
  one `<normalisation>_BR_tAge_sd` column per clock. Pass these to the
  `se_columns` argument of
  [`tage_compare_groups`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md)
  for the Bayesian ridge statistics.

- age_units:

  Units of chronological-age predictions: `"auto"` (default; months for
  rodents, years for primates), `"months"` or `"years"`.

- normalized_age:

  Scale of normalized-age clocks: `"fraction"` of the expected maximum
  lifespan (default) or `"percent"` (x100, as in the paper).

## Value

A data frame containing the predicted transcriptomic age results for all
provided ExpressionSet objects, with appropriately named columns. The
attribute `"tage_units"` is a named character vector giving the unit of
every prediction column (e.g. `"months"`, `"log10 hazard ratio"`).
