# Predict transcriptomic age with several clocks

Applies clocks to the preprocessed representations returned by
[`tAge_preprocessing`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
and returns the sample metadata with one prediction column per clock.

## Usage

``` r
predict_tAge(
  tAge_eset,
  model_paths,
  species = NULL,
  mode = NULL,
  return_std = NULL,
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

  A clock table with a `path` column, or a named list of model paths per
  representation; see Details.

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

  `"EN"` or `"BR"`. Default `NULL` takes each model's type from the
  clock table, the registry or its file name.

- return_std:

  Logical. Whether to keep the per-sample predictive standard deviation
  of Bayesian Ridge clocks. Default `NULL` keeps it for every Bayesian
  ridge clock, as a `<column>_sd` column. Pass these to the `se_columns`
  argument of
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

## Details

`model_paths` is either

- a clock table from
  [`list_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md)
  or
  [`list_module_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_module_clocks.md)
  with a `path` column (e.g. from
  [`download_clocks`](https://gladyshev-lab.github.io/tAge/reference/download_clocks.md)):
  each row is applied to the representation its `scaling` names
  (`Scaled` -\> `scaled_diff`, `YuGene` -\> `yugene_diff`), and the
  column is named after the model file, as in the Python package; or

- a named list, representation -\> model path(s). With one path per
  representation the column is `<representation>_<mode>_tAge`; with
  several, each column is named after its model file.

## Examples

``` r
if (FALSE) { # \dontrun{
clocks <- download_clocks(list_clocks(type = "EN", species = "Rodents",
                                      tissue = "Multi-Tissue"))
res <- predict_tAge(tAge_eset, clocks)       # six columns, named by model file
attr(res, "tage_units")
} # }
```
