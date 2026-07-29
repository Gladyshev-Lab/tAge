# Predict transcriptomic age for multiple processed ExpressionSet objects

This function predicts transcriptomic age for multiple processed
ExpressionSet objects using pre-trained models. It supports different
normalization methods and model types.

## Usage

``` r
predict_tAge(tAge_eset, model_paths, species, mode)
```

## Arguments

- tAge_eset:

  A named list of ExpressionSet objects, each representing a different
  normalization method (e.g., "scaled", "scaled_diff", "yugene",
  "yugene_diff").

- model_paths:

  A named list of model paths corresponding to each normalization
  method.

- species:

  Character string specifying the species for the models.

- mode:

  Character string specifying the model type. Must be either "EN" for
  Elastic Net or "BR" for Bayesian Ridge.

## Value

A data frame containing the predicted transcriptomic age results for all
provided ExpressionSet objects, with appropriately named columns.
