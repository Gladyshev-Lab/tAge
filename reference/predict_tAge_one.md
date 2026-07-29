# Predict transcriptomic age using pre-trained models

This function predicts transcriptomic age using pre-trained Elastic Net
(EN) or Bayesian Ridge (BR) models. It interfaces with Python models
through the reticulate package to perform the predictions.

## Usage

``` r
predict_tAge_one(eset, model_path, species, mode)
```

## Arguments

- eset:

  An ExpressionSet object containing processed expression data.

- model_path:

  Character string specifying the path to the pre-trained model file.
  The file must exist and be compatible with the specified mode.

- species:

  Character string specifying the species for the model. This should
  match the species the model was trained on.

- mode:

  Character string specifying the model type. Must be either "EN" for
  Elastic Net or "BR" for Bayesian Ridge.

## Value

A data frame containing the predicted transcriptomic age results with
sample information and predicted ages.
