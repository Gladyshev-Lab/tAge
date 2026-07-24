# Prediction tests — require Python modules and download clock models from
# Zenodo. They skip gracefully when either is unavailable (e.g. on CRAN or
# offline). Golden values come from a bit-for-bit head-to-head against the
# TACO reference application on the shipped mouse example data.

test_that("mortality clock returns log10(HR), not lifespan-scaled values", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
  skip_if(is.null(model), "mortality clock model could not be downloaded")

  eset <- .tage_example_eset()
  processed <- tAge_preprocessing(eset, species = "mouse",
                                  gene_mapping_type = "Ensembl", verbose = FALSE)
  res <- predict_tAge(processed, list(scaled_diff = model),
                      species = "mouse", mode = "EN")
  vals <- res$scaled_diff_EN_tAge

  # Regression guard for the "mortality x maxlifespan" bug: log10(HR) values are
  # small; multiplying by 48 (mouse) would push them well outside this range.
  expect_true(all(abs(vals) < 5))

  # Golden values (match TACO exactly).
  expect_equal(head(vals, 3),
               c(-0.593993, -0.668232, -0.530628),
               tolerance = 1e-3)
})

test_that("chronological clock is rescaled to age units", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl")
  skip_if(is.null(model), "chronological clock model could not be downloaded")

  eset <- .tage_example_eset()
  processed <- tAge_preprocessing(eset, species = "mouse",
                                  gene_mapping_type = "Ensembl", verbose = FALSE)
  res <- predict_tAge(processed, list(scaled_diff = model),
                      species = "mouse", mode = "EN")
  vals <- res$scaled_diff_EN_tAge

  # Chronological output is normalised age x max lifespan (48 months for mouse),
  # so magnitudes are much larger than the raw normalised prediction.
  expect_equal(head(vals, 3),
               c(-9.7581, -7.7424, -6.3348),
               tolerance = 1e-2)
})
