# Prediction tests — require Python modules and download clock models from
# Zenodo. They skip gracefully when either is unavailable (e.g. on CRAN or
# offline). Golden values were cross-checked bit-for-bit against an independent
# implementation of the same models on the shipped mouse example data.

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

  # Golden values (cross-checked against an independent implementation).
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

test_that("Bayesian ridge clocks return per-sample standard deviations", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("BR_Chronoage_Multispecies_Multitissue_yugenediff.pkl")
  skip_if(is.null(model), "BR chronological clock model could not be downloaded")

  eset <- .tage_example_eset()
  processed <- tAge_preprocessing(eset, species = "mouse",
                                  gene_mapping_type = "Ensembl", verbose = FALSE)
  res <- predict_tAge(processed, list(yugene_diff = model),
                      species = "mouse", mode = "BR")

  # return_std defaults to TRUE for BR, adding one _sd column per clock.
  expect_true("yugene_diff_BR_tAge_sd" %in% colnames(res))
  sds <- res$yugene_diff_BR_tAge_sd
  expect_true(all(is.finite(sds)))
  expect_true(all(sds > 0))

  # The SD must share the scale of the prediction. This clock is chronological,
  # so both are multiplied by the mouse maximum lifespan (48 months); an SD left
  # on the normalised 0-1 scale would be ~1/48 of the prediction's spread.
  expect_gt(stats::median(sds), 0.1)
  expect_lt(stats::median(sds), 48)

  # And they are usable as meta-regression weights.
  stats_tbl <- tage_compare_groups(
    res,
    value_columns   = "yugene_diff_BR_tAge",
    group_column    = "Genotype",
    reference_group = res$Genotype[1],
    se_columns      = "yugene_diff_BR_tAge_sd",
    p_adjust        = "none"
  )
  expect_gt(nrow(stats_tbl), 0)
  expect_true(all(is.finite(stats_tbl$p_value)))
  expect_true(all(is.infinite(stats_tbl$df)))
})

test_that("elastic net clocks refuse return_std instead of silently ignoring it", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
  skip_if(is.null(model), "mortality clock model could not be downloaded")

  eset <- .tage_example_eset()
  processed <- tAge_preprocessing(eset, species = "mouse",
                                  gene_mapping_type = "Ensembl", verbose = FALSE)

  expect_warning(
    res <- predict_tAge(processed, list(scaled_diff = model),
                        species = "mouse", mode = "EN", return_std = TRUE),
    "do not provide predictive standard deviations"
  )
  expect_false("scaled_diff_EN_tAge_sd" %in% colnames(res))
})

test_that("species comes from the preprocessed object, and units can be chosen", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl")
  skip_if(is.null(model), "chronological clock model could not be downloaded")

  processed <- tAge_preprocessing(.tage_example_eset(), species = "mouse", verbose = FALSE)

  implicit <- predict_tAge(processed, list(scaled_diff = model), mode = "EN")
  explicit <- predict_tAge(processed, list(scaled_diff = model), species = "mouse", mode = "EN")
  expect_equal(implicit$scaled_diff_EN_tAge, explicit$scaled_diff_EN_tAge)
  expect_equal(attr(implicit, "tage_units"), c(scaled_diff_EN_tAge = "months"))

  years <- predict_tAge(processed, list(scaled_diff = model), mode = "EN", age_units = "years")
  expect_equal(years$scaled_diff_EN_tAge * 12, implicit$scaled_diff_EN_tAge)
  expect_equal(attr(years, "tage_units"), c(scaled_diff_EN_tAge = "years"))

  expect_error(predict_tAge(processed, list(scaled_diff = model), species = "rhesus", mode = "EN"),
               "Unknown species")
})

test_that("normalized-age clocks can be reported in percent", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("EN_NormalizedAge_Multispecies_Multitissue_scaleddiff.pkl")
  skip_if(is.null(model), "normalized-age clock model could not be downloaded")

  processed <- tAge_preprocessing(.tage_example_eset(), species = "mouse", verbose = FALSE)
  frac <- predict_tAge(processed, list(scaled_diff = model), mode = "EN")
  pct  <- predict_tAge(processed, list(scaled_diff = model), mode = "EN", normalized_age = "percent")
  expect_equal(pct$scaled_diff_EN_tAge, frac$scaled_diff_EN_tAge * 100)
  expect_true(all(abs(frac$scaled_diff_EN_tAge) < 2))
  expect_equal(unname(attr(frac, "tage_units")), "fraction of maximum lifespan")
  expect_equal(unname(attr(pct, "tage_units")), "% of maximum lifespan")
})

test_that("tAge_by_group predicts per-tissue preprocessed data", {
  skip_on_cran()
  skip_if_offline()
  skip_if_no_python()
  model <- .tage_test_model("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl")
  skip_if(is.null(model), "mortality clock model could not be downloaded")

  res <- tAge_by_group(.tage_example_eset(), split_by = "Tissue",
                       model_paths = list(scaled_diff = model), species = "mouse", mode = "EN",
                       control_group_column = "Genotype", control_group_label = "WT",
                       verbose = FALSE)
  expect_equal(nrow(res), 24)
  expect_true(all(c("Tissue", "scaled_diff_EN_tAge") %in% names(res)))
  expect_true(all(is.finite(res$scaled_diff_EN_tAge)))
})
