# Clock registry tests — no network required.

test_that("list_clocks returns the registry and filters correctly", {
  all_clocks <- list_clocks()
  expect_s3_class(all_clocks, "data.frame")
  expect_true(all(c("filename", "type", "outcome", "species", "tissue",
                    "scaling", "lifespan_scaled") %in% colnames(all_clocks)))
  expect_gt(nrow(all_clocks), 0)

  en_mort <- list_clocks(type = "EN", outcome = "Mortality")
  expect_true(all(en_mort$type == "EN"))
  expect_true(all(en_mort$outcome == "Mortality"))
  expect_lt(nrow(en_mort), nrow(all_clocks))
})

test_that("only chronological clocks are lifespan-scaled", {
  reg <- list_clocks()
  expect_true(all(reg$lifespan_scaled[reg$outcome == "Chronological"]))
  expect_false(any(reg$lifespan_scaled[reg$outcome != "Chronological"]))
})

test_that(".clock_lifespan_scaled resolves clock type from name", {
  expect_true(tAge:::.clock_lifespan_scaled(
    "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl"))
  expect_false(tAge:::.clock_lifespan_scaled(
    "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"))
  expect_false(tAge:::.clock_lifespan_scaled(
    "EN_NormalizedAge_Multispecies_Multitissue_scaleddiff.pkl"))
  # Unknown name -> NA (Python side falls back to its own heuristic)
  expect_true(is.na(tAge:::.clock_lifespan_scaled("some_unknown_model.pkl")))
})
