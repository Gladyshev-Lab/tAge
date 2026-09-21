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

test_that("downloaded files are checked to be pickles", {
  good <- tempfile(fileext = ".pkl")
  writeBin(as.raw(c(0x80, 0x05, 0x95)), good)
  expect_true(.tage_check_download(good))
  expect_true(file.exists(good))

  html <- tempfile(fileext = ".pkl")
  writeLines("<!DOCTYPE html><html>Not found</html>", html)
  expect_error(.tage_check_download(html), "looks like an HTML page")
  expect_false(file.exists(html))

  empty <- tempfile(fileext = ".pkl")
  file.create(empty)
  expect_error(.tage_check_download(empty), "not a model file")
  expect_false(file.exists(empty))
})

test_that("a failed download leaves no partial file behind", {
  dest <- tempfile(fileext = ".pkl")
  expect_error(
    suppressWarnings(.tage_download_file("file:///nonexistent/path/model.pkl", dest, quiet = TRUE)),
    "failed"
  )
  expect_false(file.exists(dest))
})

test_that("download_clocks raises the timeout only while it runs", {
  old <- getOption("timeout")
  on.exit(options(timeout = old))
  options(timeout = 60)
  # Nothing to download: the file is already present, so no network is used.
  dir <- tempfile(); dir.create(dir)
  fn  <- "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"
  writeBin(as.raw(0x80), file.path(dir, fn))
  res <- download_clocks(fn, dest_dir = dir, quiet = TRUE, timeout = 7200)
  expect_equal(basename(res$path), fn)
  expect_equal(getOption("timeout"), 60)
})

test_that("clock outcomes are recognised from registry and internal names", {
  co <- tAge:::.clock_outcome
  expect_equal(co("EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"), "Mortality")
  expect_equal(co("BR_NormalizedAge_Rodents_Liver_yugenediff.pkl"), "Normalized age")
  expect_equal(co("/x/EN_Hazard_All_All_All_scaleddiff.pkl"), "Mortality")
  expect_equal(co("EN_Relage999_Mouse_All_All_yugenediff.pkl"), "Normalized age")
  expect_equal(co("EN_Chronoage_blue_Rodents_All_WT_scaleddiff.pkl"), "Chronological")
  expect_equal(co("EN_Lifespan_Rodents_All_All_scaleddiff.pkl"), "Lifespan")
  expect_true(is.na(co("Ovarian_clock_scaleddiff_EN.pkl")))
})

test_that("output scaling follows outcome, species and units", {
  sc <- tAge:::.tage_output_scaling
  chrono <- "EN_Chronoage_Mouse_Multitissue_scaleddiff.pkl"
  expect_equal(sc(chrono, "mouse", "auto", "fraction"), list(factor = 48, units = "months"))
  expect_equal(sc(chrono, "mouse", "years", "fraction"), list(factor = 4, units = "years"))
  expect_equal(sc(chrono, "human", "auto", "fraction"), list(factor = 122.5, units = "years"))
  expect_equal(sc(chrono, "human", "months", "fraction")$factor, 122.5 * 12)
  na <- "EN_NormalizedAge_Mouse_Multitissue_scaleddiff.pkl"
  expect_equal(sc(na, "mouse", "auto", "fraction")$factor, 1)
  expect_equal(sc(na, "mouse", "auto", "percent"), list(factor = 100, units = "% of maximum lifespan"))
  expect_equal(sc("EN_Mortality_Mouse_Multitissue_scaleddiff.pkl", "rat", "auto", "fraction"),
               list(factor = 1, units = "log10 hazard ratio"))
  expect_warning(res <- sc("Ovarian_clock_scaleddiff_EN.pkl", "mouse", "auto", "fraction"),
                 "Could not tell")
  expect_equal(res$factor, 1)
  expect_error(sc(chrono, "rhesus", "auto", "fraction"), "Unknown species")
})

test_that("tage_species lists the four species with lifespans and units", {
  sp <- tage_species()
  expect_setequal(sp$species, c("mouse", "rat", "human", "monkey"))
  expect_equal(sp$default_units[sp$species == "mouse"], "months")
  expect_equal(sp$default_units[sp$species == "human"], "years")
})
