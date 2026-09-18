# Figure tests. These check the contract -- what is drawn, from which
# statistics, with which labels -- rather than pixels, so they stay useful when
# the styling changes.

.figure_fixture <- function(seed = 7, n = 96) {
  set.seed(seed)
  d <- data.frame(
    Genotype = factor(rep(c("WT", "KO", "HET"), each = n / 3),
                      levels = c("WT", "KO", "HET")),
    Tissue   = rep(c("Kidney", "Muscle"), times = n / 2),
    Sex      = sample(c("M", "F"), n, replace = TRUE)
  )
  meta <- expand.grid(scaling = c("Scaled", "YuGene"),
                      species = c("Mouse", "Rodent"),
                      outcome = c("Chronological", "Mortality"),
                      stringsAsFactors = FALSE)
  meta$filename <- paste(meta$outcome, meta$species, meta$scaling, sep = "_")
  meta$name <- paste(meta$species, "· EN", meta$scaling)
  for (i in seq_len(nrow(meta))) {
    eff <- if (meta$outcome[i] == "Chronological") 3 else 0.4
    d[[meta$filename[i]]] <- eff * (d$Genotype == "KO") + stats::rnorm(n, sd = eff)
  }
  modules <- c("black", "blue", "brown", "cyan")
  for (m in modules) d[[m]] <- 0.5 * (d$Genotype == "KO") + stats::rnorm(n, sd = 0.6)
  list(data = d, meta = meta, modules = modules)
}

# Rendering is what catches most mistakes (missing aesthetics, bad factors),
# and it is cheap on these sizes.
.render <- function(p) {
  skip_if_not_installed("ggplot2")
  f <- tempfile(fileext = ".png")
  suppressMessages(ggplot2::ggsave(f, p, width = 9, height = 6, dpi = 72))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  unlink(f)
}

test_that("the forest plot carries the statistics it drew", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()

  p <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                         reference_group = "WT", compare_groups = "KO")
  expect_s3_class(p, "ggplot")

  st <- attr(p, "tage_stats")
  expect_true(is.data.frame(st))
  expect_setequal(st$value_column, fx$meta$filename)
  expect_true(all(c("ci_low", "ci_high", "p_adjusted") %in% names(st)))
  # The interval must bracket the estimate it describes.
  expect_true(all(st$ci_low <= st$estimate))
  expect_true(all(st$estimate <= st$ci_high))
  .render(p)
})

test_that("the forest plot can reuse a precomputed table", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()
  pre <- tage_compare_groups(fx$data, fx$meta$filename, "Genotype", "WT",
                             compare_groups = "KO", p_adjust = "none")
  p <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                         stats = pre)
  expect_equal(attr(p, "tage_stats")$estimate, pre$estimate)
  .render(p)
})

test_that("the forest plot keeps the significance stars", {
  # tage_compare_groups() returns a `label` column of stars, and the clock
  # display label must not overwrite it.
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()
  p <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                         reference_group = "WT", compare_groups = "KO")
  st <- attr(p, "tage_stats")
  expect_true(all(st$label %in% c("***", "**", "*", "^", "")))
  expect_false("label.x" %in% names(st))
})

test_that("the forest plot handles strata, covariates and several groups", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()

  p_split <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                               reference_group = "WT", compare_groups = "KO",
                               split_by = "Tissue")
  expect_setequal(attr(p_split, "tage_stats")$split, c("Kidney", "Muscle"))
  .render(p_split)

  p_cov <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                             reference_group = "WT", compare_groups = "KO",
                             covariates = "Sex")
  .render(p_cov)

  p_multi <- tage_clock_forest(fx$data, clocks_meta = fx$meta,
                               group_column = "Genotype", reference_group = "WT")
  expect_setequal(attr(p_multi, "tage_stats")$group2, c("KO", "HET"))
  .render(p_multi)
})

test_that("the forest plot works without clock metadata", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()
  p <- tage_clock_forest(fx$data, value_columns = fx$meta$filename,
                         group_column = "Genotype", reference_group = "WT",
                         compare_groups = "KO")
  .render(p)
})

test_that("the module heatmap carries the statistics it drew", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()

  p <- tage_module_heatmap(fx$data, fx$modules, "Genotype", "WT",
                           compare_groups = "KO", split_by = "Tissue")
  expect_s3_class(p, "ggplot")
  st <- attr(p, "tage_stats")
  expect_setequal(st$module, fx$modules)
  expect_setequal(st$split, c("Kidney", "Muscle"))
  .render(p)
})

test_that("module standardisation makes differently scaled modules comparable", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()
  d <- fx$data
  d$blue <- d$blue * 100                     # same signal, different units

  std <- attr(tage_module_heatmap(d, fx$modules, "Genotype", "WT",
                                  compare_groups = "KO", standardize = TRUE),
              "tage_stats")
  raw <- attr(tage_module_heatmap(d, fx$modules, "Genotype", "WT",
                                  compare_groups = "KO", standardize = FALSE),
              "tage_stats")
  blue_std <- abs(std$estimate[std$module == "blue"])
  blue_raw <- abs(raw$estimate[raw$module == "blue"])
  expect_lt(blue_std, blue_raw / 10)
})

test_that("figures reject columns that are not in the data", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()
  expect_error(
    tage_clock_forest(fx$data, value_columns = "nope", group_column = "Genotype",
                      reference_group = "WT"),
    "None of the requested"
  )
  expect_error(
    tage_module_heatmap(fx$data, "nope", "Genotype", "WT"),
    "None of the requested"
  )
})

test_that("module functions are looked up from the bundled annotation", {
  funcs <- load_module_functions("4.6")
  skip_if(length(funcs) == 0, "module annotation not installed")
  expect_true("blue" %in% names(funcs))
  expect_match(unname(funcs[["blue"]]), "Respiration|Mitochondrial|Muscle")

  expect_true(length(load_module_functions("5.4")) > 0)
  expect_true(length(load_module_functions("human")) > 0)
})

test_that("figures record the size they were designed for", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()

  p <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                         reference_group = "WT", compare_groups = "KO",
                         split_by = "Tissue")
  size <- tage_fig_size(p)
  expect_named(size, c("width", "height"))
  expect_true(all(size > 0))

  # Two strata at the default panel width.
  expect_equal(unname(size[["width"]]), 2 * 4.9)

  # Explicit size wins over the per-panel sizing.
  fixed <- tage_clock_forest(fx$data, clocks_meta = fx$meta, group_column = "Genotype",
                             reference_group = "WT", compare_groups = "KO",
                             width = 7, height = 4)
  expect_equal(unname(tage_fig_size(fixed)), c(7, 4))
})

test_that("the auto height grows with the number of rows", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()

  few <- tage_clock_forest(fx$data, clocks_meta = fx$meta[1:2, ],
                           group_column = "Genotype", reference_group = "WT",
                           compare_groups = "KO")
  many <- tage_clock_forest(fx$data, clocks_meta = fx$meta,
                            group_column = "Genotype", reference_group = "WT",
                            compare_groups = "KO")
  expect_gte(tage_fig_size(many)[["height"]], tage_fig_size(few)[["height"]])

  small <- tage_module_heatmap(fx$data, fx$modules[1:2], "Genotype", "WT",
                               compare_groups = "KO")
  big <- tage_module_heatmap(fx$data, fx$modules, "Genotype", "WT",
                             compare_groups = "KO")
  expect_gt(tage_fig_size(big)[["height"]], tage_fig_size(small)[["height"]])
})

test_that("tage_save_plot uses the recorded size and honours overrides", {
  skip_if_not_installed("ggplot2")
  fx <- .figure_fixture()
  p <- tage_module_heatmap(fx$data, fx$modules, "Genotype", "WT",
                           compare_groups = "KO")
  size <- tage_fig_size(p)

  f <- tempfile(fileext = ".png")
  suppressMessages(tage_save_plot(p, f, dpi = 72))
  expect_true(file.exists(f))
  dim_recorded <- dim(png::readPNG(f))
  expect_equal(dim_recorded[2] / 72, unname(size[["width"]]), tolerance = 0.02)

  f2 <- tempfile(fileext = ".png")
  suppressMessages(tage_save_plot(p, f2, width = 4, dpi = 72))
  expect_equal(dim(png::readPNG(f2))[2] / 72, 4, tolerance = 0.02)
  unlink(c(f, f2))
})

test_that("the forest plot accepts the clock registry as clocks_meta", {
  skip_if_not_installed("ggplot2")
  # predict_tAge() names its columns <normalisation>_<mode>_tAge, so a
  # list_clocks() table has to be matched through scaling + type.
  set.seed(3)
  n <- 60
  d <- data.frame(Genotype = factor(rep(c("WT", "KO"), each = n / 2), levels = c("WT", "KO")))
  d$scaled_diff_EN_tAge <- 0.5 * (d$Genotype == "KO") + stats::rnorm(n, sd = 0.4)
  d$yugene_diff_EN_tAge <- 0.6 * (d$Genotype == "KO") + stats::rnorm(n, sd = 0.4)

  clocks <- list_clocks(type = "EN", outcome = "Mortality",
                        species = "Multispecies", tissue = "Multi-Tissue")
  expect_false(any(clocks$filename %in% names(d)))

  p <- tage_clock_forest(d, clocks_meta = clocks, group_column = "Genotype",
                         reference_group = "WT")
  st <- attr(p, "tage_stats")
  expect_setequal(st$value_column, c("scaled_diff_EN_tAge", "yugene_diff_EN_tAge"))
  # Labels come from the registry fields, and the outcome drives the panel.
  labs <- levels(ggplot2::layer_data(p)$group |> factor())
  expect_true(all(grepl("Multispecies", .tage_clock_labels(clocks, clocks$filename))))
  expect_true(all(st$value_column %in% names(d)))
  .render(p)
})

test_that("registry outcomes are canonicalised to the figure keys", {
  expect_equal(.tage_canonical_outcome(c("Normalized age", "Mortality", "chronological", "other")),
               c("NormalizedAge", "Mortality", "Chronological", "other"))
})

test_that("two registry rows on one prediction column is an error", {
  skip_if_not_installed("ggplot2")
  set.seed(4)
  d <- data.frame(Genotype = factor(rep(c("WT", "KO"), each = 10), levels = c("WT", "KO")))
  d$scaled_diff_EN_tAge <- stats::rnorm(20)
  clocks <- list_clocks(type = "EN", scaling = "Scaled", species = "Multispecies",
                        tissue = "Multi-Tissue")   # Chronological + Mortality + Normalized age
  expect_gt(nrow(clocks), 1)
  expect_error(
    tage_clock_forest(d, clocks_meta = clocks, group_column = "Genotype", reference_group = "WT"),
    "same prediction column"
  )
})
