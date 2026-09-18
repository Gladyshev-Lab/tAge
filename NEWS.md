# tAge 1.3.0

## New features

* Two publication-style figures built directly on the statistics, so a figure
  and the table behind it cannot drift apart: `tage_clock_forest()` draws one
  clock per row with its confidence interval, filled when it survives the
  multiplicity correction; `tage_module_heatmap()` draws module effects as
  modules × strata with a star per significant cell. Both return the statistics
  as the `"tage_stats"` attribute, and both accept a precomputed table via
  `stats =`. `load_module_functions()` reads the bundled module-to-function
  annotation used for the row labels.

* The statistics gain `ci_low` / `ci_high` and a `conf_level` argument. The
  critical value follows the test: normal for the Bayesian ridge
  meta-regression, t otherwise.

* Statistical tests matching the TACO / tClock reference application:
  `tage_compare_groups()` for marginal-mean contrasts between groups,
  `tage_regress_continuous()` for slopes against a numeric predictor,
  `tage_module_stats()` for module-clock heatmaps, `tage_adjust_covariates()`
  for the matching plotting values, and `tage_significance_stars()` for the
  `*** ** * ^` labels. Elastic net clocks use `lm` + `emmeans`; Bayesian ridge
  clocks use `metafor::rma.uni` weighted by the per-sample prediction standard
  deviation, reporting z-tests. `emmeans` and `metafor` are new imports.

* `predict_tAge()` and `predict_tAge_one()` gain `return_std`, on by default for
  `mode = "BR"`, adding a `<normalisation>_BR_tAge_sd` column. These are the
  weights the Bayesian ridge tests need, and were previously discarded.

* `tage_boxplot()` defaults to `stat_method = "emmeans"`, annotating brackets
  from `tage_compare_groups()` and supporting covariates, per-stratum models and
  Bayesian ridge weighting. Passing any other `stat_method` keeps the previous
  `ggpubr::stat_compare_means()` behaviour.

## Bug fixes

* The predictive standard deviation of chronological-age clocks is now rescaled
  by the species maximum lifespan along with the prediction itself. It was left
  on the normalised scale, which would have made meta-regression weights wrong
  by the square of the species factor.

# tAge 1.1.0

Corrects three bugs that produced **wrong predictions** in 1.0.0 / 1.0.1.
Analyses run with those versions should be repeated with 1.1.0.

## Bug fixes

* Species rescaling is now applied **only to chronological-age clocks**.
  Previously every prediction was multiplied by the species maximum-lifespan
  factor, which turned mortality output (`log10` hazard ratio) and
  normalized-age output into meaningless numbers. The clock type is taken from
  the model file name, matching the TACO reference application.

* Reference centring is **never skipped**. `control_subtraction()` used to
  return the data unchanged when no reference group was given. All distributed
  clocks are relative (`_scaleddiff` / `_yugenediff`) models trained on
  reference-centred expression, so predictions made without centring were
  invalid. `column_name` and `control_label` now default to `NULL`, which
  centres on all samples (overall per-gene median) — the TACO default.

* Genes absent from the input are padded with `NA` instead of `0` when aligning
  to the clock gene list. The trained model's imputer then fills them with the
  training-set median for that gene, which is the correct neutral value.

* `predict_tAge()` coerces the reticulate result to a `data.frame`, fixing a
  failure on reticulate/pandas versions that return a bare vector for a
  single-column result.

* `tage_boxplot()` is exported.

## New features

* Clock registry: `list_clocks()` browses the pre-trained models (filter by
  type, outcome, species, tissue and scaling) and `download_clocks()` fetches
  them from Zenodo record 18763485, returning the table with a `path` column.
  The registry ships with the package as `inst/extdata/clocks_metadata.csv`.

## Documentation

* Two vignettes replace the former Jupyter notebooks: `vignette("tage-bulk")`
  for bulk RNA-seq and `vignette("tage-singlecell")` for the pseudobulk
  single-cell workflow.

* pkgdown site at <https://gladyshev-lab.github.io/tAge/>.

* README rewritten: units of each clock outcome, the role of reference groups
  in relative clocks, supported species, and licensing.

* `scale_eset()` documents that scaling is **per sample** (column-wise, each
  sample scaled across genes). Per-gene standardisation happens separately
  inside the trained model, using training-set statistics.

## Testing

* `testthat` suite covering the clock registry, prediction and preprocessing.

# tAge 1.0.1

* README and LICENSE updated for the MGB Open Access License 1.0.
* Model availability notice and placeholder paths corrected in the tutorials.

# tAge 1.0.0

* First release, accompanying
  [Tyshkovskiy et al. (2026), *Nature*](https://doi.org/10.1038/s41586-026-10542-3).
* Superseded by 1.1.0 — see the bug fixes above before using any results from
  this version.
