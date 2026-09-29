# tAge (development version)

## Breaking changes

* `tage_module_heatmap()` and `load_module_functions()` no longer default to
  the rodent module set: pass `module_set` (or `module_functions`). A module
  colour names a different module in each set -- 7 of the 8 colours shared by
  the rodent and multispecies sets differ, "blue" being
  Respiration/Mitochondrial translation in one and Myogenesis/Muscle
  contraction in the other -- so multispecies and human module clocks were
  labelled with rodent functions whenever the argument was left out.
  `module_functions = character(0)` labels the rows by module name only.

## Bug fixes

* `tage_adjust_covariates()` gains `group_column`. With it, the covariate
  effects removed are those of the model the statistics fit (`value ~ group +
  covariates`, one model per stratum, `lm` or the weighted meta-regression), so
  the difference between group means of the adjusted values is the tested
  estimate. The covariate-only model, still the default, attributes part of
  the group effect to covariates that are unevenly distributed across groups
  (tested KO - WT 1.22, plotted 0.91 in an unbalanced example).
  `tage_boxplot()` passes its groups, so its covariate-adjusted points now show
  what its brackets test.

* `aggregate_pseudobulk()` adds cells until a pseudobulk sample reaches
  `coverage_threshold` and then starts the next one, as
  `aggregate_on_obs_columns()` and the Python package do. It used to cut the
  cumulative read count at multiples of the threshold, so the first sample and
  about half of the others held fewer reads than the threshold, and a cell with
  more reads than the threshold left empty samples (0 cells, 0 counts) behind.
  Both functions gain `drop_incomplete` (as in Python) to drop the leftover last
  sample, and share one implementation that sums with a sparse indicator
  matrix. `aggregate_on_obs_columns()` now respects `verbose`.

* `tAge_preprocessing()` refuses input that is not raw counts -- missing,
  negative or non-integer values (TPM, CPM, log data), with the tolerance the
  Python package uses -- and a `control_group_column` that does not exist or a
  `control_group_label` that no sample carries; both used to fall back to
  centring on all samples with a warning. A stratum of `split_by` without
  reference samples is still centred on itself with a warning.
  `control_subtraction()` refuses a column that does not exist.

* `map_genes()` keeps Entrez IDs as text from integers when summing
  identifiers that collapse onto one gene. Grouping on the numbers named round
  IDs in scientific notation (rat 500000 became `"5e+05"`), so they matched
  nothing in the ortholog table and were dropped; rat 500000 is the ortholog of
  the clock gene 64945. Duplicates are summed with `rowsum()`.

* The statistics leave out a categorical covariate that has a single level in
  the data a model is fitted on (e.g. `Sex` in an all-male stratum) instead of
  skipping the stratum on `lm()`'s "contrasts can be applied only to factors
  with 2 or more levels". The Python package already fitted these strata; the
  two now agree.

# tAge 1.5.0

## Clock registry

* The registry is the Zenodo record 22166800: the 60 composite clocks and,
  new, the published module clocks. `list_module_clocks()` lists the 78
  elastic net module clocks of the rodent and multispecies sets (one per
  co-expression module plus `allmodulegenes`, chronological and mortality,
  Scaled normalisation) with their annotated function and the archive each
  set is published as. `download_clocks()` accepts its output: the archive is
  downloaded once, unpacked into `dest_dir`, and `path` points inside it.
  Downloads are checked to be a pickle or a zip archive. Module clocks are
  scaled like the composite clock of the same outcome.

* Module sets are named `"rodent"`, `"multispecies"` and `"human"`:
  `load_module_functions(module_set = )` and
  `tage_module_heatmap(module_set = )` replace the `version` /
  `modules_version` arguments. The bundled annotation files are
  `Module_to_function_map_<set>.csv`.

* Bundled data files renamed: `Gene_list_rodent_clocks.txt` (the gene list of
  the rodent clocks, read by `load_gene_list()`) and
  `metadata/Orthologs_monkey_to_mouse.csv`.

# tAge 1.4.0

## New features

* Species maximum lifespans are the values the clocks were trained with
  (AnAge): rat 3.8 years (was 4.2) and human 122 years (was 122.5); mouse
  4 and macaque 39 years are unchanged. Chronological-age predictions for
  rat samples are therefore ~10% lower than before, human ones ~0.4%.

* One figure style, shared with the Python package: `theme_tage()`, the
  colour tokens in `tage_colors()`, `tage_series_colors()` for group colours
  (the reference group in neutral grey, comparisons in colour-vision-safe
  slots) and `scale_fill_tage_diverging()` for signed effects. `tage_boxplot()`,
  `tage_clock_forest()`, `tage_module_heatmap()`, the outlier PCA plot and
  `plot_eset_density()` all draw with it: left-aligned title with subtitle
  and caption, recessive axes, one colour per clock outcome, blue-grey-red
  scale centred on zero. `tage_boxplot()` gains `subtitle` and `caption`;
  its `theme_type` is deprecated.

* `tAge_preprocessing(split_by = )` preprocesses each level of a phenoData
  column (tissue, dataset, cell type) on its own -- gene filtering,
  normalisation and reference centring within the stratum, on the stratum's
  own controls -- and combines the results. This is how the clocks were
  trained and applied in the paper; a single reference pooled across tissues
  mixes tissue differences into the signal. `tAge_by_group()` is now this
  plus `predict_tAge()`, and no longer swallows errors per stratum. The
  bulk vignette preprocesses the two-tissue Klotho example this way, with
  the paper's 25% gene-detection threshold.

* Gene identifiers are detected: `gene_mapping_type = "auto"` (the new
  default of `map_genes()` and `tAge_preprocessing()`) picks Ensembl, gene
  symbol or Entrez as the type with the most matches in the species' gene
  table; Ensembl version suffixes are stripped when that is what makes the
  IDs match. Entrez input is new. Nothing matching is an error that lists the
  match counts.

* The species is recorded in the ExpressionSets returned by
  `tAge_preprocessing()`; `predict_tAge()` / `predict_tAge_one()` take it
  from there (`species = NULL`). `species` is documented as the species of
  the samples, used only to rescale chronological-age clocks; an unknown
  species is an error instead of a silent factor of 1. `tage_species()`
  lists the supported species with their maximum lifespans and default units.

* `predict_tAge()` gains `age_units` ("auto": months for rodents, years for
  primates; or "months" / "years") and `normalized_age` ("fraction", the
  default, or "percent" as in the paper and TACO). The result carries a
  `"tage_units"` attribute naming the unit of every prediction column. Clocks
  outside the registry are recognised from their file names (Chronoage,
  Hazard / Mortality, Relage / NormalizedAge); unknown names are left on the
  native scale with a warning.

* `remove_outliers(method = )` adds the paper's two rules next to the
  Mahalanobis default: `"pca_iqr"` (PC1 or PC2 beyond 1.5 IQR, bulk data) and
  `"spearman_median"` (Spearman correlation with the group's median profile
  below 0.5, meta-dataset).

## Bug fixes

* `control_subtraction()` warns when the requested control label matches no
  sample (it used to fall back to all samples silently unless `verbose`).

## Housekeeping

* The Python bridge no longer silences every `UserWarning` for the whole
  session; the scikit-learn version warning is suppressed around the model
  load only.
* Base-package functions are imported explicitly (`R CMD check` NOTEs);
  `tage_boxplot()` no longer calls `library()`. `png` and `robustbase` are
  listed in Suggests. The test helper downloads models through
  `download_clocks()`.
* `tage_boxplot()` needs `ggpubr` only for the bracket layer it draws with
  it; the pkgdown workflow installs it for the vignettes.
* `download_clocks()` writes to `<file>.part` and renames only once the
  transfer is complete and checked, so an interrupted session cannot leave a
  truncated model under the real name. Non-ASCII characters in R sources are
  written as `\u` escapes (`R CMD check` warning).

# tAge 1.3.1

## Bug fixes

* `map_genes(species = "monkey")` failed on every input with "subscript out
  of bounds": macaque genes without a mouse ortholog were looked up with `[[`.
  Mapping is now vectorised and drops those genes, as for the other species.

* `tage_clock_forest(clocks_meta = list_clocks(...))` could not find the
  prediction columns: `predict_tAge()` names them `<normalisation>_<mode>_tAge`,
  not by model file. Registry rows are now matched through `scaling` and
  `type` (Scaled + EN -> `scaled_diff_EN_tAge`), labelled from the registry
  fields when there is no `name` column, and the registry's "Normalized age"
  outcome is recognised. Two rows landing on one column is an error with an
  explanation.

* `tage_adjust_covariates()` with a Bayesian ridge `se_column` and no
  `split_by` returned zero-centred residuals; the values are now put back on
  the tAge scale like every other branch.

* `tage_compare_groups()`, `tage_regress_continuous()`,
  `tage_module_stats()` and `tage_adjust_covariates()` no longer drop a
  stratum in silence. Every skipped clock/stratum is reported with a
  warning naming it and the reason (reference group absent, fewer than two
  groups, collinear covariates, a failed `lm` / `rma.uni` fit, ...).

* `download_clocks()` raises R's download timeout while it runs (new
  `timeout` argument, default 3600 s; the default 60 s aborted every
  Bayesian ridge model, 0.9-2.4 GB each), removes partial files instead of
  leaving them to be reported as "already present", and rejects files that
  are not pickles (Zenodo error pages).

* Normalized-age panels of `tage_clock_forest()` are labelled "fraction of
  maximum lifespan", the scale the package actually returns.

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
