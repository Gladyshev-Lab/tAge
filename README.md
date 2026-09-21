# tAge

<p align="center">
    <img src="man/figures/logo-full.svg" alt="tAge logo" width="300"/>
</p>

R package for transcriptomic biological age prediction from gene expression
data, implementing the transcriptomic clocks from
[Tyshkovskiy et al. (2026), *Nature*](https://doi.org/10.1038/s41586-026-10542-3).

The clocks predict two kinds of biological age from bulk RNA-seq (or
pseudobulk single-cell) expression: **chronological age** and **expected
mortality**. tAge handles preprocessing, normalisation, gene mapping across
species, and prediction with the pre-trained Elastic Net (EN) and Bayesian
Ridge (BR) models.

## Installation

```r
# install.packages("devtools")
devtools::install_github("Gladyshev-Lab/tAge")
```

### Python dependency

tAge calls Python (via [`reticulate`](https://rstudio.github.io/reticulate/))
to run the model predictions. Create a virtual environment with the packages
listed in `inst/python/requirements.txt`:

```bash
python -m venv .venv
.venv/bin/pip install -r inst/python/requirements.txt   # joblib pandas scikit-learn numpy
```

Point R at that interpreter **before** calling `predict_tAge()`:

```r
# Linux / macOS
Sys.setenv(RETICULATE_PYTHON = ".venv/bin/python")
# Windows
Sys.setenv(RETICULATE_PYTHON = ".venv/Scripts/python.exe")
```

## Clock models

Models are published on Zenodo
([record 18763485](https://zenodo.org/records/18763485)). You can list and
download them directly from R:

```r
library(tAge)

# Browse available clocks (filter by type / outcome / species / tissue / scaling)
list_clocks(type = "EN", outcome = "Mortality")

# Download a chosen set; returns the table with a `path` column
clocks <- list_clocks(type = "EN", outcome = "Mortality",
                      species = "Multispecies", tissue = "Multi-Tissue")
clocks <- download_clocks(clocks, dest_dir = "clocks")

model_paths <- list(
  scaled_diff = clocks$path[clocks$scaling == "Scaled"],
  yugene_diff = clocks$path[clocks$scaling == "YuGene"]
)
```

## Quick start

### Bulk RNA-seq

```r
library(tAge)
Sys.setenv(RETICULATE_PYTHON = ".venv/bin/python")

# 1. Load raw counts (genes x samples) and sample metadata (samples x variables)
exprs_data <- read.csv("expression_matrix.csv", row.names = 1)
metadata   <- read.csv("metadata.csv", row.names = 1)
eset <- make_ExpressionSet(exprs_data, metadata)

# 2. Preprocess: filter, map genes, RLE-normalise, log, scale, YuGene, centre.
#    Gene identifiers (Ensembl / Gene.Symbol / Entrez) are detected automatically.
#    With several tissues, preprocess and centre each on its own controls:
tAge_eset <- tAge_preprocessing(
  eset,
  species = "mouse",
  split_by = "Tissue",
  control_group_column = "Genotype", control_group_label = "WT"
)

# 3. Predict (the species was recorded by tAge_preprocessing)
results <- predict_tAge(tAge_eset, model_paths, mode = "EN")
head(results)
attr(results, "tage_units")   # unit of every prediction column
```

`species` is the species of the *samples* and only sets the maximum lifespan
used to express chronological-age clocks in months (rodents) or years
(primates); it is not the species group a clock was trained on. `age_units`
forces months or years, and `normalized_age = "percent"` reports normalized-age
clocks as a percentage of the expected maximum lifespan, as the paper does.

### Single-cell (pseudobulk)

```r
library(Seurat)

eset <- aggregate_on_obs_columns(
  seurat_obj,
  obs_column_names = c("sample_id", "tissue"),
  coverage_threshold = 1e7
)
eset_clean <- remove_outliers(eset, split_by = "tissue")

results <- tAge_by_group(
  eset_clean,
  split_by = "tissue",
  model_paths = model_paths,
  species = "mouse",
  mode = "EN",
  control_group_column = "age",
  control_group_label  = "young"
)
```

## Statistics

`tage_compare_groups()` and friends implement the statistical models of the
paper's clock analyses; the R and Python packages agree numerically.
Elastic net clocks are compared with estimated marginal-mean contrasts of
`value ~ group + covariates` (`emmeans` on an `lm`); Bayesian ridge clocks go
through a REML meta-regression (`metafor::rma.uni`) that weights each sample by
its own prediction standard deviation and reports z-tests.

```r
# BR predictions carry a matching <normalisation>_BR_tAge_sd column.
results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "BR")

tage_compare_groups(
  results,
  value_columns   = "yugene_diff_BR_tAge",
  group_column    = "Genotype",
  reference_group = "WT",
  covariates      = "Sex",
  split_by        = "Tissue",
  se_columns      = "yugene_diff_BR_tAge_sd",   # BR clocks only
  p_adjust        = "BH"
)
```

One row per contrast, with `estimate` always meaning `group2 - group1`, plus
`se`, a 95% interval in `ci_low` / `ci_high`, `statistic`, `df`, `p_value`,
`p_adjusted` and a `label` holding the usual `*** ** * ^` stars.

| Function | Description |
|---|---|
| `tage_compare_groups()` | marginal-mean contrasts between groups |
| `tage_regress_continuous()` | slope of tAge against a numeric predictor |
| `tage_module_stats()` | per-module effect sizes and p-values for module clocks |
| `tage_adjust_covariates()` | covariate-adjusted values matching what the model tested |
| `tage_significance_stars()` | `*** ** * ^` labels |

Two figures are built directly on those tests, and return the statistics as the
`"tage_stats"` attribute:

| Function | Description |
|---|---|
| `tage_clock_forest()` | one clock per row, effect with 95% CI, filled = survives the correction |
| `tage_module_heatmap()` | module effects as a heatmap, modules × strata |

```r
p <- tage_clock_forest(
  results, clocks_meta = clocks,
  group_column = "Genotype", reference_group = "WT", split_by = "Tissue"
)
```

A `ggplot` carries no size of its own, so the figure functions record the size
they were designed for and grow it with the number of clocks, modules and
panels. `tage_save_plot()` uses that size; `tage_fig_size()` reads it back for a
knitr chunk:

```r
tage_fig_size(p)                       # width and height in inches
tage_save_plot(p, "forest.png")        # saved at the recorded size
tage_save_plot(p, "forest.pdf", width = 7)   # override just the width

# or fix the size up front
tage_clock_forest(..., width = 7, height = 4.5)
tage_clock_forest(..., row_height = 0.22, panel_width = 4.0)
tage_module_heatmap(..., cell_height = 0.28, cell_width = 0.9)
```

`variance_strata` chooses whether the residual variance comes from the two
compared groups (`"subset"`) or from every group in the stratum
(`"all_data"`); `p_adjust_scope` chooses the correction family —
`"within_column"` across comparisons and strata of one clock,
`"across_columns"` across clocks or modules within a comparison and stratum, or
`"global"`.

`tage_boxplot()` uses this engine by default (`stat_method = "emmeans"`); pass
`stat_method = "t.test"` (or another `ggpubr::stat_compare_means` method) for a
plain two-sample test.

## Interpreting the output

The units of the prediction depend on the clock **outcome**:

| Outcome | Column | Units | Notes |
|---|---|---|---|
| Chronological | `*_tAge` | months (rodents) or years (primates); `age_units` overrides | normalised age rescaled by species max lifespan (`tage_species()`) |
| Mortality | `*_tAge` | `log10(hazard ratio)` | **not** an age; higher = higher expected mortality |
| Normalized age | `*_tAge` | fraction of max lifespan, or `%` with `normalized_age = "percent"` | age / expected maximum lifespan of the experimental model |

Only **chronological** clocks are rescaled to age units. Mortality clocks
output a `log10` hazard ratio and must **not** be interpreted in years —
tAge derives this automatically from the clock registry (and from the file
name for models outside it). The `"tage_units"` attribute of the result
records the unit of every column.

## Relative (differential) clocks and centring

All distributed clocks are **relative** (`_scaleddiff` / `_yugenediff`)
models: they operate on expression *centred against a reference group*.
`tAge_preprocessing()` centres the data automatically:

- **No reference group given (default):** centres on all samples (overall
  per-gene median). Use this for a general age/mortality readout of a cohort.
- **Reference group given:** pass `control_group_column` and
  `control_group_label` to centre on age-/sex-matched controls within your
  design (recommended when comparing treatment vs control).
- **Several tissues, datasets or cell types:** pass `split_by`. The clocks
  were trained on expression centred *within each dataset and tissue*, so
  filtering, normalisation and centring must happen per stratum; a single
  reference pooled across tissues mixes tissue differences into the signal.
  `tAge_by_group()` does the same and predicts in one call.

The bundled example is the Klotho-knockout dataset of the paper (kidney and
skeletal muscle, WT vs KO). Note that the distributed clocks contain the
*Klotho* gene itself, whereas the paper's Klotho analysis used clocks
retrained without it; expect the KO effect here to be somewhat larger.

## Supported species

Mouse, human, rat, monkey (`tage_species()`). Non-mouse species are mapped to
mouse orthologs internally, so all clocks operate in the mouse Entrez gene
space. Input identifiers may be Ensembl (with or without version suffix), gene
symbols or Entrez IDs; `map_genes()` detects which.

## Citation

[Tyshkovskiy et al. (2026) Universal Transcriptomic Hallmarks of Mammalian Ageing and Mortality. *Nature*.](https://doi.org/10.1038/s41586-026-10542-3)

## License

This package and the associated transcriptomic clock models are distributed
under the MGB Open Access License 1.0. These materials may be used for
non-commercial academic purposes, subject to the license terms. Commercial
uses require a separate commercial license or agreement with Mass General
Brigham.

For questions about commercial use, please contact MGB Innovation / Business
Development & Licensing:

- innovation@partners.org
- licensing@partners.org
