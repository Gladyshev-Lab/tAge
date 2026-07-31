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

# 2. Preprocess: filter, map genes, RLE-normalise, log, scale, YuGene, centre
tAge_eset <- tAge_preprocessing(
  eset,
  species = "mouse",
  gene_mapping_type = "Gene.Symbol"    # or "Ensembl"
)

# 3. Predict
results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "EN")
head(results)
```

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

## Interpreting the output

The units of the prediction depend on the clock **outcome**:

| Outcome | Column | Units | Notes |
|---|---|---|---|
| Chronological | `*_tAge` | age (years for human, months for rodents) | normalised age rescaled by species max lifespan |
| Mortality | `*_tAge` | `log10(hazard ratio)` | **not** an age; higher = higher expected mortality |
| Normalized age | `*_tAge` | fraction of max lifespan | reported on its native scale |

Only **chronological** clocks are rescaled to age units. Mortality clocks
output a `log10` hazard ratio and must **not** be interpreted in years —
tAge derives this automatically from the clock registry.

## Relative (differential) clocks and centring

All distributed clocks are **relative** (`_scaleddiff` / `_yugenediff`)
models: they operate on expression *centred against a reference group*.
`tAge_preprocessing()` centres the data automatically:

- **No reference group given (default):** centres on all samples (overall
  per-gene median). Use this for a general age/mortality readout of a cohort.
- **Reference group given:** pass `control_group_column` and
  `control_group_label` to centre on age-/sex-matched controls within your
  design (recommended when comparing treatment vs control).

## Supported species

Mouse, human, rat, monkey. Non-mouse species are mapped to mouse orthologs
internally, so all clocks operate in the mouse Entrez gene space.

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
