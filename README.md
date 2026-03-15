# tAge

<p align="center">
    <img src="tAge_logo.png" alt="tAge logo" width="300"/>
</p>

R package for transcriptomic biological age prediction from gene expression data.

## Installation

```r
devtools::install_github("Gladyshev-Lab/tAge")
```

### Python dependency

tAge uses Python for model prediction. Set up a virtual environment with required packages:

```bash
python -m venv .venv
.venv/bin/pip install joblib pandas scikit-learn
```

Then in R, before running predictions:

```r
Sys.setenv(RETICULATE_PYTHON = ".venv/bin/python")
```

## Quick start

### Bulk RNA-seq

```r
library(tAge)

# Load data
exprs_data <- read.csv("expression_matrix.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)
eset <- make_ExpressionSet(exprs_data, metadata)

# Preprocess
tAge_eset <- tAge_preprocessing(
  eset,
  species = "mouse",
  gene_mapping_type = "Gene.Symbol",
  control_group_column = "treatment",
  control_group_label = "control"
)

# Predict
model_paths <- list(
  scaled_diff = "path/to/EN_scaleddiff.pkl",
  yugene_diff = "path/to/EN_yugenediff.pkl"
)
Sys.setenv(RETICULATE_PYTHON = ".venv/bin/python")

results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "EN")
```

### Single-cell (pseudobulk)

```r
library(Seurat)

# Pseudobulk aggregation by sample and tissue
eset <- aggregate_on_obs_columns(
  seurat_obj,
  obs_column_names = c("sample_id", "tissue"),
  coverage_threshold = 1e7
)

# Remove outliers
eset_clean <- remove_outliers(eset, split_by = "tissue")

# Run tAge per tissue
results <- tAge_by_group(
  eset_clean,
  split_by = "tissue",
  model_paths = model_paths,
  species = "mouse",
  mode = "EN",
  control_group_column = "age",
  control_group_label = "young"
)
```

## Supported species

Mouse, human, rat, monkey. Non-mouse species are mapped to mouse orthologs internally.

## Citation

[Transcriptomic Hallmarks of Mortality Reveal Universal and Specific Mechanisms of Aging, Chronic Disease, and Rejuvenation (Tyshkovskiy et al. 2024)](https://www.biorxiv.org/content/10.1101/2024.07.04.601982v1)

## License

MIT
