# Transcriptomic age from single-cell RNA-seq (pseudobulk) with tAge

## Overview

The clocks are trained on bulk RNA-seq. To apply them to single-cell
data, aggregate cells into **pseudobulk** samples, then run the standard
bulk pipeline. `tAge` aggregates by *coverage* — each pseudobulk sample
pools enough cells to reach a target read depth — which is more robust
than a fixed cells-per-sample cut-off.

This walkthrough uses Tabula Muris Senis
([GSM4505404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4505404)).
Requires `Seurat` (and `SeuratDisk` or `reticulate`+`anndata` to read
`.h5ad`).

## 1. Load single-cell data into a Seurat object

``` r

library(tAge)
library(Seurat)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = ".venv/bin/python")

# Read an .h5ad via anndata and build a Seurat object
ad    <- import("anndata")
adata <- ad$read_h5ad("GSM4505404_tabula-muris-senis-droplet-official-raw-obj.h5ad")

counts <- Matrix::t(adata$X)
rownames(counts) <- adata$var_names$to_list()
colnames(counts) <- adata$obs_names$to_list()

seurat_obj <- CreateSeuratObject(
  counts    = counts,
  meta.data = as.data.frame(adata$obs),
  project   = "TabulaMurisSenis"
)
```

`tAge` also provides
[`load_h5ad_to_seurat()`](https://gladyshev-lab.github.io/tAge/reference/load_h5ad_to_seurat.md)
as a convenience wrapper.

## 2. (Optional) subset

``` r

selected_tissues <- c("Tongue", "Kidney", "Liver", "Heart")
seurat_subset <- subset(seurat_obj, subset = tissue %in% selected_tissues)
```

## 3. Coverage-based pseudobulk aggregation

Aggregate within each unique combination of grouping columns, pooling
cells until each pseudobulk sample reaches `coverage_threshold` reads.

``` r

eset <- aggregate_on_obs_columns(
  seurat_obj       = seurat_subset,
  obs_column_names = c("age", "mouse.id", "tissue"),
  coverage_threshold = 1e7,
  verbose = TRUE
)

# Simpler variant: aggregate the whole object by coverage only
# eset <- aggregate_pseudobulk(seurat_subset, coverage_threshold = 1e6)
```

## 4. Remove outliers

``` r

eset_clean <- remove_outliers(eset, split_by = "tissue")
plot_eset_density(eset_clean)
```

## 5. Predict per tissue

[`tAge_by_group()`](https://gladyshev-lab.github.io/tAge/reference/tAge_by_group.md)
runs preprocessing and prediction separately within each level of
`split_by` (recommended for multi-tissue data).

``` r

# Download clocks (see the bulk vignette / ?download_clocks)
clocks <- download_clocks(
  list_clocks(type = "EN", outcome = "Chronological",
              species = "Multispecies", tissue = "Multi-Tissue"),
  dest_dir = "clocks"
)
model_paths <- list(
  scaled_diff = clocks$path[clocks$scaling == "Scaled"],
  yugene_diff = clocks$path[clocks$scaling == "YuGene"]
)

results <- tAge_by_group(
  eset        = eset_clean,
  split_by    = "tissue",
  model_paths = model_paths,
  species     = "mouse",
  mode        = "EN"
)
```

## 6. Visualise

``` r

results$age_numeric <- as.numeric(gsub("m", "", results$age))

tage_boxplot(
  results,
  x_var        = "age",
  y_var        = "scaled_diff_EN_tAge",
  subgroup_var = "tissue",
  stat_method  = "wilcox.test",
  theme_type   = "bw",
  title        = "Predicted tAge by age and tissue",
  xlab         = "Chronological age",
  ylab         = "Predicted tAge (months)"
)
```

See the bulk vignette
([`vignette("tage-bulk")`](https://gladyshev-lab.github.io/tAge/articles/tage-bulk.md))
for how to interpret the output and for the role of reference groups
(centring).
