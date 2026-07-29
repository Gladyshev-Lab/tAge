# Package index

## Pipelines

End-to-end entry points.
[`tAge_preprocessing()`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
returns the list of ExpressionSets that
[`predict_tAge()`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md)
consumes.

- [`tAge_preprocessing()`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
  : Complete preprocessing pipeline for tAge analysis
- [`predict_tAge()`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge.md)
  : Predict transcriptomic age for multiple processed ExpressionSet
  objects
- [`predict_tAge_one()`](https://gladyshev-lab.github.io/tAge/reference/predict_tAge_one.md)
  : Predict transcriptomic age using pre-trained models
- [`tAge_by_group()`](https://gladyshev-lab.github.io/tAge/reference/tAge_by_group.md)
  : Run tAge pipeline separately per group factor (e.g., "tissue")

## Clock catalogue

Browse the registry of pre-trained models and fetch them from Zenodo.

- [`list_clocks()`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md)
  : List available transcriptomic clock models
- [`download_clocks()`](https://gladyshev-lab.github.io/tAge/reference/download_clocks.md)
  : Download clock models from Zenodo

## Preprocessing steps

The individual stages run by
[`tAge_preprocessing()`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md),
for when you need a non-standard order or want to inspect an
intermediate.

- [`filter_genes()`](https://gladyshev-lab.github.io/tAge/reference/filter_genes.md)
  : Filter genes based on expression thresholds
- [`map_genes()`](https://gladyshev-lab.github.io/tAge/reference/map_genes.md)
  : Map genes in an ExpressionSet using local CSV mapping tables
- [`RLE_normalization()`](https://gladyshev-lab.github.io/tAge/reference/RLE_normalization.md)
  : Perform RLE (Relative Log Expression) normalization
- [`log_transform()`](https://gladyshev-lab.github.io/tAge/reference/log_transform.md)
  : Apply log10 transformation to expression data
- [`scale_eset()`](https://gladyshev-lab.github.io/tAge/reference/scale_eset.md)
  : Scale expression data to have zero mean and unit variance
- [`YuGene()`](https://gladyshev-lab.github.io/tAge/reference/YuGene.md)
  : Perform YuGene normalization
- [`control_subtraction()`](https://gladyshev-lab.github.io/tAge/reference/control_subtraction.md)
  : Subtract reference group median from expression data (relative
  expression)

## Single-cell and pseudobulk

Pool single cells into pseudobulk samples, drop outliers, and move data
in from Seurat or h5ad.

- [`aggregate_pseudobulk()`](https://gladyshev-lab.github.io/tAge/reference/aggregate_pseudobulk.md)
  : Aggregate single-cell data into pseudobulk samples based on read
  coverage
- [`aggregate_on_obs_columns()`](https://gladyshev-lab.github.io/tAge/reference/aggregate_on_obs_columns.md)
  : Aggregate single-cell data into pseudobulk samples within obs column
  groups
- [`pseudobulk_summary()`](https://gladyshev-lab.github.io/tAge/reference/pseudobulk_summary.md)
  : Get summary statistics for pseudobulk samples
- [`remove_outliers()`](https://gladyshev-lab.github.io/tAge/reference/remove_outliers.md)
  : Remove outlier pseudobulk samples using PCA-based Mahalanobis
  distance
- [`load_h5ad_simple()`](https://gladyshev-lab.github.io/tAge/reference/load_h5ad_simple.md)
  : Load h5ad file using reticulate and anndata
- [`load_h5ad_to_seurat()`](https://gladyshev-lab.github.io/tAge/reference/load_h5ad_to_seurat.md)
  : Load AnnData h5ad file and convert to Seurat object
- [`subset_seurat_by_metadata()`](https://gladyshev-lab.github.io/tAge/reference/subset_seurat_by_metadata.md)
  : Subset Seurat object by metadata criteria

## Plotting

- [`tage_boxplot()`](https://gladyshev-lab.github.io/tAge/reference/tage_boxplot.md)
  : Box plot of tAge predictions with pairwise significance annotation
- [`plot_eset_density()`](https://gladyshev-lab.github.io/tAge/reference/plot_eset_density.md)
  : Plot density curves for ExpressionSet data

## Data and helpers

Build an ExpressionSet and reach the example data shipped with the
package.

- [`make_ExpressionSet()`](https://gladyshev-lab.github.io/tAge/reference/make_ExpressionSet.md)
  : Create ExpressionSet object from expression data and phenotype data
- [`load_example_expression_data()`](https://gladyshev-lab.github.io/tAge/reference/load_example_expression_data.md)
  : Load example expression data
- [`load_example_metadata()`](https://gladyshev-lab.github.io/tAge/reference/load_example_metadata.md)
  : Load example metadata
- [`load_gene_list()`](https://gladyshev-lab.github.io/tAge/reference/load_gene_list.md)
  : Load gene list
- [`get_package_data()`](https://gladyshev-lab.github.io/tAge/reference/get_package_data.md)
  : Get package data file paths
- [`get_metadata_dir()`](https://gladyshev-lab.github.io/tAge/reference/get_metadata_dir.md)
  : Get path to the package metadata directory
