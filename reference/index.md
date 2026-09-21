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
  : Predict transcriptomic age with one pre-trained model
- [`tAge_by_group()`](https://gladyshev-lab.github.io/tAge/reference/tAge_by_group.md)
  : Run tAge pipeline separately per group factor (e.g., "tissue")

## Clock catalogue

Browse the registry of pre-trained models and fetch them from Zenodo.

- [`list_clocks()`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md)
  : List available transcriptomic clock models
- [`download_clocks()`](https://gladyshev-lab.github.io/tAge/reference/download_clocks.md)
  : Download clock models from Zenodo
- [`tage_species()`](https://gladyshev-lab.github.io/tAge/reference/tage_species.md)
  : Species supported by the clocks

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
  : Remove outlier samples
- [`load_h5ad_simple()`](https://gladyshev-lab.github.io/tAge/reference/load_h5ad_simple.md)
  : Load h5ad file using reticulate and anndata
- [`load_h5ad_to_seurat()`](https://gladyshev-lab.github.io/tAge/reference/load_h5ad_to_seurat.md)
  : Load AnnData h5ad file and convert to Seurat object
- [`subset_seurat_by_metadata()`](https://gladyshev-lab.github.io/tAge/reference/subset_seurat_by_metadata.md)
  : Subset Seurat object by metadata criteria

## Statistics

The group comparisons of the TACO / tClock application: marginal-mean
contrasts for elastic net clocks, REML meta-regression for Bayesian
ridge clocks, module-clock effect sizes.

- [`tage_compare_groups()`](https://gladyshev-lab.github.io/tAge/reference/tage_compare_groups.md)
  : Compare tAge between experimental groups
- [`tage_regress_continuous()`](https://gladyshev-lab.github.io/tAge/reference/tage_regress_continuous.md)
  : Regress tAge on a continuous predictor
- [`tage_module_stats()`](https://gladyshev-lab.github.io/tAge/reference/tage_module_stats.md)
  : Module-clock effect sizes and p-values
- [`tage_adjust_covariates()`](https://gladyshev-lab.github.io/tAge/reference/tage_adjust_covariates.md)
  : Covariate-adjusted tAge values for plotting
- [`tage_significance_stars()`](https://gladyshev-lab.github.io/tAge/reference/tage_significance_stars.md)
  : Significance stars for tAge statistics

## Plotting

Figures built directly on the statistics, plus QC plots.

- [`tage_clock_forest()`](https://gladyshev-lab.github.io/tAge/reference/tage_clock_forest.md)
  : Forest plot of clock effects
- [`tage_module_heatmap()`](https://gladyshev-lab.github.io/tAge/reference/tage_module_heatmap.md)
  : Heatmap of module-clock effects
- [`tage_boxplot()`](https://gladyshev-lab.github.io/tAge/reference/tage_boxplot.md)
  : Box plot of tAge predictions with pairwise significance annotation
- [`plot_eset_density()`](https://gladyshev-lab.github.io/tAge/reference/plot_eset_density.md)
  : Plot density curves for ExpressionSet data
- [`load_module_functions()`](https://gladyshev-lab.github.io/tAge/reference/load_module_functions.md)
  : Module colour to biological function map
- [`tage_fig_size()`](https://gladyshev-lab.github.io/tAge/reference/tage_fig_size.md)
  : Size a tAge figure was designed for
- [`tage_save_plot()`](https://gladyshev-lab.github.io/tAge/reference/tage_save_plot.md)
  : Save a tAge figure at the size it was designed for

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
