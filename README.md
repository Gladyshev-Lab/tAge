# tAge: Transcriptomic Biological Age Analysis Package

<p align="center">
    <img src="tAge_logo.png" alt="drawing" width="300"/>
</p>

A comprehensive R package for transcriptomic age prediction and analysis, developed for the Gladyshev Lab at Harvard Medical School.

## Installation

Install from GitHub:

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install tAge
devtools::install_github("Gladyshev-Lab/tAge")
```

## Quick Start

```r
library(tAge)

# Load and preprocess data
expr_data <- load_example_expression_data()
meta_data <- load_example_metadata()
gene_list <- load_gene_list()

# Create ExpressionSet
eset <- make_ExpressionSet(expr_data, meta_data)

# Preprocess with tAge pipeline
processed <- tAge_preprocessing(eset, gene_list, species = "mouse")

# Predict transcriptomic age
results <- predict_tAge(processed$scaled, 
                        model_path = "path/to/model.pkl",
                        species = "mouse", 
                        mode = "EN")
```

## Features

- **Preprocessing**: Gene filtering, normalization (scaled/YuGene), differential expression
- **Prediction**: Elastic Net (EN) and Bayesian Ridge (BR) models
- **Visualization**: Density distribution and boxplots
- **Species support**: Mouse, human, rat, and custom gene sets

## Main Functions

- `tAge_preprocessing()` - Complete preprocessing pipeline
- `predict_tAge()` - Age prediction with pre-trained models
- `plot_tAge_results()` - Visualize predictions
- `filter_genes()` - Gene filtering and conversion

## Citation

If you use tAge in your research, please cite the Gladyshev Lab.

## License

MIT License. See `LICENSE` file for details.
