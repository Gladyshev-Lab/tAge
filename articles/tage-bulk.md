# Transcriptomic age from bulk RNA-seq with tAge

## Overview

`tAge` applies the transcriptomic clocks from [Tyshkovskiy et
al. (2026), *Nature*](https://doi.org/10.1038/s41586-026-10542-3) to
bulk RNA-seq data, predicting **chronological age** (in age units) and
**expected mortality** (as a `log10` hazard ratio).

The workflow: build an `ExpressionSet`, preprocess it, and predict.

## 1. Build an ExpressionSet

We use the mouse example data bundled with the package (genes × samples
raw counts, plus sample metadata).

``` r

exprs_data <- load_example_expression_data()
metadata   <- load_example_metadata()

eset <- make_ExpressionSet(exprs_data, metadata, verbose = FALSE)
eset
#> ExpressionSet (storageMode: lockedEnvironment)
#> assayData: 57010 features, 24 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: Klo93K Klo94K ... RNA_122M (24 total)
#>   varLabels: Mouse.ID Genotype Sex Tissue
#>   varMetadata: labelDescription
#> featureData: none
#> experimentData: use 'experimentData(object)'
#> Annotation:
```

## 2. Preprocess

[`tAge_preprocessing()`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
filters low-expressed genes, maps them to the mouse Entrez gene space,
RLE-normalises, log-transforms, scales (per sample), applies the YuGene
transform, aligns to the clock gene list, and centres the data.

``` r

tAge_eset <- tAge_preprocessing(
  eset,
  species = "mouse",
  gene_mapping_type = "Ensembl",
  verbose = FALSE
)
#> calcNormFactors has been renamed to normLibSizes
names(tAge_eset)
#> [1] "RLE_normalized"  "log_transformed" "scaled"          "scaled_diff"    
#> [5] "yugene"          "yugene_diff"
```

Prediction uses the `scaled_diff` and `yugene_diff` representations.

## 3. Browse and download clocks

The published clocks live on Zenodo; list them from R (reads the bundled
registry):

``` r

head(list_clocks(type = "EN", outcome = "Mortality"))
#>                                               filename type   outcome
#> 1        EN_Mortality_Mouse_Multitissue_scaleddiff.pkl   EN Mortality
#> 2        EN_Mortality_Mouse_Multitissue_yugenediff.pkl   EN Mortality
#> 3 EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl   EN Mortality
#> 4 EN_Mortality_Multispecies_Multitissue_yugenediff.pkl   EN Mortality
#> 5            EN_Mortality_Rodents_Liver_scaleddiff.pkl   EN Mortality
#> 6            EN_Mortality_Rodents_Liver_yugenediff.pkl   EN Mortality
#>        species       tissue scaling lifespan_scaled
#> 1        Mouse Multi-Tissue  Scaled           FALSE
#> 2        Mouse Multi-Tissue  YuGene           FALSE
#> 3 Multispecies Multi-Tissue  Scaled           FALSE
#> 4 Multispecies Multi-Tissue  YuGene           FALSE
#> 5      Rodents        Liver  Scaled           FALSE
#> 6      Rodents        Liver  YuGene           FALSE
```

[`download_clocks()`](https://gladyshev-lab.github.io/tAge/reference/download_clocks.md)
fetches a chosen set and returns a `path` column:

``` r

clocks <- list_clocks(type = "EN", outcome = "Mortality",
                      species = "Multispecies", tissue = "Multi-Tissue")
clocks <- download_clocks(clocks, dest_dir = "clocks")

model_paths <- list(
  scaled_diff = clocks$path[clocks$scaling == "Scaled"],
  yugene_diff = clocks$path[clocks$scaling == "YuGene"]
)
```

## 4. Predict

``` r

Sys.setenv(RETICULATE_PYTHON = Sys.getenv("RETICULATE_PYTHON"))  # use your env

results <- predict_tAge(
  tAge_eset,
  model_paths = list(scaled_diff = .model_path),
  species = "mouse",
  mode = "EN"
)
head(results[, c("Genotype", "Tissue", "scaled_diff_EN_tAge")])
#>         Genotype Tissue scaled_diff_EN_tAge
#> Klo93K        WT Kidney          -0.5939931
#> Klo94K        WT Kidney          -0.6682321
#> Klo95K        WT Kidney          -0.5306280
#> Klo96K        WT Kidney          -0.1638941
#> Klo99K        WT Kidney          -0.5889994
#> Klo100K       WT Kidney          -0.6160492
```

These mortality values are `log10(hazard ratio)` — small numbers centred
near zero, **not** ages. See below.

## Interpreting the output

| Outcome | Units | Notes |
|----|----|----|
| Chronological | age (years / months) | normalised age rescaled by species max lifespan |
| Mortality | `log10(hazard ratio)` | **not** an age; higher = higher expected mortality |
| Normalized age | fraction of max lifespan | native scale |

`tAge` looks up each clock in the registry and rescales **only**
chronological clocks to age units.

## Reference groups (centring)

All distributed clocks are *relative*: they operate on expression
centred against a reference group.

- **Default (no reference group):** centre on all samples (overall
  per-gene median) — a general cohort readout.
- **Matched controls:** pass `control_group_column` and
  `control_group_label` to
  [`tAge_preprocessing()`](https://gladyshev-lab.github.io/tAge/reference/tAge_preprocessing.md)
  to centre on matched controls (recommended for treatment-vs-control
  designs),
  e.g. `control_group_column = "Genotype", control_group_label = "WT"`.

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] tAge_1.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6           knitr_1.51          rlang_1.3.0        
#>  [4] xfun_0.60           otel_0.2.0          png_0.1-9          
#>  [7] generics_0.1.4      textshaping_1.0.5   jsonlite_2.0.0     
#> [10] statmod_1.5.2       htmltools_0.5.9     ragg_1.5.2         
#> [13] sass_0.4.10         locfit_1.5-9.12     Biobase_2.72.0     
#> [16] rmarkdown_2.31      grid_4.6.1          evaluate_1.0.5     
#> [19] jquerylib_0.1.4     fastmap_1.2.0       yaml_2.3.12        
#> [22] lifecycle_1.0.5     compiler_4.6.1      fs_2.1.0           
#> [25] limma_3.68.4        Rcpp_1.1.2          edgeR_4.10.1       
#> [28] lattice_0.22-9      systemfonts_1.3.2   digest_0.6.39      
#> [31] R6_2.6.1            reticulate_1.46.0   bslib_0.11.0       
#> [34] Matrix_1.7-5        tools_4.6.1         BiocGenerics_0.58.1
#> [37] pkgdown_2.2.1       cachem_1.1.0        desc_1.4.3
```
