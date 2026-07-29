# List available transcriptomic clock models

Returns the registry of published clock models available on Zenodo, with
optional filtering. Pass a returned (filtered) data frame to
[`download_clocks`](https://gladyshev-lab.github.io/tAge/reference/download_clocks.md)
to download the corresponding files.

## Usage

``` r
list_clocks(
  type = NULL,
  outcome = NULL,
  species = NULL,
  tissue = NULL,
  scaling = NULL
)
```

## Arguments

- type:

  Character. Filter by model type: "EN" (Elastic Net) or "BR" (Bayesian
  Ridge). Default NULL (no filter).

- outcome:

  Character. Filter by prediction outcome: "Chronological", "Mortality",
  or "Normalized age". Default NULL.

- species:

  Character. Filter by species group: "Mouse", "Rodents", or
  "Multispecies". Default NULL.

- tissue:

  Character. Filter by tissue (e.g. "Multi-Tissue", "Liver"). Default
  NULL.

- scaling:

  Character. Filter by normalisation: "Scaled" or "YuGene". Default
  NULL.

## Value

A data frame with columns `filename`, `type`, `outcome`, `species`,
`tissue`, `scaling` and `lifespan_scaled` (whether the output is
rescaled to age units).

## Examples

``` r
# All Elastic Net mortality clocks
list_clocks(type = "EN", outcome = "Mortality")
#>                                               filename type   outcome
#> 1        EN_Mortality_Mouse_Multitissue_scaleddiff.pkl   EN Mortality
#> 2        EN_Mortality_Mouse_Multitissue_yugenediff.pkl   EN Mortality
#> 3 EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl   EN Mortality
#> 4 EN_Mortality_Multispecies_Multitissue_yugenediff.pkl   EN Mortality
#> 5            EN_Mortality_Rodents_Liver_scaleddiff.pkl   EN Mortality
#> 6            EN_Mortality_Rodents_Liver_yugenediff.pkl   EN Mortality
#> 7      EN_Mortality_Rodents_Multitissue_scaleddiff.pkl   EN Mortality
#> 8      EN_Mortality_Rodents_Multitissue_yugenediff.pkl   EN Mortality
#>        species       tissue scaling lifespan_scaled
#> 1        Mouse Multi-Tissue  Scaled           FALSE
#> 2        Mouse Multi-Tissue  YuGene           FALSE
#> 3 Multispecies Multi-Tissue  Scaled           FALSE
#> 4 Multispecies Multi-Tissue  YuGene           FALSE
#> 5      Rodents        Liver  Scaled           FALSE
#> 6      Rodents        Liver  YuGene           FALSE
#> 7      Rodents Multi-Tissue  Scaled           FALSE
#> 8      Rodents Multi-Tissue  YuGene           FALSE
```
