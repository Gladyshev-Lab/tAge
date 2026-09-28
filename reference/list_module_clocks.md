# List available module clock models

Returns the registry of published module clocks: one elastic net clock
per co-expression module (plus `"allmodulegenes"`, trained on the genes
of all modules), for the rodent and the multispecies module sets. Pass a
returned (filtered) data frame to
[`download_clocks`](https://gladyshev-lab.github.io/tAge/reference/download_clocks.md),
which fetches the archive each set is published as.

## Usage

``` r
list_module_clocks(outcome = NULL, species = NULL, color = NULL)
```

## Arguments

- outcome:

  Character. "Chronological" or "Mortality". Default NULL.

- species:

  Character. Module set: "Rodents" or "Multispecies". Default NULL.

- color:

  Character. Module colour(s), e.g. "blue". Default NULL.

## Value

A data frame with columns `filename`, `type`, `outcome`, `species`,
`tissue`, `scaling`, `color`, `function` (the module's annotated
biological process) and `archive` (the Zenodo archive holding the file).

## Examples

``` r
list_module_clocks(outcome = "Mortality", species = "Rodents")
#>                                                                      filename
#> 1  EN_Moduleclock_Mortality_Rodents_Multitissue_allmodulegenes_scaleddiff.pkl
#> 2           EN_Moduleclock_Mortality_Rodents_Multitissue_black_scaleddiff.pkl
#> 3            EN_Moduleclock_Mortality_Rodents_Multitissue_blue_scaleddiff.pkl
#> 4           EN_Moduleclock_Mortality_Rodents_Multitissue_brown_scaleddiff.pkl
#> 5            EN_Moduleclock_Mortality_Rodents_Multitissue_cyan_scaleddiff.pkl
#> 6       EN_Moduleclock_Mortality_Rodents_Multitissue_darkgreen_scaleddiff.pkl
#> 7        EN_Moduleclock_Mortality_Rodents_Multitissue_darkgrey_scaleddiff.pkl
#> 8      EN_Moduleclock_Mortality_Rodents_Multitissue_darkorange_scaleddiff.pkl
#> 9         EN_Moduleclock_Mortality_Rodents_Multitissue_darkred_scaleddiff.pkl
#> 10          EN_Moduleclock_Mortality_Rodents_Multitissue_green_scaleddiff.pkl
#> 11    EN_Moduleclock_Mortality_Rodents_Multitissue_greenyellow_scaleddiff.pkl
#> 12         EN_Moduleclock_Mortality_Rodents_Multitissue_grey60_scaleddiff.pkl
#> 13     EN_Moduleclock_Mortality_Rodents_Multitissue_lightgreen_scaleddiff.pkl
#> 14    EN_Moduleclock_Mortality_Rodents_Multitissue_lightyellow_scaleddiff.pkl
#> 15        EN_Moduleclock_Mortality_Rodents_Multitissue_magenta_scaleddiff.pkl
#> 16   EN_Moduleclock_Mortality_Rodents_Multitissue_midnightblue_scaleddiff.pkl
#> 17         EN_Moduleclock_Mortality_Rodents_Multitissue_orange_scaleddiff.pkl
#> 18           EN_Moduleclock_Mortality_Rodents_Multitissue_pink_scaleddiff.pkl
#> 19            EN_Moduleclock_Mortality_Rodents_Multitissue_red_scaleddiff.pkl
#> 20         EN_Moduleclock_Mortality_Rodents_Multitissue_salmon_scaleddiff.pkl
#> 21        EN_Moduleclock_Mortality_Rodents_Multitissue_skyblue_scaleddiff.pkl
#> 22            EN_Moduleclock_Mortality_Rodents_Multitissue_tan_scaleddiff.pkl
#> 23      EN_Moduleclock_Mortality_Rodents_Multitissue_turquoise_scaleddiff.pkl
#> 24         EN_Moduleclock_Mortality_Rodents_Multitissue_yellow_scaleddiff.pkl
#>    type   outcome species       tissue scaling          color
#> 1    EN Mortality Rodents Multi-Tissue  Scaled allmodulegenes
#> 2    EN Mortality Rodents Multi-Tissue  Scaled          black
#> 3    EN Mortality Rodents Multi-Tissue  Scaled           blue
#> 4    EN Mortality Rodents Multi-Tissue  Scaled          brown
#> 5    EN Mortality Rodents Multi-Tissue  Scaled           cyan
#> 6    EN Mortality Rodents Multi-Tissue  Scaled      darkgreen
#> 7    EN Mortality Rodents Multi-Tissue  Scaled       darkgrey
#> 8    EN Mortality Rodents Multi-Tissue  Scaled     darkorange
#> 9    EN Mortality Rodents Multi-Tissue  Scaled        darkred
#> 10   EN Mortality Rodents Multi-Tissue  Scaled          green
#> 11   EN Mortality Rodents Multi-Tissue  Scaled    greenyellow
#> 12   EN Mortality Rodents Multi-Tissue  Scaled         grey60
#> 13   EN Mortality Rodents Multi-Tissue  Scaled     lightgreen
#> 14   EN Mortality Rodents Multi-Tissue  Scaled    lightyellow
#> 15   EN Mortality Rodents Multi-Tissue  Scaled        magenta
#> 16   EN Mortality Rodents Multi-Tissue  Scaled   midnightblue
#> 17   EN Mortality Rodents Multi-Tissue  Scaled         orange
#> 18   EN Mortality Rodents Multi-Tissue  Scaled           pink
#> 19   EN Mortality Rodents Multi-Tissue  Scaled            red
#> 20   EN Mortality Rodents Multi-Tissue  Scaled         salmon
#> 21   EN Mortality Rodents Multi-Tissue  Scaled        skyblue
#> 22   EN Mortality Rodents Multi-Tissue  Scaled            tan
#> 23   EN Mortality Rodents Multi-Tissue  Scaled      turquoise
#> 24   EN Mortality Rodents Multi-Tissue  Scaled         yellow
#>                                     function              archive
#> 1                                All Modules Rodent module clocks
#> 2              Myogenesis/Muscle contraction Rodent module clocks
#> 3      Respiration/Mitochondrial translation Rodent module clocks
#> 4       Chromatin modification/Transcription Rodent module clocks
#> 5                            Heme metabolism Rodent module clocks
#> 6               Cell adhesion/VEGF signaling Rodent module clocks
#> 7          Protein processing/mTOR signaling Rodent module clocks
#> 8                   Lipid met/PPAR signaling Rodent module clocks
#> 9             Cholesterol met/mTOR signaling Rodent module clocks
#> 10                                Cell cycle Rodent module clocks
#> 11                      ECM organization/EMT Rodent module clocks
#> 12       Apoptosis/Nrf2 signaling/Proteasome Rodent module clocks
#> 13                      Interferon signaling Rodent module clocks
#> 14                 Nrf2 signaling/Proteasome Rodent module clocks
#> 15 Amino acid met/Coagulation/Xenobiotic met Rodent module clocks
#> 16                             mRNA splicing Rodent module clocks
#> 17                 Fatty acid met/Peroxisome Rodent module clocks
#> 18        Adaptive immunity/T cell signaling Rodent module clocks
#> 19                      Translation/Ribosome Rodent module clocks
#> 20    Cholesterol met/Platelet degranulation Rodent module clocks
#> 21                      Heat stress response Rodent module clocks
#> 22       Oxidative phosphorylation/TCA cycle Rodent module clocks
#> 23              Innate immunity/Inflammation Rodent module clocks
#> 24    Chromatin modification/Transcription 2 Rodent module clocks
```
