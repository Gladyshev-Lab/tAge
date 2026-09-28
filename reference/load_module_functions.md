# Module colour to biological function map

Reads the bundled annotation that names each co-expression module, e.g.
`"blue"` -\> `"Respiration/Mitochondrial translation"`.

## Usage

``` r
load_module_functions(module_set = c("rodent", "multispecies", "human"))
```

## Arguments

- module_set:

  Module set: `"rodent"` (default) or `"multispecies"`, the two sets
  published with the clocks (see
  [`list_module_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_module_clocks.md)),
  or `"human"`.

## Value

A named character vector, module colour to function, or an empty vector
when the annotation is unavailable.

## Examples

``` r
head(load_module_functions("rodent"))
#>                                   black                                    blue 
#>         "Myogenesis/Muscle contraction" "Respiration/Mitochondrial translation" 
#>                                   brown                                    cyan 
#>  "Chromatin modification/Transcription"                       "Heme metabolism" 
#>                               darkgreen                                darkgrey 
#>          "Cell adhesion/VEGF signaling"     "Protein processing/mTOR signaling" 
```
