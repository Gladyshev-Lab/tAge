# Module colour to biological function map

Reads the bundled annotation that names each co-expression module, e.g.
`"blue"` -\> `"Respiration/Mitochondrial translation"`.

## Usage

``` r
load_module_functions(version = c("4.6", "5.4", "human"))
```

## Arguments

- version:

  Module set: `"4.6"` (default), `"5.4"` or `"human"`.

## Value

A named character vector, module colour to function, or an empty vector
when the annotation is unavailable.

## Examples

``` r
head(load_module_functions("4.6"))
#>                                   black                                    blue 
#>         "Myogenesis/Muscle contraction" "Respiration/Mitochondrial translation" 
#>                                   brown                                    cyan 
#>  "Chromatin modification/Transcription"                       "Heme metabolism" 
#>                               darkgreen                                darkgrey 
#>          "Cell adhesion/VEGF signaling"     "Protein processing/mTOR signaling" 
```
