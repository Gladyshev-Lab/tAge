# Download clock models from Zenodo

Downloads the given clock model files from the Zenodo record and returns
the input augmented with a `path` column pointing to the local files.

## Usage

``` r
download_clocks(
  clocks,
  dest_dir = "clocks",
  record = .TAGE_ZENODO_RECORD,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- clocks:

  Either a data frame returned by
  [`list_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md)
  (the `filename` column is used) or a character vector of model file
  names.

- dest_dir:

  Directory to save the models into. Created if needed. Default
  "clocks".

- record:

  Character Zenodo record id. Defaults to the published record.

- overwrite:

  Logical. Re-download files that already exist. Default FALSE.

- quiet:

  Logical. Suppress progress messages. Default FALSE.

## Value

The `clocks` data frame with an added `path` column.

## Examples

``` r
if (FALSE) { # \dontrun{
clocks <- list_clocks(type = "EN", outcome = "Mortality",
                      species = "Multispecies", tissue = "Multi-Tissue")
clocks <- download_clocks(clocks, dest_dir = "clocks")
model_paths <- list(
  scaled_diff = clocks$path[clocks$scaling == "Scaled"],
  yugene_diff = clocks$path[clocks$scaling == "YuGene"]
)
} # }
```
