# Download clock models from Zenodo

Downloads the given clock models from the Zenodo record and returns the
input augmented with a `path` column pointing to the local files.
Composite clocks
([`list_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md))
are single files; module clocks
([`list_module_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_module_clocks.md))
come in one archive per module set, which is downloaded once and
unpacked into `dest_dir`.

## Usage

``` r
download_clocks(
  clocks,
  dest_dir = "clocks",
  record = .TAGE_ZENODO_RECORD,
  overwrite = FALSE,
  quiet = FALSE,
  timeout = 3600
)
```

## Arguments

- clocks:

  A data frame returned by
  [`list_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_clocks.md)
  or
  [`list_module_clocks`](https://gladyshev-lab.github.io/tAge/reference/list_module_clocks.md),
  or a character vector of composite clock file names.

- dest_dir:

  Directory to save the models into. Created if needed. Default
  "clocks".

- record:

  Character Zenodo record id. Defaults to the published record.

- overwrite:

  Logical. Re-download files that already exist. Default FALSE.

- quiet:

  Logical. Suppress progress messages. Default FALSE.

- timeout:

  Seconds allowed per file. R's default of 60 s aborts the Bayesian
  ridge models, which are 0.9-2.4 GB each; elastic net models are about
  1 MB and the module clock archives under 1 MB. Default 3600.

## Value

The `clocks` data frame with an added `path` column.

## Details

The whole record is about 64 GB, nearly all of it the 30 Bayesian ridge
models. A partial or failed download is removed rather than left on
disk, and every file is checked to be a pickle or a zip archive (Zenodo
answers some errors with an HTML page, which would otherwise be saved
under the model's name).

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

modules <- download_clocks(list_module_clocks(outcome = "Mortality", species = "Rodents"),
                           dest_dir = "clocks")
} # }
```
