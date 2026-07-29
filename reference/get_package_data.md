# Get package data file paths

This function provides access to example data files included with the
tAge package. The data files are located in the inst/extdata directory
and include: - Exprs_example.csv: Example expression data -
Metadata_example.csv: Example metadata - Gene_list_all_4.6.txt: Gene
list for analysis

## Usage

``` r
get_package_data(filename)
```

## Arguments

- filename:

  Character string specifying the name of the data file to retrieve.

## Value

Character string with the full path to the requested data file.

## Examples

``` r
# Get path to example expression data
expr_path <- get_package_data("Exprs_example.csv")

# Get path to example metadata
meta_path <- get_package_data("Metadata_example.csv")

# Get path to gene list
gene_path <- get_package_data("Gene_list_all_4.6.txt")
```
