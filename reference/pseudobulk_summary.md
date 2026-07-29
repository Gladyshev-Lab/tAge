# Get summary statistics for pseudobulk samples

Calculate and display summary statistics for pseudobulk ExpressionSet
objects.

## Usage

``` r
pseudobulk_summary(eset, group_vars = NULL)
```

## Arguments

- eset:

  An ExpressionSet object created from pseudobulk aggregation.

- group_vars:

  Character vector of metadata columns to group by for summary. Default
  is NULL (no grouping).

## Value

A data frame with summary statistics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic summary
summary_stats <- pseudobulk_summary(eset)

# Summary by tissue
summary_stats <- pseudobulk_summary(eset, group_vars = "tissue")

# Summary by multiple variables
summary_stats <- pseudobulk_summary(eset, group_vars = c("tissue", "age"))
} # }
```
