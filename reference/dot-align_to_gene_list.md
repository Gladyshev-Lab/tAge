# Align an ExpressionSet to a reference gene list

Reorders and subsets rows to match the reference gene list exactly.
Genes missing from the ExpressionSet are padded with `NA`. This is
intentional: at prediction time the trained clock model's imputer fills
these with the training-set median for each gene, which is the correct
neutral value (padding with zeros would not be). This matches the TACO
reference application, where absent genes remain `NA`.

## Usage

``` r
.align_to_gene_list(eset, gene_list)
```

## Arguments

- eset:

  An ExpressionSet object.

- gene_list:

  Character vector of reference gene identifiers.

## Value

ExpressionSet with rows ordered and padded to match gene_list exactly.
