# Map genes in an ExpressionSet using local CSV mapping tables

Map genes in an ExpressionSet using local CSV mapping tables

## Usage

``` r
map_genes(eset, species, gene_mapping_type, verbose = TRUE)
```

## Arguments

- eset:

  An ExpressionSet object.

- species:

  One of "human", "mouse", "rat", "monkey".

- gene_mapping_type:

  One of "Ensembl" or "Gene.Symbol".

- verbose:

  Logical. Print progress messages. Default TRUE.

## Value

ExpressionSet with Entrez (mouse) gene IDs as rownames.
