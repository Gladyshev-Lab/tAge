# Map genes in an ExpressionSet using local CSV mapping tables

Translates the row names of `eset` to mouse Entrez IDs, the gene space
every clock operates in: identifiers of the given species are looked up
in the bundled gene table, and rat, human and macaque genes are then
carried over to mouse through 1:1 ortholog tables. Genes that do not map
are dropped; identifiers that collapse onto one Entrez ID are summed.

## Usage

``` r
map_genes(eset, species, gene_mapping_type = "auto", verbose = TRUE)
```

## Arguments

- eset:

  An ExpressionSet object.

- species:

  One of "human", "mouse", "rat", "monkey".

- gene_mapping_type:

  Identifier type of the row names: "Ensembl", "Gene.Symbol", "Entrez",
  or "auto" (default) to detect it as the type with the most matches in
  the gene table. Ensembl IDs may carry a version suffix
  (`ENSMUSG00000000001.4`); it is stripped when that is what makes them
  match.

- verbose:

  Logical. Print progress messages. Default TRUE.

## Value

ExpressionSet with Entrez (mouse) gene IDs as rownames.
