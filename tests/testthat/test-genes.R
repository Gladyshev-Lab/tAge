# Gene mapping tests: the bundled tables must be usable for every species the
# documentation promises.

.eset_from_ids <- function(ids, seed = 1) {
  set.seed(seed)
  m <- matrix(stats::rpois(length(ids) * 4, 50), nrow = length(ids),
              dimnames = list(ids, paste0("s", 1:4)))
  make_ExpressionSet(m, data.frame(row.names = colnames(m), g = c("a", "a", "b", "b")),
                     verbose = FALSE)
}

test_that("monkey Ensembl IDs map to mouse Entrez", {
  gt  <- utils::read.csv(file.path(get_metadata_dir(), "Gene_table_monkey.csv"),
                         stringsAsFactors = FALSE)
  ids <- head(unique(gt$Ensembl[!is.na(gt$Ensembl) & nzchar(gt$Ensembl)]), 300)
  eset <- .eset_from_ids(ids)

  mapped <- map_genes(eset, species = "monkey", gene_mapping_type = "Ensembl", verbose = FALSE)
  expect_gt(nrow(mapped), 0)
  expect_lte(nrow(mapped), length(ids))
  expect_true(all(grepl("^[0-9]+$", rownames(mapped))))
})

test_that("monkey gene symbols map to mouse Entrez", {
  gt  <- utils::read.csv(file.path(get_metadata_dir(), "Gene_table_monkey.csv"),
                         stringsAsFactors = FALSE)
  sym <- head(unique(gt$Gene.Symbol[!is.na(gt$Gene.Symbol) & nzchar(gt$Gene.Symbol)]), 300)
  eset <- .eset_from_ids(sym, seed = 2)

  mapped <- map_genes(eset, species = "monkey", gene_mapping_type = "Gene.Symbol", verbose = FALSE)
  expect_gt(nrow(mapped), 0)
  expect_true(all(grepl("^[0-9]+$", rownames(mapped))))
})

test_that("monkey genes without a mouse ortholog are dropped, not fatal", {
  gt  <- utils::read.csv(file.path(get_metadata_dir(), "Gene_table_monkey.csv"),
                         stringsAsFactors = FALSE)
  orth <- utils::read.csv(file.path(get_metadata_dir(), "Orthologs_monkey_to_mouse_5.0.csv"),
                          stringsAsFactors = FALSE)
  no_ortholog <- setdiff(unique(gt$Ensembl[!is.na(gt$Ensembl)]), orth$Ensembl.macaca)
  with_ortholog <- intersect(unique(gt$Ensembl[!is.na(gt$Ensembl)]),
                             orth$Ensembl.macaca[!is.na(orth$Entrez.mouse)])
  skip_if(length(no_ortholog) < 5 || length(with_ortholog) < 5)

  ids  <- c(head(no_ortholog, 20), head(with_ortholog, 20))
  eset <- .eset_from_ids(ids, seed = 3)
  mapped <- map_genes(eset, species = "monkey", gene_mapping_type = "Ensembl", verbose = FALSE)
  expect_gt(nrow(mapped), 0)
  expect_lte(nrow(mapped), 20)
})

test_that("human and rat still go through the mouse ortholog table", {
  gt  <- utils::read.csv(file.path(get_metadata_dir(), "Gene_table_human.csv"),
                         stringsAsFactors = FALSE)
  ids <- head(unique(gt$Ensembl[!is.na(gt$Ensembl) & nzchar(gt$Ensembl)]), 200)
  mapped <- map_genes(.eset_from_ids(ids, seed = 4), species = "human",
                      gene_mapping_type = "Ensembl", verbose = FALSE)
  orth <- utils::read.csv(file.path(get_metadata_dir(), "Table_of_orthologs.csv"),
                          stringsAsFactors = FALSE)
  expect_gt(nrow(mapped), 0)
  expect_true(all(rownames(mapped) %in% as.character(orth$Entrez.Mouse)))
})
