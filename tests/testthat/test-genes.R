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

test_that("the identifier type is detected from the gene table", {
  eset <- .tage_example_eset()                     # Ensembl row names
  gt <- tAge:::.read_gene_table(get_metadata_dir(), "mouse")

  d <- tAge:::.detect_gene_mapping_type(rownames(eset), gt)
  expect_equal(d$type, "Ensembl")
  expect_false(d$strip_version)
  expect_gt(d$n_matched, 10000)

  symbols <- head(gt$Gene.Symbol[!is.na(gt$Gene.Symbol)], 500)
  expect_equal(tAge:::.detect_gene_mapping_type(symbols, gt)$type, "Gene.Symbol")

  entrez <- head(gt$Entrez[!is.na(gt$Entrez)], 500)
  expect_equal(tAge:::.detect_gene_mapping_type(entrez, gt)$type, "Entrez")

  expect_error(tAge:::.detect_gene_mapping_type(c("foo", "bar"), gt), "None of the row names match")
  expect_error(tAge:::.detect_gene_mapping_type(symbols, gt, requested = "RefSeq"),
               "not available")
})

test_that("Ensembl version suffixes are stripped when that is what matches", {
  eset <- .tage_example_eset()
  eset <- filter_genes(eset, verbose = FALSE)
  versioned <- eset
  rownames(versioned) <- paste0(rownames(eset), ".", seq_len(nrow(eset)) %% 7 + 1)

  gt <- tAge:::.read_gene_table(get_metadata_dir(), "mouse")
  d <- tAge:::.detect_gene_mapping_type(rownames(versioned), gt)
  expect_true(d$strip_version)
  expect_equal(d$type, "Ensembl")

  plain  <- map_genes(eset, "mouse", verbose = FALSE)
  strip  <- map_genes(versioned, "mouse", verbose = FALSE)
  expect_equal(rownames(strip), rownames(plain))
  expect_equal(Biobase::exprs(strip), Biobase::exprs(plain))
})

test_that("Entrez input maps onto itself and auto-detection agrees with the explicit type", {
  eset <- filter_genes(.tage_example_eset(), verbose = FALSE)
  explicit <- map_genes(eset, "mouse", gene_mapping_type = "Ensembl", verbose = FALSE)
  auto     <- map_genes(eset, "mouse", verbose = FALSE)
  expect_equal(Biobase::exprs(auto), Biobase::exprs(explicit))

  entrez_in <- explicit
  again <- map_genes(entrez_in, "mouse", gene_mapping_type = "Entrez", verbose = FALSE)
  expect_equal(rownames(again), rownames(explicit))
  expect_equal(Biobase::exprs(again), Biobase::exprs(explicit))
})

test_that("an unknown species is an error", {
  eset <- .tage_example_eset()
  expect_error(map_genes(eset, "rhesus", verbose = FALSE), "Unknown species")
  expect_error(tAge_preprocessing(eset, species = "rhesus", verbose = FALSE), "Unknown species")
})
