# Preprocessing tests — no Python or clock models required.

test_that("filter_genes removes low-expressed genes", {
  eset <- .tage_example_eset()
  n_before <- nrow(eset)
  filtered <- filter_genes(eset, count_threshold = 10, percent_threshold = 20,
                           verbose = FALSE)
  expect_true(nrow(filtered) < n_before)
  expect_s4_class(filtered, "ExpressionSet")
})

test_that("scale_eset scales per sample (each column mean 0, sd 1)", {
  eset <- .tage_example_eset()
  eset <- filter_genes(eset, verbose = FALSE)
  eset <- log_transform(eset, verbose = FALSE)
  scaled <- scale_eset(eset, verbose = FALSE)
  m <- Biobase::exprs(scaled)

  # Per-sample (column-wise) standardisation: R's scale() operates on columns.
  expect_equal(unname(colMeans(m)), rep(0, ncol(m)), tolerance = 1e-8)
  expect_equal(unname(apply(m, 2, stats::sd)), rep(1, ncol(m)), tolerance = 1e-8)
})

test_that("YuGene is invariant to per-sample scaling", {
  # YuGene is invariant to a per-sample positive-affine transform, so applying
  # it to scaled vs log-transformed data must give identical results.
  eset <- .tage_example_eset()
  eset <- filter_genes(eset, verbose = FALSE)
  eset <- log_transform(eset, verbose = FALSE)
  scaled <- scale_eset(eset, verbose = FALSE)

  yg_from_log    <- Biobase::exprs(YuGene(eset,   verbose = FALSE))
  yg_from_scaled <- Biobase::exprs(YuGene(scaled, verbose = FALSE))

  expect_equal(yg_from_scaled, yg_from_log, tolerance = 1e-10)
})

test_that("control_subtraction centres on all samples by default", {
  # With no reference group, each gene's median across samples must become ~0.
  eset <- .tage_example_eset()
  eset <- filter_genes(eset, verbose = FALSE)
  eset <- log_transform(eset, verbose = FALSE)
  eset <- scale_eset(eset, verbose = FALSE)

  centred <- control_subtraction(eset, verbose = FALSE)
  m <- Biobase::exprs(centred)
  gene_medians <- apply(m, 1, stats::median, na.rm = TRUE)
  expect_equal(unname(gene_medians), rep(0, nrow(m)), tolerance = 1e-8)
})

test_that(".align_to_gene_list pads missing genes with NA and orders to gene_list", {
  eset <- .tage_example_eset()
  eset <- filter_genes(eset, verbose = FALSE)
  eset <- map_genes(eset, species = "mouse", gene_mapping_type = "Ensembl",
                    verbose = FALSE)
  # Alignment runs after scaling in the real pipeline (double-typed matrix).
  eset <- log_transform(eset, verbose = FALSE)
  gene_list <- load_gene_list()

  aligned <- tAge:::.align_to_gene_list(eset, gene_list)

  # Rows are exactly the gene list, in order.
  expect_identical(rownames(aligned), as.character(gene_list))
  # Genes absent from the data are padded with NA (not zeros).
  missing <- setdiff(as.character(gene_list), rownames(eset))
  if (length(missing) > 0) {
    expect_true(all(is.na(Biobase::exprs(aligned)[missing[1], ])))
  }
})
