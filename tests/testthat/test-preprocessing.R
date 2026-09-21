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

test_that("split_by preprocesses and centres each stratum on its own controls", {
  eset <- .tage_example_eset()
  proc <- tAge_preprocessing(
    eset, species = "mouse", split_by = "Tissue",
    control_group_column = "Genotype", control_group_label = "WT",
    verbose = FALSE
  )
  sd <- proc$scaled_diff
  expect_s4_class(sd, "ExpressionSet")
  expect_equal(colnames(sd), colnames(eset))            # sample order kept
  expect_equal(rownames(sd), as.character(load_gene_list()))

  pd <- Biobase::pData(sd)
  m  <- Biobase::exprs(sd)
  for (tissue in unique(pd$Tissue)) {
    ctrl <- pd$Tissue == tissue & pd$Genotype == "WT"
    med  <- apply(m[, ctrl, drop = FALSE], 1, stats::median, na.rm = TRUE)
    med  <- med[is.finite(med)]
    # Within each tissue the WT medians are the reference, so they are ~0 ...
    expect_equal(unname(med), rep(0, length(med)), tolerance = 1e-8)
  }
  # ... whereas a reference pooled across tissues is not.
  pooled <- tAge_preprocessing(eset, species = "mouse", control_group_column = "Genotype",
                               control_group_label = "WT", verbose = FALSE)$scaled_diff
  ctrl_k <- Biobase::pData(pooled)$Tissue == "Kidney" & Biobase::pData(pooled)$Genotype == "WT"
  med_k  <- apply(Biobase::exprs(pooled)[, ctrl_k, drop = FALSE], 1, stats::median, na.rm = TRUE)
  expect_gt(stats::median(abs(med_k), na.rm = TRUE), 0.01)
})

test_that("the species is recorded in the preprocessed objects", {
  proc <- tAge_preprocessing(.tage_example_eset(), species = "mouse", verbose = FALSE)
  for (e in proc) expect_equal(tAge:::.tage_get_species(e), "mouse")
  expect_equal(tAge:::.tage_resolve_species(NULL, proc$scaled_diff), "mouse")
  expect_equal(tAge:::.tage_resolve_species("Mouse", proc$scaled_diff), "mouse")
  expect_error(tAge:::.tage_resolve_species(NULL, .tage_example_eset()), "species")
})

test_that("a control label that matches nothing is reported", {
  eset <- log_transform(filter_genes(.tage_example_eset(), verbose = FALSE), verbose = FALSE)
  expect_warning(control_subtraction(eset, "Genotype", "Mutant", verbose = FALSE),
                 "No sample has Genotype == 'Mutant'")
})

test_that("gene filtering uses the paper's default thresholds", {
  eset <- .tage_example_eset()
  f20 <- filter_genes(eset, verbose = FALSE)
  f25 <- filter_genes(eset, percent_threshold = 25, verbose = FALSE)
  expect_lte(nrow(f25), nrow(f20))
})
