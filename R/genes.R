#' Map genes in an ExpressionSet using local CSV mapping tables
#'
#' @param eset An ExpressionSet object.
#' @param species One of "human", "mouse", "rat", "monkey".
#' @param gene_mapping_type One of "Ensembl" or "Gene.Symbol".
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return ExpressionSet with Entrez (mouse) gene IDs as rownames.
#' @export
map_genes <- function(eset,
                      species,
                      gene_mapping_type,
                      verbose = TRUE) {
  metadata_dir <- get_metadata_dir()

  valid_gene_types <- c("Ensembl", "Gene.Symbol")
  if (!gene_mapping_type %in% valid_gene_types) {
    stop(sprintf("gene_mapping_type must be one of {%s}, got '%s'",
                 paste(valid_gene_types, collapse = ", "), gene_mapping_type))
  }

  gene_map <- .load_gene_mapping(metadata_dir, species, gene_mapping_type)
  eset     <- .apply_gene_mapping(eset, gene_map)

  if (!species %in% c("mouse", "monkey")) {
    eset <- .apply_ortholog_mapping(eset, metadata_dir, species, verbose)
  }

  eset <- .set_var_names(eset, species)

  return(eset)
}


.load_gene_mapping <- function(metadata_dir, species, gene_mapping_type) {
  gene_table_path <- file.path(metadata_dir, sprintf("Gene_table_%s.csv", species))

  if (!file.exists(gene_table_path)) {
    stop(sprintf("Gene table not found at %s", gene_table_path))
  }

  gene_table <- read.csv(gene_table_path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!gene_mapping_type %in% colnames(gene_table)) {
    stop(sprintf("'%s' not found in %s", gene_mapping_type, gene_table_path))
  }

  if (species == "monkey") {
    return(.create_monkey_gene_mapping(metadata_dir, gene_table, gene_mapping_type))
  }

  gene_table <- gene_table[!is.na(gene_table[["Entrez"]]), ]

  # Deduplicate by gene_mapping_type, keeping first occurrence
  gene_table <- gene_table[!duplicated(gene_table[[gene_mapping_type]]), ]

  setNames(gene_table[["Entrez"]], gene_table[[gene_mapping_type]])
}


.create_monkey_gene_mapping <- function(metadata_dir, gene_table, gene_mapping_type) {
  if (gene_mapping_type == "Ensembl") {
    monkey_ens_map <- setNames(gene_table[["Ensembl"]], gene_table[["Ensembl"]])
  } else {
    monkey_ens_map <- setNames(gene_table[["Ensembl"]], gene_table[["Gene.Symbol"]])
  }

  orthologs_path <- file.path(metadata_dir, "Orthologs_monkey_to_mouse_5.0.csv")
  orthologs      <- read.csv(orthologs_path, stringsAsFactors = FALSE, check.names = FALSE)
  orthologs      <- orthologs[!is.na(orthologs[["Entrez.mouse"]]), ]
  ortholog_map   <- setNames(orthologs[["Entrez.mouse"]], orthologs[["Ensembl.macaca"]])

  # Chain: original ID -> Ensembl -> mouse Entrez
  gene_map <- list()
  for (orig_name in names(monkey_ens_map)) {
    mouse_entrez <- ortholog_map[[ monkey_ens_map[[orig_name]] ]]
    if (!is.null(mouse_entrez) && !is.na(mouse_entrez)) {
      gene_map[[orig_name]] <- mouse_entrez
    }
  }

  unlist(gene_map)
}


.apply_gene_mapping <- function(eset, gene_map) {
  genes  <- rownames(eset)
  mapped <- gene_map[genes]
  valid  <- !is.na(mapped)

  expr_valid <- Biobase::exprs(eset)[valid, , drop = FALSE]
  mapped_ids <- unname(mapped[valid])

  # Aggregate duplicates by summing
  expr_agg <- apply(expr_valid, 2, function(x) {
    tapply(x, mapped_ids, sum)
  })

  # Track which original genes mapped to each Entrez ID
  original_genes_map <- tapply(
    names(mapped[valid]),
    mapped_ids,
    paste, collapse = ";"
  )

  agg_ids <- rownames(expr_agg)

  fdata <- data.frame(
    mapped_genes   = agg_ids,
    original_genes = unname(original_genes_map[agg_ids]),
    row.names      = agg_ids,
    stringsAsFactors = FALSE
  )

  new_eset <- Biobase::ExpressionSet(
    assayData   = as.matrix(expr_agg),
    phenoData   = Biobase::AnnotatedDataFrame(Biobase::pData(eset)),
    featureData = Biobase::AnnotatedDataFrame(fdata)
  )

  return(new_eset)
}

.apply_ortholog_mapping <- function(eset, metadata_dir, species, verbose) {
  orthologs_path <- file.path(metadata_dir, "Table_of_orthologs.csv")

  if (!file.exists(orthologs_path)) {
    stop(sprintf("Ortholog table not found at %s", orthologs_path))
  }

  ortholog_key <- paste0("Entrez.", .capitalize(species))
  orthologs    <- read.csv(orthologs_path, stringsAsFactors = FALSE, check.names = FALSE)
  orthologs    <- orthologs[!is.na(orthologs[["Entrez.Mouse"]]), c(ortholog_key, "Entrez.Mouse")]

  # Deduplicate by source Entrez
  orthologs    <- orthologs[!duplicated(orthologs[[ortholog_key]]), ]
  ortholog_map <- setNames(orthologs[["Entrez.Mouse"]], orthologs[[ortholog_key]])

  mapped_entrez <- as.character(Biobase::fData(eset)[["mapped_genes"]])
  mouse_ids     <- unname(ortholog_map[mapped_entrez])

  ortho_valid <- !is.na(mouse_ids)
  ortholog_unmapped <- sum(!ortho_valid)

  if (verbose) {
    if (ortholog_unmapped > 0) {
      message(sprintf("Warning: %d orthologs could not be mapped and will be dropped",
                      ortholog_unmapped))
    }
    message(sprintf("Mapped %d genes to mouse orthologs", sum(ortho_valid)))
  }

  # Filter valid, then drop duplicate Mouse Entrez
  expr_valid  <- Biobase::exprs(eset)[ortho_valid, , drop = FALSE]
  mouse_valid <- mouse_ids[ortho_valid]

  dup_mouse <- duplicated(mouse_valid)
  expr_valid  <- expr_valid[!dup_mouse, , drop = FALSE]
  mouse_valid <- mouse_valid[!dup_mouse]

  # Assign Mouse Entrez as rownames
  rownames(expr_valid) <- as.character(mouse_valid)

  fdata <- data.frame(
    mapped_genes   = as.character(mouse_valid),
    ortholog_genes = as.character(mouse_valid),
    row.names      = as.character(mouse_valid),
    stringsAsFactors = FALSE
  )

  Biobase::ExpressionSet(
    assayData   = as.matrix(expr_valid),
    phenoData   = Biobase::AnnotatedDataFrame(Biobase::pData(eset)),
    featureData = Biobase::AnnotatedDataFrame(fdata)
  )
}


.set_var_names <- function(eset, species) {
  index_column <- if (species %in% c("mouse", "monkey")) {
    "mapped_genes"
  } else {
    "ortholog_genes"
  }

  rownames(eset) <- as.character(as.integer(Biobase::fData(eset)[[index_column]]))
  return(eset)
}


# Capitalizes first letter, lowercases the rest (e.g. "human" -> "Human")
.capitalize <- function(s) {
  paste0(toupper(substring(s, 1, 1)), tolower(substring(s, 2)))
}