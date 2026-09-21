#' Map genes in an ExpressionSet using local CSV mapping tables
#'
#' Translates the row names of \code{eset} to mouse Entrez IDs, the gene space
#' every clock operates in: identifiers of the given species are looked up in
#' the bundled gene table, and rat, human and macaque genes are then carried
#' over to mouse through 1:1 ortholog tables. Genes that do not map are
#' dropped; identifiers that collapse onto one Entrez ID are summed.
#'
#' @param eset An ExpressionSet object.
#' @param species One of "human", "mouse", "rat", "monkey".
#' @param gene_mapping_type Identifier type of the row names: "Ensembl",
#'   "Gene.Symbol", "Entrez", or "auto" (default) to detect it as the type
#'   with the most matches in the gene table. Ensembl IDs may carry a version
#'   suffix (\code{ENSMUSG00000000001.4}); it is stripped when that is what
#'   makes them match.
#' @param verbose Logical. Print progress messages. Default TRUE.
#' @return ExpressionSet with Entrez (mouse) gene IDs as rownames.
#' @export
map_genes <- function(eset,
                      species,
                      gene_mapping_type = "auto",
                      verbose = TRUE) {
  metadata_dir <- get_metadata_dir()
  species <- .tage_species_row(species)$species

  valid_gene_types <- c("auto", "Ensembl", "Gene.Symbol", "Entrez")
  if (!gene_mapping_type %in% valid_gene_types) {
    stop(sprintf("gene_mapping_type must be one of {%s}, got '%s'",
                 paste(valid_gene_types, collapse = ", "), gene_mapping_type))
  }

  gene_table <- .read_gene_table(metadata_dir, species)
  detected   <- .detect_gene_mapping_type(rownames(eset), gene_table, gene_mapping_type)
  gene_mapping_type <- detected$type
  if (verbose) {
    cat(sprintf("\u2713 Gene identifiers: %s (%d of %d row names found in the %s gene table%s)\n",
                gene_mapping_type, detected$n_matched, nrow(eset), species,
                if (detected$strip_version) ", after removing Ensembl version suffixes" else ""))
  }

  gene_map <- .load_gene_mapping(metadata_dir, species, gene_mapping_type, gene_table)
  eset     <- .apply_gene_mapping(eset, gene_map, keys = detected$keys)

  if (!species %in% c("mouse", "monkey")) {
    eset <- .apply_ortholog_mapping(eset, metadata_dir, species, verbose)
  }

  eset <- .set_var_names(eset, species)

  return(eset)
}


.read_gene_table <- function(metadata_dir, species) {
  gene_table_path <- file.path(metadata_dir, sprintf("Gene_table_%s.csv", species))
  if (!file.exists(gene_table_path)) {
    stop(sprintf("Gene table not found at %s", gene_table_path))
  }
  gene_table <- utils::read.csv(gene_table_path, stringsAsFactors = FALSE, check.names = FALSE)
  # Entrez IDs are read as numbers; the row names they are matched against
  # are text.
  gene_table[["Entrez"]] <- ifelse(is.na(gene_table[["Entrez"]]), NA_character_,
                                   as.character(as.integer(gene_table[["Entrez"]])))
  gene_table
}

.strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

# Which identifier type the row names are, as the type with the most matches
# in the gene table. Returns the type, the
# keys to look up (row names, or Ensembl IDs with their version suffix
# removed when that is what makes them match) and the match count.
.detect_gene_mapping_type <- function(genes, gene_table, requested = "auto") {
  genes <- as.character(genes)
  candidates <- intersect(c("Ensembl", "Gene.Symbol", "Entrez"), colnames(gene_table))
  if (requested != "auto") {
    if (!requested %in% colnames(gene_table)) {
      stop(sprintf("'%s' identifiers are not available in the gene table of this species", requested),
           call. = FALSE)
    }
    candidates <- requested
  }

  count <- function(keys, column) sum(keys %in% unique(gene_table[[column]]))
  matches <- vapply(candidates, function(col) count(genes, col), numeric(1))
  stripped <- .strip_ensembl_version(genes)
  matches_stripped <- if ("Ensembl" %in% candidates) count(stripped, "Ensembl") else 0

  best <- candidates[which.max(matches)]
  strip_version <- "Ensembl" %in% candidates && matches_stripped > max(matches)
  if (strip_version) best <- "Ensembl"
  n_matched <- if (strip_version) matches_stripped else max(matches)

  if (n_matched == 0) {
    detail <- paste(sprintf("%s: %d", candidates, matches[candidates]), collapse = ", ")
    stop("None of the row names match the gene table (matches by identifier type -- ", detail,
         "). Check `species` and the identifiers, e.g. ", paste(head(genes, 3), collapse = ", "),
         call. = FALSE)
  }

  list(type = best, keys = if (strip_version) stripped else genes,
       n_matched = as.integer(n_matched), strip_version = strip_version)
}


.load_gene_mapping <- function(metadata_dir, species, gene_mapping_type, gene_table = NULL) {
  if (is.null(gene_table)) gene_table <- .read_gene_table(metadata_dir, species)

  if (!gene_mapping_type %in% colnames(gene_table)) {
    stop(sprintf("'%s' identifiers are not available for %s", gene_mapping_type, species))
  }

  if (species == "monkey") {
    return(.create_monkey_gene_mapping(metadata_dir, gene_table, gene_mapping_type))
  }

  gene_table <- gene_table[!is.na(gene_table[["Entrez"]]), ]

  # Deduplicate by gene_mapping_type, keeping first occurrence
  gene_table <- gene_table[!is.na(gene_table[[gene_mapping_type]]), ]
  gene_table <- gene_table[!duplicated(gene_table[[gene_mapping_type]]), ]

  setNames(as.numeric(gene_table[["Entrez"]]), gene_table[[gene_mapping_type]])
}


.create_monkey_gene_mapping <- function(metadata_dir, gene_table, gene_mapping_type) {
  if (gene_mapping_type == "Ensembl") {
    monkey_ens_map <- setNames(gene_table[["Ensembl"]], gene_table[["Ensembl"]])
  } else {
    # Gene.Symbol or (macaque) Entrez -> Ensembl
    monkey_ens_map <- setNames(gene_table[["Ensembl"]], gene_table[[gene_mapping_type]])
  }

  orthologs_path <- file.path(metadata_dir, "Orthologs_monkey_to_mouse_5.0.csv")
  orthologs      <- read.csv(orthologs_path, stringsAsFactors = FALSE, check.names = FALSE)
  orthologs      <- orthologs[!is.na(orthologs[["Entrez.mouse"]]), ]
  ortholog_map   <- setNames(orthologs[["Entrez.mouse"]], orthologs[["Ensembl.macaca"]])

  # Chain: original ID -> Ensembl -> mouse Entrez. `[` returns NA for macaque
  # genes without a mouse ortholog; those are dropped.
  ens   <- as.character(unname(monkey_ens_map))
  mouse <- unname(ortholog_map[ens])
  keep  <- !is.na(names(monkey_ens_map)) & !is.na(ens) & !is.na(mouse)

  stats::setNames(mouse[keep], names(monkey_ens_map)[keep])
}


# `keys` are the identifiers looked up in gene_map (the row names, or Ensembl
# IDs without version suffix); the row names themselves are what is recorded
# as the original genes.
.apply_gene_mapping <- function(eset, gene_map, keys = rownames(eset)) {
  genes  <- rownames(eset)
  mapped <- gene_map[as.character(keys)]
  names(mapped) <- genes
  valid  <- !is.na(mapped)
  if (!any(valid)) stop("No gene could be mapped to Entrez IDs.", call. = FALSE)

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