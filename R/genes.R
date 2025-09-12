#' Convert gene identifiers between different formats
#'
#' This function converts gene identifiers between different formats using the biomaRt
#' package. It supports conversion between gene symbols, Ensembl IDs, RefSeq IDs,
#' UniProt IDs, and Entrez IDs for multiple species.
#'
#' @param gene_ids A character vector of gene identifiers to convert.
#' @param species Character string specifying the species. Options: "mouse", "rat", "human",
#'   "monkey", "rhesus". Default is "mouse".
#' @param from_type Character string specifying the input gene ID type. Options: "auto",
#'   "symbol", "ensembl", "refseq", "uniprot", "entrez". Default is "auto".
#' @param to_type Character string specifying the output gene ID type. Options: "entrez",
#'   "symbol", "ensembl". Default is "entrez".
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A data frame with columns:
#'   \item{original_id}{Original gene identifiers}
#'   \item{converted_id}{Converted gene identifiers}
#'   \item{gene_name}{Gene symbols}
#'   \item{biotype}{Gene biotype}
#' @export
#' @examples
#' # Convert mouse gene symbols to Entrez IDs
#' gene_symbols <- c("Actb", "Gapdh", "Tbp")
#' converted <- convert_gene_ids(gene_symbols, species = "mouse", 
#'                              from_type = "symbol", to_type = "entrez")
convert_gene_ids <- function(gene_ids, species = "mouse", 
                           from_type = "auto", to_type = "entrez",
                           verbose = TRUE) {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt package is required for gene ID conversion. Please install it.")
  }
  
  # Validate inputs
  if (length(gene_ids) == 0) {
    stop("gene_ids cannot be empty")
  }
  
  # Map species to Ensembl dataset names
  species_datasets <- list(
    "mouse" = "mmusculus_gene_ensembl",
    "rat" = "rnorvegicus_gene_ensembl", 
    "human" = "hsapiens_gene_ensembl",
    "monkey" = "mmulatta_gene_ensembl",
    "rhesus" = "mmulatta_gene_ensembl"
  )
  
  if (!species %in% names(species_datasets)) {
    stop("Unsupported species. Use: mouse, rat, human, monkey, or rhesus")
  }
  
  dataset <- species_datasets[[species]]
  
  # Auto-detect input type if requested
  if (from_type == "auto") {
    from_type <- detect_gene_id_type(gene_ids[1:min(10, length(gene_ids))])
    if (verbose) {
      cat("Auto-detected input type:", from_type, "\n")
    }
  }
  
  # Map input type to biomaRt attributes
  input_attributes <- get_input_attributes(from_type)
  
  # Map output type to biomaRt attributes  
  output_attributes <- get_output_attributes(to_type)
  
  if (verbose) {
    cat("Converting", length(gene_ids), "genes from", from_type, "to", to_type, "for", species, "\n")
  }
  
  tryCatch({
    # Connect to Ensembl
    mart <- biomaRt::useMart("ensembl", dataset = dataset)
    
    # Get gene information
    gene_info <- biomaRt::getBM(
      attributes = c(input_attributes, output_attributes, "external_gene_name", "gene_biotype"),
      filters = input_attributes,
      values = gene_ids,
      mart = mart
    )
    
    if (nrow(gene_info) == 0) {
      warning("No genes found in database")
      return(data.frame(
        original_id = gene_ids,
        converted_id = NA,
        gene_name = NA,
        biotype = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    # Clean up the results
    colnames(gene_info)[colnames(gene_info) == input_attributes] <- "original_id"
    colnames(gene_info)[colnames(gene_info) == output_attributes] <- "converted_id"
    
    # Remove duplicates and keep best matches
    gene_info <- gene_info[!duplicated(gene_info$original_id), ]
    
    # Add missing genes
    missing_genes <- setdiff(gene_ids, gene_info$original_id)
    if (length(missing_genes) > 0) {
      missing_df <- data.frame(
        original_id = missing_genes,
        converted_id = NA,
        external_gene_name = NA,
        gene_biotype = NA,
        stringsAsFactors = FALSE
      )
      gene_info <- rbind(gene_info, missing_df)
    }
    
    # Reorder to match input order
    gene_info <- gene_info[match(gene_ids, gene_info$original_id), ]
    
    # Rename columns for clarity
    colnames(gene_info)[colnames(gene_info) == "external_gene_name"] <- "gene_name"
    
    if (verbose) {
      cat("✓ Conversion completed\n")
      cat("  - Successfully converted:", sum(!is.na(gene_info$converted_id)), "genes\n")
      cat("  - Failed to convert:", sum(is.na(gene_info$converted_id)), "genes\n")
    }
    
    return(gene_info)
    
  }, error = function(e) {
    warning("Gene ID conversion failed: ", e$message)
    return(data.frame(
      original_id = gene_ids,
      converted_id = NA,
      gene_name = NA,
      biotype = NA,
      stringsAsFactors = FALSE
    ))
  })
}


#' Automatically detect gene identifier type
#'
#' This function automatically detects the type of gene identifiers based on their format.
#' It recognizes Ensembl IDs, RefSeq IDs, UniProt IDs, Entrez IDs, and gene symbols.
#'
#' @param gene_ids A character vector of gene identifiers to analyze.
#' @return A character string indicating the detected ID type: "ensembl", "refseq",
#'   "uniprot", "entrez", or "symbol".
#' @export
#' @examples
#' # Detect ID type for different gene identifiers
#' detect_gene_id_type(c("ENSG00000000003", "ENSG00000000005"))  # "ensembl"
#' detect_gene_id_type(c("NM_000001", "NM_000002"))              # "refseq"
#' detect_gene_id_type(c("12345", "67890"))                      # "entrez"
#' detect_gene_id_type(c("ACTB", "GAPDH"))                       # "symbol"
detect_gene_id_type <- function(gene_ids) {
  # Remove NA and empty values
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  
  if (length(gene_ids) == 0) return("symbol")
  
  # Check for Ensembl IDs (start with ENS)
  if (all(grepl("^ENS", gene_ids))) {
    return("ensembl")
  }
  
  # Check for RefSeq IDs (start with NM_, NR_, XM_, XR_)
  if (all(grepl("^[NX][MR]_", gene_ids))) {
    return("refseq")
  }
  
  # Check for UniProt IDs (start with P_, Q_, O_, etc.)
  if (all(grepl("^[A-Z][0-9][A-Z0-9]{5}$", gene_ids))) {
    return("uniprot")
  }
  
  # Check for Entrez IDs (all numeric)
  if (all(grepl("^[0-9]+$", gene_ids))) {
    return("entrez")
  }
  
  # Default to gene symbols
  return("symbol")
}


#' Get biomaRt input attributes for gene ID conversion
#'
#' This is a helper function that maps input gene ID types to the corresponding
#' biomaRt attribute names used for querying the Ensembl database.
#'
#' @param from_type Character string specifying the input gene ID type.
#' @return A character string with the corresponding biomaRt attribute name.
#' @export
#' @examples
#' get_input_attributes("symbol")    # "external_gene_name"
#' get_input_attributes("ensembl")   # "ensembl_gene_id"
#' get_input_attributes("entrez")    # "entrezgene_id"
get_input_attributes <- function(from_type) {
  switch(from_type,
    "symbol" = "external_gene_name",
    "ensembl" = "ensembl_gene_id", 
    "refseq" = "refseq_mrna",
    "uniprot" = "uniprotswissprot",
    "entrez" = "entrezgene_id",
    stop("Unsupported input type: ", from_type)
  )
}


#' Get biomaRt output attributes for gene ID conversion
#'
#' This is a helper function that maps output gene ID types to the corresponding
#' biomaRt attribute names used for querying the Ensembl database.
#'
#' @param to_type Character string specifying the output gene ID type.
#' @return A character string with the corresponding biomaRt attribute name.
#' @export
#' @examples
#' get_output_attributes("entrez")   # "entrezgene_id"
#' get_output_attributes("symbol")   # "external_gene_name"
#' get_output_attributes("ensembl")  # "ensembl_gene_id"
get_output_attributes <- function(to_type) {
  switch(to_type,
    "entrez" = "entrezgene_id",
    "symbol" = "external_gene_name",
    "ensembl" = "ensembl_gene_id",
    stop("Unsupported output type: ", to_type)
  )
}


#' Convert mouse gene symbols to NCBI Entrez IDs
#'
#' This is a convenience function that converts mouse gene symbols to NCBI Entrez IDs
#' using the convert_gene_ids function with species set to "mouse".
#'
#' @param gene_symbols A character vector of mouse gene symbols to convert.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A data frame with converted gene information (see convert_gene_ids for details).
#' @export
#' @examples
#' # Convert mouse gene symbols to Entrez IDs
#' mouse_genes <- c("Actb", "Gapdh", "Tbp", "Hprt1")
#' converted <- mouse_symbols_to_ncbi(mouse_genes)
mouse_symbols_to_ncbi <- function(gene_symbols, verbose = TRUE) {
  convert_gene_ids(gene_symbols, species = "mouse", from_type = "symbol", 
                  to_type = "entrez", verbose = verbose)
}


#' Convert human gene symbols to NCBI Entrez IDs
#'
#' This is a convenience function that converts human gene symbols to NCBI Entrez IDs
#' using the convert_gene_ids function with species set to "human".
#'
#' @param gene_symbols A character vector of human gene symbols to convert.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A data frame with converted gene information (see convert_gene_ids for details).
#' @export
#' @examples
#' # Convert human gene symbols to Entrez IDs
#' human_genes <- c("ACTB", "GAPDH", "TBP", "HPRT1")
#' converted <- human_symbols_to_ncbi(human_genes)
human_symbols_to_ncbi <- function(gene_symbols, verbose = TRUE) {
  convert_gene_ids(gene_symbols, species = "human", from_type = "symbol", 
                  to_type = "entrez", verbose = verbose)
}


#' Convert rat gene symbols to NCBI Entrez IDs
#'
#' This is a convenience function that converts rat gene symbols to NCBI Entrez IDs
#' using the convert_gene_ids function with species set to "rat".
#'
#' @param gene_symbols A character vector of rat gene symbols to convert.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A data frame with converted gene information (see convert_gene_ids for details).
#' @export
#' @examples
#' # Convert rat gene symbols to Entrez IDs
#' rat_genes <- c("Actb", "Gapdh", "Tbp", "Hprt1")
#' converted <- rat_symbols_to_ncbi(rat_genes)
rat_symbols_to_ncbi <- function(gene_symbols, verbose = TRUE) {
  convert_gene_ids(gene_symbols, species = "rat", from_type = "symbol", 
                  to_type = "entrez", verbose = verbose)
}


#' Convert monkey gene symbols to NCBI Entrez IDs
#'
#' This is a convenience function that converts monkey gene symbols to NCBI Entrez IDs
#' using the convert_gene_ids function with species set to "monkey".
#'
#' @param gene_symbols A character vector of monkey gene symbols to convert.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A data frame with converted gene information (see convert_gene_ids for details).
#' @export
#' @examples
#' # Convert monkey gene symbols to Entrez IDs
#' monkey_genes <- c("ACTB", "GAPDH", "TBP", "HPRT1")
#' converted <- monkey_symbols_to_ncbi(monkey_genes)
monkey_symbols_to_ncbi <- function(gene_symbols, verbose = TRUE) {
  convert_gene_ids(gene_symbols, species = "monkey", from_type = "symbol", 
                  to_type = "entrez", verbose = verbose)
}


#' Convert gene IDs for multiple species in batch
#'
#' This function converts gene identifiers for multiple species in a single call.
#' It processes each species separately and returns a list with results for each species.
#'
#' @param gene_ids A list of character vectors, where each element contains gene
#'   identifiers for one species. Names should correspond to species names.
#' @param species A character vector of species names corresponding to each element
#'   in gene_ids. Must have the same length as gene_ids.
#' @param from_type Character string specifying the input gene ID type. Options: "auto",
#'   "symbol", "ensembl", "refseq", "uniprot", "entrez". Default is "auto".
#' @param to_type Character string specifying the output gene ID type. Options: "entrez",
#'   "symbol", "ensembl". Default is "entrez".
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A list with converted gene information for each species (see convert_gene_ids
#'   for details of individual results).
#' @export
#' @examples
#' # Convert genes for multiple species
#' gene_lists <- list(
#'   mouse = c("Actb", "Gapdh"),
#'   human = c("ACTB", "GAPDH")
#' )
#' species <- c("mouse", "human")
#' converted <- batch_convert_gene_ids(gene_lists, species)
batch_convert_gene_ids <- function(gene_ids, species, from_type = "auto", 
                                 to_type = "entrez", verbose = TRUE) {
  
  if (length(gene_ids) != length(species)) {
    stop("gene_ids and species must have the same length")
  }
  
  results <- list()
  
  for (i in seq_along(gene_ids)) {
    if (verbose) {
      cat("\n--- Converting", species[i], "genes ---\n")
    }
    
    results[[species[i]]] <- convert_gene_ids(
      gene_ids[[i]], 
      species = species[i], 
      from_type = from_type, 
      to_type = to_type, 
      verbose = verbose
    )
  }
  
  return(results)
}

#' Process genes in ExpressionSet for tAge analysis
#'
#' This function processes genes in an ExpressionSet by converting gene IDs, filtering
#' to common genes, collapsing duplicates, and reordering to match a reference gene list.
#' It is designed to prepare data for transcriptomic age analysis.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param gene_list A character vector of reference gene identifiers to match against.
#' @param species Character string specifying the species. Options: "mouse", "rat", "human",
#'   "monkey", "rhesus". Default is "mouse".
#' @param from_type Character string specifying the input gene ID type. Options: "auto",
#'   "symbol", "ensembl", "refseq", "uniprot", "entrez". Default is "auto".
#' @param to_type Character string specifying the output gene ID type. Options: "entrez",
#'   "symbol", "ensembl". Default is "entrez".
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return An ExpressionSet object with processed genes matching the reference gene list.
#'   Missing genes from the reference list are added with NA values.
#' @export
#' @examples
#' # Load example data
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' gene_list <- load_gene_list()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Process genes for tAge analysis
#' processed_eset <- process_genes(eset, gene_list, species = "mouse")
process_genes <- function(eset, gene_list, species = "mouse", 
                          from_type = "auto", to_type = "entrez", 
                          verbose = TRUE) {
  if (!inherits(eset, "ExpressionSet")) {
    stop("eset must be an ExpressionSet object")
  }

  if (length(gene_list) == 0) {
    stop("gene_list cannot be empty")
  }

  if (verbose) {
    cat("Converting gene IDs...\n")
  }

  # 1) Convert gene IDs
  conv <- convert_gene_ids(rownames(eset), species = species, 
                           from_type = from_type, to_type = to_type, 
                           verbose = verbose)
  ids  <- conv$converted_id

  if (verbose) {
    cat("  - Available genes:", sum(!is.na(ids)), "\n")
    cat("  - Coverage:", round(sum(!is.na(ids)) / length(gene_list) * 100, 1), "%\n")
    cat("  - First genes in the list:", paste(head(gene_list), collapse = ", "), "\n")
  }

  # 2) Select common genes
  common_genes <- intersect(gene_list, ids)
  ids[!ids %in% common_genes] <- NA
  if (verbose) {
    cat("  - Common genes in dataset:", length(common_genes), "\n")
    cat("  - Coverage of gene_list:", round(length(common_genes) / length(gene_list) * 100, 1), "%\n")
  }

  # 3) Drop rows with no mapping
  keep <- !is.na(ids)
  eset2 <- eset[keep, ]
  ids2  <- ids[keep]
  
  # 4) Collapse duplicate rows that share the same gene id
  #    (mean; you can use max, median, or custom if you prefer)
  E       <- Biobase::exprs(eset2)
  E_coll  <- avereps(E, ids2, fun = mean)   # rownames(E_coll) are unique gene IDs
  
  # 5) Carry forward feature data for one representative per gene (first occurrence)
  fvars      <- Biobase::fData(eset2)
  rep_idx    <- match(rownames(E_coll), ids2)
  fvars_coll <- fvars[rep_idx, , drop = FALSE]
  rownames(fvars_coll) <- rownames(E_coll)

  # 6) Rebuild an ExpressionSet with collapsed matrix
  eset_coll <- Biobase::ExpressionSet(
    assayData   = E_coll,
    phenoData   = Biobase::phenoData(eset2),
    featureData = methods::new("AnnotatedDataFrame", data = fvars_coll),
    experimentData = Biobase::experimentData(eset2),
    annotation  = Biobase::annotation(eset2),
  )

  # 7) Add gene IDs from gene_list if missing
  missing_genes <- setdiff(gene_list, rownames(eset_coll))
  if (length(missing_genes) > 0) {
    if (verbose) {
      cat("  - Adding", length(missing_genes), "missing genes with NA expressions\n")
    }
    n_samples <- ncol(eset_coll)

    ## (a) NA-matrix
    na_matrix <- matrix(NA_real_, nrow = length(missing_genes), ncol = n_samples,
                        dimnames = list(missing_genes, colnames(eset_coll)))

    storage.mode(na_matrix) <- storage.mode(Biobase::exprs(eset_coll))

    ## (b) featureData
    fvars_template <- Biobase::fData(eset_coll)
    if (ncol(fvars_template) == 0) {
      na_fvars <- data.frame(row.names = missing_genes)
    } else {
      na_fvars <- as.data.frame(
        setNames(
          replicate(ncol(fvars_template), NA, simplify = FALSE),
          colnames(fvars_template)
        ),
        stringsAsFactors = FALSE
      )
      rownames(na_fvars) <- missing_genes
      for (j in seq_len(ncol(na_fvars))) {
        mode(na_fvars[[j]]) <- mode(fvars_template[[j]])
      }
    }

    na_eset <- Biobase::ExpressionSet(
      assayData      = na_matrix,
      phenoData      = Biobase::phenoData(eset_coll),
      featureData    = methods::new("AnnotatedDataFrame", data = na_fvars),
      experimentData = Biobase::experimentData(eset_coll),
      annotation     = Biobase::annotation(eset_coll)
    )

    eset_coll <- Biobase::combine(eset_coll, na_eset)
  } else {
    if (verbose) cat("  - All genes from gene_list are present in the dataset\n")
  }

  # 8) Reorder rows to match gene_list
  final_order <- match(gene_list, rownames(eset_coll))
  eset_coll   <- eset_coll[final_order, ]

  # 9) Make sure gene IDs are integers
  rownames(eset_coll) <- as.character(as.integer(rownames(eset_coll)))

  return(eset_coll)
}

#' Map gene IDs between species using orthologs
#'
#' This function maps gene identifiers between different species using ortholog
#' relationships from the Ensembl database. It uses biomaRt to find orthologs
#' between the source and target species.
#'
#' @param gene_ids A character vector of gene identifiers to map.
#' @param species_from Character string specifying the source species dataset name
#'   (e.g., "hsapiens_gene_ensembl" for human). Default is "hsapiens_gene_ensembl".
#' @param species_to Character string specifying the target species dataset name
#'   (e.g., "mmusculus_gene_ensembl" for mouse). Default is "mmusculus_gene_ensembl".
#' @param id_type Character string specifying the gene ID type to use for mapping.
#'   Default is "entrezgene_id".
#' @return A data frame with ortholog mapping information:
#'   \item{From_ID}{Original gene identifiers}
#'   \item{From_Name}{Gene names in source species}
#'   \item{To_ID}{Ortholog gene identifiers in target species}
#'   \item{To_Name}{Gene names in target species}
#' @export
#' @examples
#' # Map human genes to mouse orthologs
#' human_genes <- c("12345", "67890", "11111")
#' orthologs <- map_to_orthologs(human_genes, 
#'                              species_from = "hsapiens_gene_ensembl",
#'                              species_to = "mmusculus_gene_ensembl")
map_to_orthologs <- function(gene_ids, species_from = "hsapiens_gene_ensembl", 
                            species_to = "mmusculus_gene_ensembl", 
                            id_type = "entrezgene_id") {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt package is required for ortholog mapping. Please install it.")
  }
  
  tryCatch({
    # Connect to source species
    source_mart <- biomaRt::useMart("ensembl", dataset = species_from)
    
    # Connect to target species
    target_mart <- biomaRt::useMart("ensembl", dataset = species_to)
    
    # Get ortholog mapping
    mapping <- biomaRt::getLDS(
      attributes = c(id_type, "external_gene_name"),
      filters = id_type,
      values = gene_ids,
      mart = source_mart,
      attributesL = c(id_type, "external_gene_name"),
      martL = target_mart
    )
    
    colnames(mapping) <- c("From_ID", "From_Name", "To_ID", "To_Name")
    return(mapping)
    
  }, error = function(e) {
    warning("Ortholog mapping failed: ", e$message)
    return(data.frame(From_ID = gene_ids, From_Name = NA, To_ID = NA, To_Name = NA))
  })
}
