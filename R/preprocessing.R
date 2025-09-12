#' Create ExpressionSet object from expression data and phenotype data
#'
#' This function creates an ExpressionSet object from expression data and phenotype data,
#' which is the standard format used in Bioconductor for microarray and RNA-seq data.
#'
#' @param exprs_data A data frame or matrix containing expression data where rows are genes
#'   and columns are samples.
#' @param phenodata A data frame containing phenotype data where rows are samples and
#'   columns are phenotypic variables.
#' @param verbose Logical indicating whether to print progress messages and create
#'   density plots. Default is TRUE.
#' @return An ExpressionSet object containing the expression data and phenotype information.
#' @export
#' @examples
#' # Load example data
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' 
#' # Create ExpressionSet
#' eset <- make_ExpressionSet(expr_data, meta_data)
make_ExpressionSet <- function(exprs_data, phenodata, verbose = TRUE) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package is required. Please install it.")
  }
  
  # Ensure phenodata has the same samples as exprs_data
  phenodata <- phenodata[colnames(exprs_data), , drop = FALSE]
  meta.info <- data.frame(labelDescription = colnames(phenodata))
  
  # Create AnnotatedDataFrame
  pheno <- methods::new("AnnotatedDataFrame", data = phenodata, varMetadata = meta.info)

  # Create ExpressionSet
  eset <- methods::new("ExpressionSet", exprs = as.matrix(exprs_data), phenoData = pheno)

  if (verbose) {
    cat("✓ ExpressionSet created successfully\n")
    cat("  - Number of genes:", nrow(eset), "\n")
    cat("  - Number of samples:", ncol(eset), "\n")
    plot_eset_density(eset, title = "Raw Data")
  }

  return(eset)
}


#' Filter genes based on expression thresholds
#'
#' This function filters genes from an ExpressionSet based on count and percentage thresholds.
#' Genes are retained if they have expression values above the count threshold in at least
#' the specified percentage of samples.
#'
#' @param exprs_set An ExpressionSet object containing expression data.
#' @param count_threshold Numeric threshold for minimum expression count. Default is 10.
#' @param percent_threshold Numeric threshold for minimum percentage of samples that must
#'   have expression above count_threshold. Default is 20.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A filtered ExpressionSet object with only genes passing the thresholds.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Filter genes
#' filtered_eset <- filter_genes(eset, count_threshold = 10, percent_threshold = 20)
filter_genes <- function(exprs_set, count_threshold = 10, percent_threshold = 20, verbose = TRUE) {
  exprs_matrix <- Biobase::exprs(exprs_set)
  
  # Calculate which genes pass the threshold
  genes_passing <- apply(exprs_matrix, 1, function(x) {
    sum(x >= count_threshold, na.rm = TRUE)
  }) >= ncol(exprs_matrix) * (percent_threshold / 100)
    
  # Filter the ExpressionSet
  filtered_set <- exprs_set[genes_passing, ]

  if (verbose) {
    cat("✓ Gene filtering completed\n")
    cat("  - Number of genes before filtering:", nrow(exprs_set), "\n")
    cat("  - Number of genes after filtering:", nrow(filtered_set), "\n")
    cat("  - Percentage of genes retained:", round(nrow(filtered_set) / nrow(exprs_set) * 100, 1), "%\n")
  }

  return(filtered_set)
}




#' Perform RLE (Relative Log Expression) normalization
#'
#' This function performs RLE normalization on an ExpressionSet using the edgeR package.
#' RLE normalization is commonly used for RNA-seq data to correct for library size differences
#' and composition bias.
#'
#' @param original_dataset An ExpressionSet object containing raw expression data.
#' @param verbose Logical indicating whether to print progress messages and create
#'   density plots. Default is FALSE.
#' @return An ExpressionSet object with RLE-normalized expression data.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Perform RLE normalization
#' rle_eset <- RLE_normalization(eset, verbose = TRUE)
RLE_normalization  <- function(original_dataset, verbose=FALSE){
  library(edgeR)
  #First, we need to get library sizes vector
  lib.size = apply(Biobase::exprs(original_dataset), 2, sum)
  
  #Now we should obtain normalization factors for each samples using RLE technique
  edger.rle = lib.size * calcNormFactors(Biobase::exprs(original_dataset), method="RLE")
  
  #Finally we should divide expression values on normalization factors
  RLE_dataset <- original_dataset
  Biobase::exprs(RLE_dataset) = sweep(Biobase::exprs(original_dataset), 2, edger.rle, "/")*10^7

  if (verbose) {
    cat("✓ RLE normalization completed\n")
    plot_eset_density(RLE_dataset, title = "RLE Normalized Data")
  }
  
  return(RLE_dataset)
}


#' Perform YuGene normalization
#'
#' This function performs YuGene normalization on an ExpressionSet. YuGene is a
#' non-parametric normalization method that transforms expression data to a uniform
#' distribution by computing cumulative proportions and applying a specific transformation.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param verbose Logical indicating whether to print progress messages and create
#'   density plots. Default is TRUE.
#' @return An ExpressionSet object with YuGene-normalized expression data.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Perform YuGene normalization
#' yugene_eset <- YuGene(eset, verbose = TRUE)
YuGene <- function(eset, verbose = TRUE) {
  expr_data <- Biobase::exprs(eset)
  
  # Convert to matrix for faster operations
  counts <- as.matrix(expr_data)
  
  # 1. Shift counts to non-negative range
  counts <- sweep(counts, 2, apply(counts, 2, min), "-")
  
  # 2. Clip negative values to 0
  if (any(counts < 0)) {
    warning("Some negative values were set to 0")
    counts[counts < 0] <- 0
  }
  
  # 3. Process each sample
  result <- apply(counts, 2, function(col_data) {
    # Sort in descending order
    sort_idx <- order(col_data, decreasing = TRUE)
    sorted_vals <- col_data[sort_idx]
    
    # Compute cumulative proportions
    cumsum_vals <- cumsum(sorted_vals)
    total <- sum(sorted_vals)
    
    if (total == 0) {
      return(rep(1.0, length(col_data)))
    }
    
    cumprop <- cumsum_vals / total
    
    # Handle duplicates
    for (i in 2:length(cumprop)) {
      if (sorted_vals[i] == sorted_vals[i - 1]) {
        cumprop[i] <- cumprop[i - 1]
      }
    }
    
    # YuGene transform
    final_col <- 1.0 - cumprop
    
    # Return to original order
    sample_result <- numeric(length(col_data))
    sample_result[sort_idx] <- final_col
    
    return(sample_result)
  })
  
  # Convert to data frame with original names
  result_df <- as.data.frame(result)
  rownames(result_df) <- rownames(expr_data)
  colnames(result_df) <- colnames(expr_data)
  
  # Create new ExpressionSet
  eset_yugene <- eset
  Biobase::exprs(eset_yugene) <- as.matrix(result_df)

  if (verbose) {
    cat("✓ YuGene normalization completed\n")
    plot_eset_density(eset_yugene, title = "YuGene Normalized Data")
  }

  return(eset_yugene)
}


#' Apply log10 transformation to expression data
#'
#' This function applies a log10(x + 1) transformation to the expression data in an
#' ExpressionSet. This transformation is commonly used to stabilize variance and
#' make the data more suitable for downstream analysis.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param verbose Logical indicating whether to print progress messages and create
#'   density plots. Default is TRUE.
#' @return An ExpressionSet object with log10-transformed expression data.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Apply log transformation
#' log_eset <- log_transform(eset, verbose = TRUE)
log_transform <- function(eset, verbose = TRUE) {
  eset_log_transformed <- eset
  expr_data <- Biobase::exprs(eset)
  log_data <- log10(expr_data + 1)
  Biobase::exprs(eset_log_transformed) <- log_data
  if (verbose) {
    cat("✓ Log transformation completed\n")
    plot_eset_density(eset_log_transformed, title = "Log-Transformed Data", log_transform = FALSE)
  }
  return(eset_log_transformed)
}


#' Scale expression data to have zero mean and unit variance
#'
#' This function applies z-score scaling to the expression data in an ExpressionSet,
#' where each gene is scaled to have zero mean and unit variance across samples.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param verbose Logical indicating whether to print progress messages and create
#'   density plots. Default is TRUE.
#' @return An ExpressionSet object with scaled expression data.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Scale expression data
#' scaled_eset <- scale_eset(eset, verbose = TRUE)
scale_eset <- function(eset, verbose = TRUE) {
  eset_scaled <- eset
  expr_data <- Biobase::exprs(eset)
  scaled_data <- scale(expr_data)
  Biobase::exprs(eset_scaled) <- scaled_data
  if (verbose) {
    cat("✓ Scaling completed\n")
    plot_eset_density(eset_scaled, title = "Scaled Data", log_transform = FALSE)
  }
  return(eset_scaled)
}


#' Subtract control group median from expression data
#'
#' This function subtracts the median expression of a control group from all samples
#' in the ExpressionSet. If no control samples are found, it uses the overall median
#' across all samples.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param column_name Character string specifying the column name in phenoData that
#'   contains the group labels.
#' @param control_label Character string specifying the label for control samples.
#' @param verbose Logical indicating whether to print progress messages and create
#'   density plots. Default is TRUE.
#' @return An ExpressionSet object with control-subtracted expression data.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Subtract control group (assuming 'Group' column has 'Control' label)
#' control_eset <- control_subtraction(eset, "Group", "Control", verbose = TRUE)
control_subtraction <- function(eset, column_name, control_label, verbose = TRUE) {
  X <- Biobase::exprs(eset)

  control_samples <- Biobase::pData(eset)[[column_name]]
  control_idx <- which(!is.na(control_samples) & control_samples == control_label)

  if (length(control_idx) == 0) {
    Xc <- apply(X, 1, median, na.rm = TRUE)
  } else {
    Xc <- apply(X[, control_idx, drop = FALSE], 1, median, na.rm = TRUE)
  }

  X_adjusted <- sweep(X, 1, Xc, "-")

  eset_final <- eset
  Biobase::exprs(eset_final) <- X_adjusted

  if (verbose) {
    if (length(control_idx) == 0) {
      cat("✓ No control samples found for label '", control_label, "'. Using overall median for subtraction.\n", sep = "")
      plot_eset_density(eset_final, title = "Subtraction using Overall Median")
    } else {
      cat("✓ Control samples found for label '", control_label, "'. Using control group median for subtraction.\n", sep = "")
      plot_eset_density(eset_final, title = paste("Subtraction using Control Group (", control_label, ")", sep = ""))
    }
  }

  return(eset_final)
}

#' Complete preprocessing pipeline for tAge analysis
#'
#' This function performs a complete preprocessing pipeline for transcriptomic age analysis,
#' including gene filtering, normalization, transformation, scaling, control subtraction,
#' and gene ID conversion. It returns multiple versions of the processed data suitable
#' for different analysis approaches.
#'
#' @param eset An ExpressionSet object containing raw expression data.
#' @param gene_list A character vector of gene identifiers to be used for analysis.
#' @param species Character string specifying the species. Options: "mouse", "rat", "human",
#'   "monkey", "rhesus". Default is "mouse".
#' @param from_type Character string specifying the input gene ID type. Options: "auto",
#'   "symbol", "ensembl", "refseq", "uniprot", "entrez". Default is "auto".
#' @param to_type Character string specifying the output gene ID type. Options: "entrez",
#'   "symbol", "ensembl". Default is "entrez".
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @param control_group_column Character string specifying the column name in phenoData
#'   that contains control group labels. Default is NULL.
#' @param control_group_label Character string specifying the label for control samples.
#'   Default is NULL.
#' @param count_threshold Numeric threshold for minimum expression count in gene filtering.
#'   Default is 10.
#' @param percent_threshold Numeric threshold for minimum percentage of samples that must
#'   have expression above count_threshold. Default is 20.
#' @return A list containing six processed ExpressionSet objects:
#'   \item{RLE_normalized}{RLE-normalized data}
#'   \item{log_transformed}{Log-transformed data}
#'   \item{scaled}{Scaled data with gene ID conversion}
#'   \item{scaled_diff}{Scaled data with control subtraction and gene ID conversion}
#'   \item{yugene}{YuGene-normalized data with gene ID conversion}
#'   \item{yugene_diff}{YuGene-normalized data with control subtraction and gene ID conversion}
#' @export
#' @examples
#' # Load example data
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' gene_list <- load_gene_list()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Run complete preprocessing pipeline
#' processed_data <- tAge_preprocessing(eset, gene_list, species = "mouse")
tAge_preprocessing <- function(
  eset,
  gene_list,
  species = "mouse",
  from_type = "auto",
  to_type = "entrez",
  verbose = TRUE,
  control_group_column = NULL,
  control_group_label = NULL,
  count_threshold = 10,
  percent_threshold = 20
) {
  eset_filtered <- filter_genes(eset, count_threshold = count_threshold, percent_threshold = percent_threshold, verbose = verbose)
  eset_RLE  <- RLE_normalization(eset_filtered, verbose)
  eset_log_transformed <- log_transform(eset_RLE, verbose = verbose)
  eset_scaled <- scale_eset(eset_log_transformed, verbose = verbose)
  eset_scaled_diff <- control_subtraction(eset_scaled, column_name = control_group_column, control_label = control_group_label, verbose = verbose)
  eset_yugene <- YuGene(eset_scaled, verbose = verbose)
  eset_yugene_diff <- control_subtraction(eset_yugene, column_name = control_group_column, control_label = control_group_label, verbose = verbose)
  
  eset_scaled_gene_ids <- process_genes(eset_scaled, gene_list, species, from_type, to_type, verbose)
  eset_scaled_diff_gene_ids <- process_genes(eset_scaled_diff, gene_list, species, from_type, to_type, verbose)
  eset_yugene_gene_ids <- process_genes(eset_yugene, gene_list, species, from_type, to_type, verbose)
  eset_yugene_diff_gene_ids <- process_genes(eset_yugene_diff, gene_list, species, from_type, to_type, verbose)

  return(list(
    RLE_normalized = eset_RLE,
    log_transformed = eset_log_transformed,
    scaled = eset_scaled_gene_ids,
    scaled_diff = eset_scaled_diff_gene_ids,
    yugene = eset_yugene_gene_ids,
    yugene_diff = eset_yugene_diff_gene_ids
  ))
}
