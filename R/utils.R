#' Load h5ad file using reticulate and anndata
#'
#' This function provides an alternative method to load h5ad files using the
#' reticulate package and Python's anndata library, then converts to Seurat.
#'
#' @param h5ad_path Character string specifying the path to the h5ad file.
#' @param python_path Character string specifying the path to Python executable.
#'   Default is NULL (uses current Python environment).
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A Seurat object containing the single-cell data.
#' @export
#' @examples
#' \dontrun{
#' # Using default Python
#' seurat_obj <- load_h5ad_simple("data.h5ad")
#'
#' # Using specific Python environment
#' seurat_obj <- load_h5ad_simple("data.h5ad", python_path = "./.venv/bin/python")
#' }
load_h5ad_simple <- function(h5ad_path, python_path = NULL, verbose = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate package is required. Please install it with: install.packages('reticulate')")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Please install it.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required. Please install it.")
  }

  if (!file.exists(h5ad_path)) {
    stop("File not found: ", h5ad_path)
  }

  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  if (verbose) {
    cat("Loading h5ad file using anndata...\n")
    cat("  File:", h5ad_path, "\n")
  }

  tryCatch({
    ad <- reticulate::import("anndata")

    if (verbose) {
      cat("Reading h5ad file...\n")
    }
    adata <- ad$read_h5ad(h5ad_path)

    if (verbose) {
      cat("  - Cells:", adata$n_obs, "\n")
      cat("  - Genes:", adata$n_vars, "\n")
      cat("Converting to Seurat object...\n")
    }

    counts <- Matrix::t(adata$X)
    rownames(counts) <- adata$var_names$to_list()
    colnames(counts) <- adata$obs_names$to_list()

    metadata <- as.data.frame(adata$obs)

    seurat_obj <- Seurat::CreateSeuratObject(
      counts = counts,
      meta.data = metadata
    )

    if (verbose) {
      cat("✓ Successfully created Seurat object\n")
      cat("  - Cells:", ncol(seurat_obj), "\n")
      cat("  - Genes:", nrow(seurat_obj), "\n")
      cat("  - Metadata columns:", ncol(seurat_obj@meta.data), "\n")
    }

    return(seurat_obj)

  }, error = function(e) {
    stop("Failed to load h5ad file: ", e$message,
         "\nMake sure anndata is installed in your Python environment: pip install anndata")
  })
}


#' Subset Seurat object by metadata criteria
#'
#' Convenience function to subset Seurat objects with multiple criteria.
#'
#' @param seurat_obj A Seurat object to subset.
#' @param filters Named list of filtering criteria. Each element should be a
#'   metadata column name with a vector of values to keep.
#' @param min_cells Minimum number of cells to retain after filtering. If fewer
#'   cells remain, returns NULL with a warning. Default is 1.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A subsetted Seurat object, or NULL if fewer than min_cells remain.
#' @export
#' @examples
#' \dontrun{
#' # Subset by tissue and age
#' seurat_subset <- subset_seurat_by_metadata(
#'   seurat_obj,
#'   filters = list(
#'     tissue = c("Kidney", "Liver"),
#'     age = c("3m", "24m")
#'   )
#' )
#'
#' # Subset by multiple criteria
#' seurat_subset <- subset_seurat_by_metadata(
#'   seurat_obj,
#'   filters = list(
#'     tissue = "Tongue",
#'     sex = "male",
#'     age = "24m"
#'   ),
#'   min_cells = 100
#' )
#' }
subset_seurat_by_metadata <- function(
  seurat_obj,
  filters,
  min_cells = 1,
  verbose = TRUE
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Please install it.")
  }

  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.list(filters) || length(filters) == 0) {
    stop("filters must be a non-empty named list")
  }

  metadata <- seurat_obj@meta.data
  keep_cells <- rep(TRUE, nrow(metadata))

  if (verbose) {
    cat("Subsetting Seurat object...\n")
    cat("  - Starting cells:", ncol(seurat_obj), "\n")
  }

  for (col_name in names(filters)) {
    if (!col_name %in% colnames(metadata)) {
      stop("Column '", col_name, "' not found in metadata")
    }

    values <- filters[[col_name]]
    keep_cells <- keep_cells & (metadata[[col_name]] %in% values)

    if (verbose) {
      cat("  - After filtering by", col_name, ":", sum(keep_cells), "cells\n")
    }
  }

  if (sum(keep_cells) < min_cells) {
    warning("Only ", sum(keep_cells), " cells remain after filtering (min_cells = ", min_cells, "). Returning NULL.")
    return(NULL)
  }

  cell_names <- rownames(metadata)[keep_cells]
  seurat_subset <- subset(seurat_obj, cells = cell_names)

  if (verbose) {
    cat("✓ Subsetting completed\n")
    cat("  - Final cells:", ncol(seurat_subset), "\n")
    cat("  - Retained:", round(ncol(seurat_subset) / ncol(seurat_obj) * 100, 1), "%\n")
  }

  return(seurat_subset)
}


#' Get summary statistics for pseudobulk samples
#'
#' Calculate and display summary statistics for pseudobulk ExpressionSet objects.
#'
#' @param eset An ExpressionSet object created from pseudobulk aggregation.
#' @param group_vars Character vector of metadata columns to group by for summary.
#'   Default is NULL (no grouping).
#' @return A data frame with summary statistics.
#' @export
#' @examples
#' \dontrun{
#' # Basic summary
#' summary_stats <- pseudobulk_summary(eset)
#'
#' # Summary by tissue
#' summary_stats <- pseudobulk_summary(eset, group_vars = "tissue")
#'
#' # Summary by multiple variables
#' summary_stats <- pseudobulk_summary(eset, group_vars = c("tissue", "age"))
#' }
pseudobulk_summary <- function(eset, group_vars = NULL) {
  if (!inherits(eset, "ExpressionSet")) {
    stop("eset must be an ExpressionSet object")
  }

  pdata <- Biobase::pData(eset)
  exprs_data <- Biobase::exprs(eset)

  if (is.null(group_vars)) {
    summary_df <- data.frame(
      n_samples = ncol(eset),
      n_genes = nrow(eset),
      mean_cells_per_sample = mean(pdata$n_cells, na.rm = TRUE),
      median_cells_per_sample = median(pdata$n_cells, na.rm = TRUE),
      min_cells_per_sample = min(pdata$n_cells, na.rm = TRUE),
      max_cells_per_sample = max(pdata$n_cells, na.rm = TRUE),
      mean_counts_per_sample = mean(colSums(exprs_data)),
      median_counts_per_sample = median(colSums(exprs_data))
    )
  } else {
    if (!all(group_vars %in% colnames(pdata))) {
      missing <- setdiff(group_vars, colnames(pdata))
      stop("The following group_vars are not in metadata: ", paste(missing, collapse = ", "))
    }

    if (length(group_vars) == 1) {
      pdata$group <- pdata[[group_vars]]
    } else {
      pdata$group <- apply(pdata[, group_vars, drop = FALSE], 1,
                          function(x) paste(x, collapse = "_"))
    }

    summary_list <- lapply(split(pdata, pdata$group), function(grp_data) {
      grp_samples <- rownames(grp_data)
      grp_exprs <- exprs_data[, grp_samples, drop = FALSE]

      data.frame(
        group = grp_data$group[1],
        n_samples = length(grp_samples),
        mean_cells = mean(grp_data$n_cells, na.rm = TRUE),
        median_cells = median(grp_data$n_cells, na.rm = TRUE),
        mean_counts = mean(colSums(grp_exprs)),
        median_counts = median(colSums(grp_exprs)),
        stringsAsFactors = FALSE
      )
    })

    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL
  }

  return(summary_df)
}
