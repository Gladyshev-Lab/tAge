#' Aggregate single-cell data into pseudobulk samples based on read coverage
#'
#' This function aggregates single-cell RNA-seq data into pseudobulk samples
#' by sequentially accumulating cells until a cumulative read coverage threshold
#' is reached. Each pseudobulk sample is guaranteed to contain at least the
#' specified number of total reads, ensuring sufficient sequencing depth for
#' downstream transcriptomic age prediction.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param coverage_threshold Integer specifying the minimum cumulative read count
#'   per pseudobulk sample. Cells are sequentially added to a group until this
#'   threshold is met, then a new group begins. Default is 1e6 (1 million reads).
#' @param assay Character string specifying which assay to use. Default is "RNA".
#' @param layer Character string specifying which layer to use. Default is "counts".
#' @param shuffle Logical indicating whether to randomly shuffle cells before
#'   aggregation. Shuffling breaks any ordering present in the data. Default is FALSE.
#' @param seed Integer seed for reproducibility when shuffle = TRUE. Default is NULL.
#' @param new_sample_prefix Character string prefix for pseudobulk sample names.
#'   Default is "".
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return An ExpressionSet object with pseudobulk expression data and metadata
#'   including cumulative_coverage (total reads) and n_cells per pseudobulk sample.
#' @export
#' @examples
#' \dontrun{
#' # Basic coverage-based aggregation
#' eset <- aggregate_pseudobulk(seurat_obj, coverage_threshold = 1e6)
#'
#' # With shuffling for randomized grouping
#' eset <- aggregate_pseudobulk(seurat_obj, coverage_threshold = 5e5,
#'                              shuffle = TRUE, seed = 42)
#' }
aggregate_pseudobulk <- function(
  seurat_obj,
  coverage_threshold = 1e6,
  assay = "RNA",
  layer = "counts",
  shuffle = FALSE,
  seed = NULL,
  new_sample_prefix = "",
  verbose = TRUE
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required.")
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object")
  if (coverage_threshold <= 0) stop("coverage_threshold must be positive")

  counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  n_cells <- ncol(counts)
  n_genes <- nrow(counts)

  if (verbose) {
    cat("Aggregating pseudobulk samples by coverage...\n")
    cat("  - Total cells:", n_cells, "\n")
    cat("  - Total genes:", n_genes, "\n")
    cat("  - Coverage threshold:", format(coverage_threshold, big.mark = ","), "reads\n")
  }

  cell_totals <- Matrix::colSums(counts)

  # Cell ordering
  if (shuffle) {
    if (!is.null(seed)) set.seed(seed)
    indices <- sample.int(n_cells)
  } else {
    indices <- seq_len(n_cells)
  }

  # Assign groups by cumulative coverage
  ordered_totals <- cell_totals[indices]
  cum_totals <- cumsum(ordered_totals)
  group_assignments <- as.integer(cum_totals %/% coverage_threshold) + 1L

  # Last group may not reach threshold — that's ok
  n_groups <- max(group_assignments)

  # Build sparse indicator matrix (cells x groups) for fast aggregation
  indicator <- Matrix::sparseMatrix(
    i = seq_len(n_cells),
    j = group_assignments,
    x = 1,
    dims = c(n_cells, n_groups)
  )

  # Single matrix multiply: (genes x cells) %*% (cells x groups) = (genes x groups)
  pseudobulk_counts <- counts[, indices] %*% indicator

  # Compute group stats
  group_cell_counts <- as.integer(Matrix::colSums(indicator))
  group_coverages <- as.numeric(Matrix::colSums(
    Matrix::Diagonal(n_cells, ordered_totals) %*% indicator
  ))

  # Convert to regular matrix
  pseudobulk_counts <- as.matrix(pseudobulk_counts)

  sample_names <- paste0(new_sample_prefix, "pseudobulk_", seq_len(n_groups))
  colnames(pseudobulk_counts) <- sample_names
  rownames(pseudobulk_counts) <- rownames(counts)

  pseudobulk_meta <- data.frame(
    sample_id = sample_names,
    cumulative_coverage = group_coverages,
    n_cells = group_cell_counts,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    n_full <- sum(group_coverages >= coverage_threshold)
    cat("  - Pseudobulk samples created:", n_groups, "\n")
    cat("  - Samples meeting threshold:", n_full, "\n")
    if (n_groups > n_full) {
      cat("  - Remainder sample (below threshold): 1 (",
          format(group_coverages[n_groups], big.mark = ","), " reads, ",
          group_cell_counts[n_groups], " cells)\n")
    }
    cat("  - Coverage per sample (median):",
        format(median(group_coverages), big.mark = ","), "reads\n")
    cat("  - Cells per sample (median):", median(group_cell_counts), "\n")
  }

  eset <- make_ExpressionSet(
    exprs_data = as.data.frame(pseudobulk_counts),
    phenodata = pseudobulk_meta,
    verbose = verbose
  )

  if (verbose) cat("\u2713 Pseudobulk aggregation completed\n")

  return(eset)
}


#' Aggregate single-cell data into pseudobulk samples within obs column groups
#'
#' This function first splits the data by all unique combinations of the specified
#' observation (metadata) columns, then applies coverage-based aggregation within
#' each group. This is useful for creating pseudobulk samples stratified by
#' biological variables such as sample, cell type, tissue, or condition.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param obs_column_names Character vector of metadata column names to stratify by.
#'   For example, c("sample_id", "cell_type") will first split cells into groups
#'   defined by each unique sample_id x cell_type combination, then aggregate
#'   within each group.
#' @param coverage_threshold Integer specifying the minimum cumulative read count
#'   per pseudobulk sample. Default is 1e6.
#' @param assay Character string specifying which assay to use. Default is "RNA".
#' @param layer Character string specifying which layer to use. Default is "counts".
#' @param shuffle Logical indicating whether to randomly shuffle cells before
#'   aggregation within each group. Default is FALSE.
#' @param seed Integer seed for reproducibility when shuffle = TRUE. Default is NULL.
#' @param new_sample_prefix Character string prefix for pseudobulk sample names.
#'   Default is "".
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return An ExpressionSet object with pseudobulk expression data. Metadata includes
#'   the stratification columns, cumulative_coverage, and n_cells.
#' @export
#' @examples
#' \dontrun{
#' # Aggregate within each sample
#' eset <- aggregate_on_obs_columns(seurat_obj,
#'                                  obs_column_names = "sample_id",
#'                                  coverage_threshold = 1e6)
#'
#' # Aggregate within each sample x cell_type combination
#' eset <- aggregate_on_obs_columns(seurat_obj,
#'                                  obs_column_names = c("sample_id", "cell_type"),
#'                                  coverage_threshold = 5e5)
#' }
aggregate_on_obs_columns <- function(
  seurat_obj,
  obs_column_names,
  coverage_threshold = 1e6,
  assay = "RNA",
  layer = "counts",
  shuffle = FALSE,
  seed = NULL,
  new_sample_prefix = "",
  verbose = TRUE
) {
  if (!all(obs_column_names %in% colnames(seurat_obj@meta.data))) {
    missing <- setdiff(obs_column_names, colnames(seurat_obj@meta.data))
    stop("Columns not found in metadata: ", paste(missing, collapse = ", "))
  }

  counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  metadata <- seurat_obj@meta.data

  # Create group key
  group_key <- apply(metadata[, obs_column_names, drop = FALSE], 1, paste, collapse = "___")
  unique_groups <- unique(group_key)

  if (verbose) {
    cat("Aggregating pseudobulk by obs columns:",
        paste(obs_column_names, collapse = ", "), "\n")
    cat("  - Unique combinations:", length(unique_groups), "\n")
  }

  all_exprs <- list()
  all_meta <- list()
  global_counter <- 0L

  for (grp in unique_groups) {
    cell_idx <- which(group_key == grp)
    if (length(cell_idx) < 1) next

    grp_counts <- counts[, cell_idx, drop = FALSE]
    grp_meta_row <- metadata[cell_idx[1], obs_column_names, drop = FALSE]

    # Coverage-based aggregation on this sub-matrix
    result <- .aggregate_matrix(
      counts_matrix = grp_counts,
      coverage_threshold = coverage_threshold,
      shuffle = shuffle,
      seed = seed
    )

    if (is.null(result)) next

    n <- ncol(result$counts)
    new_names <- paste0(new_sample_prefix, "pseudobulk_", seq(global_counter + 1, global_counter + n))
    global_counter <- global_counter + n

    colnames(result$counts) <- new_names

    sample_meta <- data.frame(
      sample_id = new_names,
      cumulative_coverage = result$coverages,
      n_cells = result$cell_counts,
      row.names = new_names,
      stringsAsFactors = FALSE
    )
    for (col in obs_column_names) {
      sample_meta[[col]] <- grp_meta_row[[col]]
    }

    all_exprs[[length(all_exprs) + 1]] <- result$counts
    all_meta[[length(all_meta) + 1]] <- sample_meta
  }

  if (length(all_exprs) == 0) stop("No valid pseudobulk samples created")

  combined_exprs <- do.call(cbind, all_exprs)
  combined_meta <- do.call(rbind, all_meta)

  eset <- make_ExpressionSet(
    exprs_data = as.data.frame(combined_exprs),
    phenodata = combined_meta,
    verbose = TRUE
  )

  if (verbose) {
    cat("✓ Combined pseudobulk aggregation completed\n")
    cat("  - Total pseudobulk samples:", ncol(eset), "\n")
    cat("  - Total genes:", nrow(eset), "\n")
  }

  return(eset)
}


# Internal: aggregate a count matrix by coverage threshold
.aggregate_matrix <- function(counts_matrix, coverage_threshold, shuffle, seed) {
  n_cells <- ncol(counts_matrix)
  n_genes <- nrow(counts_matrix)

  if (is(counts_matrix, "dgCMatrix") || is(counts_matrix, "sparseMatrix")) {
    cell_totals <- Matrix::colSums(counts_matrix)
  } else {
    cell_totals <- colSums(counts_matrix)
  }

  if (shuffle) {
    if (!is.null(seed)) set.seed(seed)
    indices <- sample.int(n_cells)
  } else {
    indices <- seq_len(n_cells)
  }

  ordered_totals <- cell_totals[indices]

  group_assignments <- integer(n_cells)
  group_coverages <- numeric()
  group_cell_counts <- integer()
  current_group <- 1L
  current_coverage <- 0
  current_n_cells <- 0L

  for (i in seq_len(n_cells)) {
    current_coverage <- current_coverage + ordered_totals[i]
    current_n_cells <- current_n_cells + 1L
    group_assignments[i] <- current_group

    if (current_coverage >= coverage_threshold) {
      group_coverages <- c(group_coverages, current_coverage)
      group_cell_counts <- c(group_cell_counts, current_n_cells)
      current_group <- current_group + 1L
      current_coverage <- 0
      current_n_cells <- 0L
    }
  }

  if (current_n_cells > 0) {
    group_coverages <- c(group_coverages, current_coverage)
    group_cell_counts <- c(group_cell_counts, current_n_cells)
  }

  n_groups <- length(group_coverages)
  if (n_groups == 0) return(NULL)

  pseudobulk <- matrix(0, nrow = n_genes, ncol = n_groups)
  rownames(pseudobulk) <- rownames(counts_matrix)

  for (g in seq_len(n_groups)) {
    member_positions <- which(group_assignments == g)
    member_indices <- indices[member_positions]
    grp <- counts_matrix[, member_indices, drop = FALSE]

    if (is(grp, "dgCMatrix") || is(grp, "sparseMatrix")) {
      pseudobulk[, g] <- Matrix::rowSums(grp)
    } else {
      pseudobulk[, g] <- rowSums(grp)
    }
  }

  list(counts = pseudobulk, coverages = group_coverages, cell_counts = group_cell_counts)
}

#' Load AnnData h5ad file and convert to Seurat object
#'
#' This function loads a single-cell dataset in h5ad format (AnnData) and
#' converts it to a Seurat object for downstream analysis.
#'
#' @param h5ad_path Character string specifying the path to the h5ad file.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @return A Seurat object containing the single-cell data.
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- load_h5ad_to_seurat("data.h5ad")
#' }
load_h5ad_to_seurat <- function(h5ad_path, verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Please install it.")
  }
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    stop("SeuratDisk package is required. Please install it with: remotes::install_github('mojaveazure/seurat-disk')")
  }

  if (!file.exists(h5ad_path)) {
    stop("File not found: ", h5ad_path)
  }

  if (verbose) {
    cat("Loading h5ad file:", h5ad_path, "\n")
  }

  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path)

  if (verbose) {
    cat("Converting h5ad to h5seurat format...\n")
  }

  SeuratDisk::Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)

  if (verbose) {
    cat("Loading Seurat object...\n")
  }

  seurat_obj <- SeuratDisk::LoadH5Seurat(h5seurat_path)

  if (file.exists(h5seurat_path)) {
    file.remove(h5seurat_path)
  }

  if (verbose) {
    cat("\u2713 Successfully loaded Seurat object\n")
    cat("  - Cells:", ncol(seurat_obj), "\n")
    cat("  - Genes:", nrow(seurat_obj), "\n")
  }

  return(seurat_obj)
}
