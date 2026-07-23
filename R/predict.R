#' Predict transcriptomic age using pre-trained models
#'
#' This function predicts transcriptomic age using pre-trained Elastic Net (EN) or
#' Bayesian Ridge (BR) models. It interfaces with Python models through the reticulate
#' package to perform the predictions.
#'
#' @param eset An ExpressionSet object containing processed expression data.
#' @param model_path Character string specifying the path to the pre-trained model file.
#'   The file must exist and be compatible with the specified mode.
#' @param species Character string specifying the species for the model. This should
#'   match the species the model was trained on.
#' @param mode Character string specifying the model type. Must be either "EN" for
#'   Elastic Net or "BR" for Bayesian Ridge.
#' @return A data frame containing the predicted transcriptomic age results with
#'   sample information and predicted ages.
#' @export
predict_tAge_one <- function(eset, model_path, species, mode) {
  if (missing(model_path) || !file.exists(model_path)) {
    stop("Model path is missing or the file does not exist.")
  }
  if (!(mode %in% c("EN", "BR"))) {
    stop("Mode must be either 'EN' or 'BR'.")
  }

  # Check if EN or BR in the model_path
  if (mode == "EN" && !grepl("EN_", basename(model_path))) {
    warning("The model path does not seem to correspond to an 'EN' model.")
  }
  if (mode == "BR" && !grepl("BR_", basename(model_path))) {
    warning("The model path does not seem to correspond to a 'BR' model.")
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required.")
  }
  mod <- reticulate::import_from_path(
    "tage_predict",
    path = system.file("python", package = "tAge"),
    convert = TRUE
  )

  expr_df <- Biobase::exprs(eset)
  meta_df <- Biobase::pData(eset)
  
  if (is.matrix(expr_df))  expr_df <- as.data.frame(expr_df, check.names = FALSE)
  if (is.matrix(meta_df))  meta_df <- as.data.frame(meta_df, check.names = FALSE)

  if (mode == "EN") {
    sample_result <- mod$predict_tAge(model_path, expr_df, meta_df, species = species, return_std = FALSE, prefix = "EN_")
  } else if (mode == "BR") {
    sample_result <- mod$predict_tAge(model_path, expr_df, meta_df, species = species, return_std = FALSE, prefix = "BR_")
  } else {
    stop("Unsupported mode. Use 'EN' or 'BR'.")
  }
  return(sample_result)
}

#' Predict transcriptomic age for multiple processed ExpressionSet objects
#'
#' This function predicts transcriptomic age for multiple processed ExpressionSet objects
#' using pre-trained models. It supports different normalization methods and model types.
#' @param tAge_eset A named list of ExpressionSet objects, each representing a different
#'   normalization method (e.g., "scaled", "scaled_diff", "yugene", "yugene_diff").
#' @param model_paths A named list of model paths corresponding to each normalization method.
#' @param species Character string specifying the species for the models.
#' @param mode Character string specifying the model type. Must be either "EN" for
#'   Elastic Net or "BR" for Bayesian Ridge.
#' @return A data frame containing the predicted transcriptomic age results for all
#'   provided ExpressionSet objects, with appropriately named columns.
#' @export
predict_tAge <- function(tAge_eset, model_paths, species, mode) {
  if (!is.list(tAge_eset) || length(tAge_eset) == 0) {
    stop("tAge_eset must be a non-empty list of ExpressionSet objects.")
  }
  valid_names <- c("scaled", "scaled_diff", "yugene", "yugene_diff")
  tAge_eset <- tAge_eset[names(tAge_eset) %in% valid_names]
  if (length(tAge_eset) == 0) {
    stop("No valid ExpressionSet objects found in tAge_eset. Valid names are: 'scaled', 'scaled_diff', 'yugene', 'yugene_diff'.")
  }
  if (!is.list(model_paths) || length(model_paths) == 0) {
    stop("model_paths must be a non-empty named list of model paths.")
  }
  model_paths <- model_paths[names(model_paths) %in% valid_names]
  if (length(model_paths) == 0) {
    stop("No valid model paths found in model_paths. Valid names are: 'scaled', 'scaled_diff', 'yugene', 'yugene_diff'.")
  }
  # Use only common names between tAge_eset and model_paths
  common_names <- intersect(names(tAge_eset), names(model_paths))
  if (length(common_names) == 0) {
    stop("No overlapping names between tAge_eset and model_paths. Ensure at least one shared name like 'scaled_diff'.")
  }

  results <- NULL
  for (name in common_names) {
    eset <- tAge_eset[[name]]
    if (!inherits(eset, "ExpressionSet")) {
      stop(paste("Element", name, "in tAge_eset is not an ExpressionSet."))
    }
  
    model_path <- model_paths[[name]]
    res <- predict_tAge_one(eset, model_path, species, mode)

    # Ensure a 2D data.frame regardless of how reticulate converts the Python
    # result (some reticulate/pandas versions can return a bare vector for a
    # single column, which breaks the colnames<- below).
    res <- as.data.frame(res, check.names = FALSE)

    # Rename 'EN_tAge' or 'BR_tAge' to name + mode + '_tAge'
    tAge_col <- if (mode == "EN") "EN_tAge" else "BR_tAge"
    new_tAge_col <- paste0(name, "_", mode, "_tAge")
    colnames(res)[colnames(res) == tAge_col] <- new_tAge_col
    if (is.null(results)) {
      results <- res
    } else {
      # Add new columns to results, only the tAge column at a time
      results <- cbind(results, res[, new_tAge_col, drop = FALSE])
    }
  }
  return(results)
}


#' Run tAge pipeline separately per group factor (e.g., "tissue")
#'
#' @param eset An ExpressionSet from pseudobulk aggregation.
#' @param split_by Character. Column in pData to split by (e.g., "tissue").
#' @param model_paths Named list of model paths.
#' @param species Character. Species name.
#' @param mode Character. "EN" or "BR".
#' @param control_group_column Character or NULL. Column for control subtraction.
#' @param control_group_label Character or NULL. Label for control group.
#' @param count_threshold Integer. Gene filtering threshold.
#' @param percent_threshold Numeric. Gene filtering percent.
#' @param min_samples Integer. Minimum samples per group to run tAge. Default 5.
#' @param verbose Logical.
#' @return Data frame with predictions from all groups combined.
#' @export
tAge_by_group <- function(
  eset,
  split_by,
  model_paths,
  species = "mouse",
  mode = "EN",
  control_group_column = NULL,
  control_group_label = NULL,
  count_threshold = 10,
  percent_threshold = 20,
  min_samples = 5,
  verbose = TRUE
) {
  if (!split_by %in% colnames(Biobase::pData(eset))) {
    stop(paste0("'", split_by, "' not found in pData"))
  }

  groups <- unique(Biobase::pData(eset)[[split_by]])

  if (verbose) {
    cat("Running tAge by", split_by, "\n")
    cat("  - Groups:", length(groups), "\n")
    cat("  - Groups:", paste(groups, collapse = ", "), "\n\n")
  }

  results_list <- list()

  for (grp in groups) {
    mask <- Biobase::pData(eset)[[split_by]] == grp
    eset_sub <- eset[, mask]

    if (ncol(eset_sub) < min_samples) {
      if (verbose) cat("Skipping", grp, "(", ncol(eset_sub), "samples < min", min_samples, ")\n\n")
      next
    }

    if (verbose) cat("=== ", grp, " (", ncol(eset_sub), " samples) ===\n")

    tryCatch({
      tAge_sub <- tAge_preprocessing(
        eset = eset_sub,
        species = species,
        verbose = verbose,
        control_group_column = control_group_column,
        control_group_label = control_group_label,
        count_threshold = count_threshold,
        percent_threshold = percent_threshold
      )

      res <- predict_tAge(
        tAge_eset = tAge_sub,
        model_paths = model_paths,
        species = species,
        mode = mode
      )

      res[[split_by]] <- grp
      results_list[[grp]] <- res

      if (verbose) cat("\n")

    }, error = function(e) {
      if (verbose) cat("Error for", grp, ":", e$message, "\n\n")
    })
  }

  if (length(results_list) == 0) {
    stop("No groups produced valid predictions")
  }

  # Combine — align columns across groups
  all_cols <- unique(unlist(lapply(results_list, colnames)))
  results_aligned <- lapply(results_list, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    for (col in missing_cols) df[[col]] <- NA
    df[, all_cols, drop = FALSE]
  })

  combined <- do.call(rbind, results_aligned)

  if (verbose) {
    cat("\u2713 Combined results:", nrow(combined), "samples from",
        length(results_list), "groups\n")
  }

  return(combined)
}