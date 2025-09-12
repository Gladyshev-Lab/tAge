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
#' @examples
#' # Load example data and preprocess
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' gene_list <- load_gene_list()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' processed_data <- tAge_preprocessing(eset, gene_list, species = "mouse")
#' 
#' # Predict using Elastic Net model (example with dummy model path)
#' # results <- predict_tAge(processed_data$scaled, 
#' #                        model_path = "path/to/model.pkl", 
#' #                        species = "mouse", mode = "EN")
predict_tAge <- function(eset, model_path, species, mode) {
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

  data <- Biobase::exprs(eset)
  metadata <- Biobase::pData(eset)
  
  if (is.matrix(data))   data    <- as.data.frame(data, check.names = FALSE)
  if (is.matrix(metadata)) metadata <- as.data.frame(metadata, check.names = FALSE)

  if (mode == "EN") {
    sample_result <- mod$predict_tAge(model_path, data, metadata, species = species, return_std = FALSE, prefix = "EN_")
  } else if (mode == "BR") {
    sample_result <- mod$predict_tAge(model_path, data, metadata, species = species, return_std = FALSE, prefix = "BR_")
  } else {
    stop("Unsupported mode. Use 'EN' or 'BR'.")
  }
  return(sample_result)
}


#' Predict transcriptomic age for multiple processed datasets
#'
#' This function predicts transcriptomic age for multiple processed ExpressionSet
#' objects (e.g., from tAge_preprocessing). It processes each dataset separately
#' and combines the results into a single data frame with predictions from all
#' normalization methods.
#'
#' @param tAge_esets_list A list of ExpressionSet objects from tAge_preprocessing.
#'   Must contain valid names: "scaled", "scaled_diff", "yugene", "yugene_diff".
#' @param model_path Character string specifying the path to the pre-trained model file.
#'   The file must exist and be compatible with the specified mode.
#' @param species Character string specifying the species for the model. This should
#'   match the species the model was trained on.
#' @param mode Character string specifying the model type. Must be either "EN" for
#'   Elastic Net or "BR" for Bayesian Ridge.
#' @return A data frame containing predicted transcriptomic age results from all
#'   processed datasets, with columns named according to the dataset and mode
#'   (e.g., "scaled_EN_tAge", "yugene_diff_BR_tAge").
#' @export
#' @examples
#' # Load example data and preprocess
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' gene_list <- load_gene_list()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' processed_data <- tAge_preprocessing(eset, gene_list, species = "mouse")
#' 
#' # Predict using all processed datasets (example with dummy model path)
#' # all_results <- predict_tAge_all(processed_data, 
#' #                                model_path = "path/to/model.pkl", 
#' #                                species = "mouse", mode = "EN")
predict_tAge_all <- function(tAge_esets_list, model_path, species, mode) {
  if (!is.list(tAge_esets_list) || length(tAge_esets_list) == 0) {
    stop("tAge_esets_list must be a non-empty list of ExpressionSet objects.")
  }

  valid_names <- c("scaled", "scaled_diff", "yugene", "yugene_diff")
  tAge_esets_list <- tAge_esets_list[names(tAge_esets_list) %in% valid_names]
  if (length(tAge_esets_list) == 0) {
    stop("No valid ExpressionSet objects found in tAge_esets_list. Valid names are: 'scaled', 'scaled_diff', 'yugene', 'yugene_diff'.")
  }
  
  results <- NULL
  for (name in names(tAge_esets_list)) {
    eset <- tAge_esets_list[[name]]
    if (!inherits(eset, "ExpressionSet")) {
      stop(paste("Element", name, "is not an ExpressionSet."))
    }
    res <- predict_tAge(eset, model_path, species, mode)
  
    # Rename 'EN_tAge' or 'BR_tAge' to name + mode + '_tAge'
    tAge_col <- if (mode == "EN") "EN_tAge" else "BR_tAge"
    new_tAge_col <- paste0(name, "_", mode, "_tAge")
    colnames(res)[colnames(res) == tAge_col] <- new_tAge_col

    if (is.null(results)) {
      results <- res
    } else {
      # Add new columns to results, only one tAge column at a time
      results <- cbind(results, res[, new_tAge_col, drop = FALSE])
    }

  }

  return(results)
}
