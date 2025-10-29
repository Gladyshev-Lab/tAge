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


predict_tAge_multiple <- function(tAge_eset, model_paths, species, mode) {
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
    res <- predict_tAge(eset, model_path, species, mode)
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
