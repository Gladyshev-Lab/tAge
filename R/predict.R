#' Predict transcriptomic age with one pre-trained model
#'
#' Applies one pre-trained Elastic Net (EN) or Bayesian Ridge (BR) clock to a
#' preprocessed ExpressionSet through the Python model (via reticulate).
#'
#' @param eset An ExpressionSet object containing processed expression data
#'   (one element of the list returned by \code{\link{tAge_preprocessing}}).
#' @param model_path Character string specifying the path to the pre-trained model file.
#' @param species Species of the \emph{samples}: \code{"mouse"}, \code{"rat"},
#'   \code{"human"} or \code{"monkey"} (see \code{\link{tage_species}}). Its
#'   only effect is the rescaling of chronological-age clocks to age units by
#'   the species maximum lifespan; mortality and normalized-age clocks ignore
#'   it. It is \emph{not} the species group the model was trained on
#'   ("Mouse" / "Rodents" / "Multispecies" in \code{\link{list_clocks}}).
#'   Default \code{NULL} takes the species recorded by
#'   \code{\link{tAge_preprocessing}}.
#' @param mode Character string specifying the model type. Must be either "EN" for
#'   Elastic Net or "BR" for Bayesian Ridge.
#' @param return_std Logical. Whether to also return the per-sample predictive
#'   standard deviation, which only Bayesian Ridge models provide. Defaults to
#'   \code{TRUE} for \code{mode = "BR"}. The standard deviations are what
#'   \code{\link{tage_compare_groups}} weights samples by, so keep them if you
#'   intend to run statistics on BR predictions.
#' @param age_units Units of chronological-age predictions: \code{"auto"}
#'   (default; months for rodents, years for primates), \code{"months"} or
#'   \code{"years"}.
#' @param normalized_age Scale of normalized-age clocks: \code{"fraction"} of
#'   the expected maximum lifespan (default) or \code{"percent"} (x100, as in
#'   the paper and the TACO application).
#' @return A data frame containing the predicted transcriptomic age results with
#'   sample information and predicted ages, plus a \code{BR_tAge_std} column
#'   when \code{return_std} is \code{TRUE}. The attribute \code{"tage_units"}
#'   names the unit of the prediction column.
#' @export
predict_tAge_one <- function(eset, model_path, species = NULL, mode,
                             return_std = identical(mode, "BR"),
                             age_units = c("auto", "months", "years"),
                             normalized_age = c("fraction", "percent")) {
  if (missing(model_path) || !file.exists(model_path)) {
    stop("Model path is missing or the file does not exist.")
  }
  if (!(mode %in% c("EN", "BR"))) {
    stop("Mode must be either 'EN' or 'BR'.")
  }
  age_units      <- match.arg(age_units)
  normalized_age <- match.arg(normalized_age)
  species        <- .tage_resolve_species(species, eset)

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

  scaling <- .tage_output_scaling(model_path, species, age_units, normalized_age)

  if (mode == "EN" && isTRUE(return_std)) {
    warning("Elastic net models do not provide predictive standard deviations; ignoring return_std.")
    return_std <- FALSE
  }

  sample_result <- mod$predict_tAge(
    model_path, expr_df, meta_df,
    species = species,
    return_std = isTRUE(return_std),
    prefix = paste0(mode, "_"),
    output_factor = scaling$factor
  )
  sample_result <- as.data.frame(sample_result, check.names = FALSE)
  attr(sample_result, "tage_units") <- stats::setNames(scaling$units, paste0(mode, "_tAge"))
  sample_result
}

# How a clock's raw output is brought onto the reported scale: species maximum
# lifespan (months or years) for chronological clocks, x100 for normalized age
# in percent, untouched otherwise. The outcome comes from the registry or, for
# models outside it, from the file name.
.tage_output_scaling <- function(model_path, species, age_units, normalized_age) {
  outcome <- .clock_outcome(model_path)
  if (is.na(outcome)) {
    warning("Could not tell what ", basename(model_path), " predicts from its name; ",
            "its output is left on the model's native scale.", call. = FALSE)
    return(list(factor = 1, units = "native"))
  }
  switch(
    outcome,
    "Chronological"  = {
      lf <- .tage_lifespan_factor(species, age_units)
      list(factor = lf$factor, units = lf$units)
    },
    "Normalized age" = if (normalized_age == "percent") {
      list(factor = 100, units = "% of maximum lifespan")
    } else {
      list(factor = 1, units = "fraction of maximum lifespan")
    },
    "Mortality"      = list(factor = 1, units = "log10 hazard ratio"),
    list(factor = 1, units = "native")
  )
}

#' Predict transcriptomic age for multiple processed ExpressionSet objects
#'
#' Applies one model per normalisation (\code{scaled_diff}, \code{yugene_diff},
#' ...) and returns the sample metadata with one prediction column per
#' normalisation, named \code{<normalisation>_<mode>_tAge}.
#'
#' @param tAge_eset A named list of ExpressionSet objects, each representing a different
#'   normalization method (e.g., "scaled", "scaled_diff", "yugene", "yugene_diff"),
#'   as returned by \code{\link{tAge_preprocessing}}.
#' @param model_paths A named list of model paths corresponding to each normalization method.
#' @inheritParams predict_tAge_one
#' @param return_std Logical. Whether to keep the per-sample predictive standard
#'   deviation of Bayesian Ridge clocks. Defaults to \code{TRUE} for
#'   \code{mode = "BR"}, adding one \code{<normalisation>_BR_tAge_sd} column per
#'   clock. Pass these to the \code{se_columns} argument of
#'   \code{\link{tage_compare_groups}} to reproduce the reference application's
#'   Bayesian ridge statistics.
#' @return A data frame containing the predicted transcriptomic age results for all
#'   provided ExpressionSet objects, with appropriately named columns. The
#'   attribute \code{"tage_units"} is a named character vector giving the unit of
#'   every prediction column (e.g. \code{"months"}, \code{"log10 hazard ratio"}).
#' @export
predict_tAge <- function(tAge_eset, model_paths, species = NULL, mode,
                         return_std = identical(mode, "BR"),
                         age_units = c("auto", "months", "years"),
                         normalized_age = c("fraction", "percent")) {
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
  age_units      <- match.arg(age_units)
  normalized_age <- match.arg(normalized_age)

  results <- NULL
  units <- character(0)
  for (name in common_names) {
    eset <- tAge_eset[[name]]
    if (!inherits(eset, "ExpressionSet")) {
      stop(paste("Element", name, "in tAge_eset is not an ExpressionSet."))
    }

    model_path <- model_paths[[name]]
    res <- predict_tAge_one(eset, model_path, species, mode, return_std = return_std,
                            age_units = age_units, normalized_age = normalized_age)
    res_units <- attr(res, "tage_units")

    # Ensure a 2D data.frame regardless of how reticulate converts the Python
    # result (some reticulate/pandas versions can return a bare vector for a
    # single column, which breaks the colnames<- below).
    res <- as.data.frame(res, check.names = FALSE)

    # Rename 'EN_tAge' or 'BR_tAge' to name + mode + '_tAge', and the matching
    # predictive standard deviation to '<name>_<mode>_tAge_sd'.
    tAge_col <- paste0(mode, "_tAge")
    new_tAge_col <- paste0(name, "_", mode, "_tAge")
    colnames(res)[colnames(res) == tAge_col] <- new_tAge_col
    units[new_tAge_col] <- unname(res_units[tAge_col])

    new_cols <- new_tAge_col
    std_col <- paste0(mode, "_tAge_std")
    if (std_col %in% colnames(res)) {
      new_sd_col <- paste0(new_tAge_col, "_sd")
      colnames(res)[colnames(res) == std_col] <- new_sd_col
      new_cols <- c(new_cols, new_sd_col)
      units[new_sd_col] <- unname(res_units[tAge_col])
    }

    if (is.null(results)) {
      results <- res
    } else {
      # Add new columns to results, only the prediction columns at a time
      results <- cbind(results, res[, new_cols, drop = FALSE])
    }
  }
  attr(results, "tage_units") <- units
  results
}


#' Run tAge pipeline separately per group factor (e.g., "tissue")
#'
#' Preprocesses every level of \code{split_by} on its own -- gene filtering,
#' normalisation and reference centring all happen within the stratum, as the
#' clocks were trained and applied in the paper -- and predicts on the
#' combined data. This is \code{\link{tAge_preprocessing}} with
#' \code{split_by} followed by \code{\link{predict_tAge}}.
#'
#' @param eset An ExpressionSet, e.g. from pseudobulk aggregation.
#' @param split_by Character. Column in pData to split by (e.g., "tissue").
#' @param model_paths Named list of model paths.
#' @inheritParams tAge_preprocessing
#' @inheritParams predict_tAge_one
#' @param min_samples Integer. Strata with fewer samples are left out, with a
#'   warning. Default 5.
#' @return Data frame with predictions for all strata combined (the
#'   \code{split_by} column is part of the sample metadata).
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
  verbose = TRUE,
  gene_mapping_type = "auto",
  return_std = identical(mode, "BR"),
  age_units = c("auto", "months", "years"),
  normalized_age = c("fraction", "percent")
) {
  if (!split_by %in% colnames(Biobase::pData(eset))) {
    stop(paste0("'", split_by, "' not found in pData"))
  }
  age_units      <- match.arg(age_units)
  normalized_age <- match.arg(normalized_age)

  groups <- Biobase::pData(eset)[[split_by]]
  counts <- table(as.character(groups))
  small  <- names(counts)[counts < min_samples]
  if (length(small)) {
    warning("Leaving out ", split_by, " level(s) with fewer than ", min_samples, " samples: ",
            paste(sprintf("%s (n = %d)", small, counts[small]), collapse = ", "), call. = FALSE)
    eset <- eset[, !(as.character(groups) %in% small)]
  }
  if (ncol(eset) == 0L) stop("No groups with at least ", min_samples, " samples.", call. = FALSE)

  if (verbose) {
    cat("Running tAge by", split_by, "\n")
    cat("  - Groups:", paste(setdiff(names(counts), small), collapse = ", "), "\n\n")
  }

  tAge_all <- tAge_preprocessing(
    eset = eset,
    species = species,
    gene_mapping_type = gene_mapping_type,
    verbose = verbose,
    control_group_column = control_group_column,
    control_group_label = control_group_label,
    count_threshold = count_threshold,
    percent_threshold = percent_threshold,
    split_by = split_by
  )

  combined <- predict_tAge(
    tAge_eset = tAge_all,
    model_paths = model_paths,
    species = species,
    mode = mode,
    return_std = return_std,
    age_units = age_units,
    normalized_age = normalized_age
  )

  if (verbose) {
    cat("\u2713 Combined results:", nrow(combined), "samples from",
        length(setdiff(names(counts), small)), "groups\n")
  }
  combined
}
