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
#'   the paper).
#' @return A data frame containing the predicted transcriptomic age results with
#'   sample information and predicted ages, plus a \code{BR_tAge_std} column
#'   when \code{return_std} is \code{TRUE}. The attribute \code{"tage_units"}
#'   names the unit of the prediction column.
#' @export
predict_tAge_one <- function(eset, model_path, species = NULL, mode,
                             return_std = identical(mode, "BR"),
                             age_units = c("auto", "months", "years"),
                             normalized_age = c("fraction", "percent")) {
  if (missing(model_path) || length(model_path) != 1L) {
    stop("`model_path` must be the path of one model; use predict_tAge() for several.",
         call. = FALSE)
  }
  if (!file.exists(model_path)) stop("Model file does not exist: ", model_path, call. = FALSE)
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

#' Predict transcriptomic age with several clocks
#'
#' Applies clocks to the preprocessed representations returned by
#' \code{\link{tAge_preprocessing}} and returns the sample metadata with one
#' prediction column per clock.
#'
#' \code{model_paths} is either
#' \itemize{
#'   \item a clock table from \code{\link{list_clocks}} or
#'   \code{\link{list_module_clocks}} with a \code{path} column (e.g. from
#'   \code{\link{download_clocks}}): each row is applied to the representation
#'   its \code{scaling} names (\code{Scaled} -> \code{scaled_diff},
#'   \code{YuGene} -> \code{yugene_diff}), and the column is named after the
#'   model file, as in the Python package; or
#'   \item a named list, representation -> model path(s). With one path per
#'   representation the column is \code{<representation>_<mode>_tAge}; with
#'   several, each column is named after its model file.
#' }
#'
#' @param tAge_eset A named list of ExpressionSet objects, each representing a different
#'   normalization method (e.g., "scaled", "scaled_diff", "yugene", "yugene_diff"),
#'   as returned by \code{\link{tAge_preprocessing}}.
#' @param model_paths A clock table with a \code{path} column, or a named list
#'   of model paths per representation; see Details.
#' @param mode \code{"EN"} or \code{"BR"}. Default \code{NULL} takes each
#'   model's type from the clock table, the registry or its file name.
#' @inheritParams predict_tAge_one
#' @param return_std Logical. Whether to keep the per-sample predictive standard
#'   deviation of Bayesian Ridge clocks. Default \code{NULL} keeps it for every
#'   Bayesian ridge clock, as a \code{<column>_sd} column. Pass these to the
#'   \code{se_columns} argument of \code{\link{tage_compare_groups}} for the
#'   Bayesian ridge statistics.
#' @return A data frame containing the predicted transcriptomic age results for all
#'   provided ExpressionSet objects, with appropriately named columns. The
#'   attribute \code{"tage_units"} is a named character vector giving the unit of
#'   every prediction column (e.g. \code{"months"}, \code{"log10 hazard ratio"}).
#' @examples
#' \dontrun{
#' clocks <- download_clocks(list_clocks(type = "EN", species = "Rodents",
#'                                       tissue = "Multi-Tissue"))
#' res <- predict_tAge(tAge_eset, clocks)       # six columns, named by model file
#' attr(res, "tage_units")
#' }
#' @export
predict_tAge <- function(tAge_eset, model_paths, species = NULL, mode = NULL,
                         return_std = NULL,
                         age_units = c("auto", "months", "years"),
                         normalized_age = c("fraction", "percent")) {
  if (!is.list(tAge_eset) || length(tAge_eset) == 0) {
    stop("tAge_eset must be a non-empty list of ExpressionSet objects.")
  }
  age_units      <- match.arg(age_units)
  normalized_age <- match.arg(normalized_age)
  if (!is.null(mode) && !mode %in% c("EN", "BR")) stop("Mode must be either 'EN' or 'BR'.")

  jobs <- .tage_prediction_jobs(model_paths, names(tAge_eset), mode)

  results <- NULL
  units <- character(0)
  for (i in seq_len(nrow(jobs))) {
    eset <- tAge_eset[[jobs$eset[i]]]
    if (!inherits(eset, "ExpressionSet")) {
      stop("Element ", jobs$eset[i], " in tAge_eset is not an ExpressionSet.", call. = FALSE)
    }
    m <- jobs$mode[i]
    # An explicit return_std = TRUE reaches the elastic net models too, which
    # say that they have no standard deviation.
    keep_std <- if (is.null(return_std)) m == "BR" else isTRUE(return_std)
    res <- predict_tAge_one(eset, jobs$path[i], species, m, return_std = keep_std,
                            age_units = age_units, normalized_age = normalized_age)
    res_units <- attr(res, "tage_units")

    # Ensure a 2D data.frame regardless of how reticulate converts the Python
    # result (some reticulate/pandas versions can return a bare vector for a
    # single column, which breaks the colnames<- below).
    res <- as.data.frame(res, check.names = FALSE)

    tAge_col <- paste0(m, "_tAge")
    column <- jobs$column[i]
    colnames(res)[colnames(res) == tAge_col] <- column
    units[column] <- unname(res_units[tAge_col])
    new_cols <- column

    std_col <- paste0(m, "_tAge_std")
    if (std_col %in% colnames(res)) {
      sd_col <- paste0(column, "_sd")
      colnames(res)[colnames(res) == std_col] <- sd_col
      new_cols <- c(new_cols, sd_col)
      units[sd_col] <- unname(res_units[tAge_col])
    }

    results <- if (is.null(results)) res else cbind(results, res[, new_cols, drop = FALSE])
  }
  attr(results, "tage_units") <- units
  results
}

# One row per prediction: which representation, which model, which type and
# which output column.
.tage_prediction_jobs <- function(model_paths, eset_names, mode) {
  valid_names <- c("scaled", "scaled_diff", "yugene", "yugene_diff")

  if (is.data.frame(model_paths)) {
    if (!"path" %in% names(model_paths)) {
      stop("The clock table needs a `path` column; add it with download_clocks() or ",
           "clocks$path <- file.path(dir, clocks$filename).", call. = FALSE)
    }
    if (!"scaling" %in% names(model_paths)) {
      stop("The clock table needs a `scaling` column (Scaled / YuGene).", call. = FALSE)
    }
    rep <- c(Scaled = "scaled_diff", YuGene = "yugene_diff")[as.character(model_paths$scaling)]
    if (anyNA(rep)) stop("Unknown `scaling` in the clock table; expected Scaled or YuGene.", call. = FALSE)
    path <- as.character(model_paths$path)
    type <- if (!is.null(mode)) rep(mode, length(path))
            else if ("type" %in% names(model_paths)) as.character(model_paths$type)
            else vapply(path, .tage_clock_type, character(1))
    jobs <- data.frame(eset = unname(rep), path = path, mode = type, column = basename(path),
                       stringsAsFactors = FALSE)
  } else {
    if (!is.list(model_paths) || length(model_paths) == 0 || is.null(names(model_paths))) {
      stop("model_paths must be a clock table with a `path` column or a named list of model paths.",
           call. = FALSE)
    }
    model_paths <- model_paths[names(model_paths) %in% valid_names]
    if (length(model_paths) == 0) {
      stop("No valid model paths found in model_paths. Valid names are: 'scaled', 'scaled_diff', 'yugene', 'yugene_diff'.")
    }
    by_file <- any(lengths(model_paths) > 1L)
    jobs <- do.call(rbind, lapply(names(model_paths), function(nm) {
      path <- as.character(model_paths[[nm]])
      type <- if (!is.null(mode)) rep(mode, length(path)) else vapply(path, .tage_clock_type, character(1))
      data.frame(eset = nm, path = path, mode = unname(type),
                 column = if (by_file) basename(path) else paste0(nm, "_", type, "_tAge"),
                 stringsAsFactors = FALSE)
    }))
  }

  missing_eset <- setdiff(jobs$eset, eset_names)
  if (length(missing_eset)) {
    stop("tAge_eset has no ", paste(unique(missing_eset), collapse = ", "),
         " element; pass the list returned by tAge_preprocessing().", call. = FALSE)
  }
  dup <- unique(jobs$column[duplicated(jobs$column)])
  if (length(dup)) stop("Several predictions would share the column ", paste(dup, collapse = ", "), ".",
                        call. = FALSE)
  jobs
}

# Model type from the registries or, for models outside them, the file name.
.tage_clock_type <- function(model_path) {
  fn <- basename(as.character(model_path))
  for (reg in list(tryCatch(.clock_registry(), error = function(e) NULL),
                   tryCatch(.module_clock_registry(), error = function(e) NULL))) {
    hit <- reg[reg$filename == fn, , drop = FALSE]
    if (nrow(hit) == 1) return(as.character(hit$type))
  }
  prefix <- sub("_.*$", "", fn)
  if (prefix %in% c("EN", "BR")) return(prefix)
  stop("Cannot tell whether ", fn, " is an elastic net or a Bayesian ridge model; pass mode.",
       call. = FALSE)
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
#' @param model_paths Clock table with a \code{path} column, or named list of
#'   model paths; see \code{\link{predict_tAge}}.
#' @inheritParams tAge_preprocessing
#' @inheritParams predict_tAge
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
  mode = NULL,
  control_group_column = NULL,
  control_group_label = NULL,
  count_threshold = 10,
  percent_threshold = 20,
  min_samples = 5,
  verbose = TRUE,
  gene_mapping_type = "auto",
  return_std = NULL,
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
