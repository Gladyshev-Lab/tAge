# Registry of published transcriptomic clock models and helpers to list and
# download them from Zenodo. The registry lists only the models that are
# publicly available in the Zenodo record referenced below.

# Zenodo record holding the published clock models.
.TAGE_ZENODO_RECORD <- "18763485"

# Load the bundled clock registry.
.clock_registry <- function() {
  path <- system.file("extdata", "clocks_metadata.csv", package = "tAge")
  if (path == "") {
    stop("clocks_metadata.csv not found in the installed package.")
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df$lifespan_scaled <- as.logical(df$lifespan_scaled)
  df
}

# Should a clock's output be multiplied by species maximum lifespan?
# Looks the model up in the registry by file name; falls back to a name-based
# heuristic, and returns NA if the type cannot be determined.
.clock_lifespan_scaled <- function(model_path) {
  fn  <- basename(as.character(model_path))
  reg <- tryCatch(.clock_registry(), error = function(e) NULL)
  if (!is.null(reg)) {
    hit <- reg[reg$filename == fn, , drop = FALSE]
    if (nrow(hit) == 1) return(isTRUE(hit$lifespan_scaled))
  }
  if (grepl("chronoage", fn, ignore.case = TRUE)) return(TRUE)
  if (grepl("mortality|normalizedage", fn, ignore.case = TRUE)) return(FALSE)
  NA
}

#' List available transcriptomic clock models
#'
#' Returns the registry of published clock models available on Zenodo, with
#' optional filtering. Pass a returned (filtered) data frame to
#' \code{\link{download_clocks}} to download the corresponding files.
#'
#' @param type Character. Filter by model type: "EN" (Elastic Net) or "BR"
#'   (Bayesian Ridge). Default NULL (no filter).
#' @param outcome Character. Filter by prediction outcome: "Chronological",
#'   "Mortality", or "Normalized age". Default NULL.
#' @param species Character. Filter by species group: "Mouse", "Rodents", or
#'   "Multispecies". Default NULL.
#' @param tissue Character. Filter by tissue (e.g. "Multi-Tissue", "Liver").
#'   Default NULL.
#' @param scaling Character. Filter by normalisation: "Scaled" or "YuGene".
#'   Default NULL.
#' @return A data frame with columns \code{filename}, \code{type},
#'   \code{outcome}, \code{species}, \code{tissue}, \code{scaling} and
#'   \code{lifespan_scaled} (whether the output is rescaled to age units).
#' @export
#' @examples
#' # All Elastic Net mortality clocks
#' list_clocks(type = "EN", outcome = "Mortality")
list_clocks <- function(type = NULL, outcome = NULL, species = NULL,
                        tissue = NULL, scaling = NULL) {
  df <- .clock_registry()
  keep <- function(df, col, val) if (is.null(val)) df else df[df[[col]] %in% val, , drop = FALSE]
  df <- keep(df, "type", type)
  df <- keep(df, "outcome", outcome)
  df <- keep(df, "species", species)
  df <- keep(df, "tissue", tissue)
  df <- keep(df, "scaling", scaling)
  rownames(df) <- NULL
  df
}

#' Download clock models from Zenodo
#'
#' Downloads the given clock model files from the Zenodo record and returns the
#' input augmented with a \code{path} column pointing to the local files.
#'
#' @param clocks Either a data frame returned by \code{\link{list_clocks}} (the
#'   \code{filename} column is used) or a character vector of model file names.
#' @param dest_dir Directory to save the models into. Created if needed.
#'   Default "clocks".
#' @param record Character Zenodo record id. Defaults to the published record.
#' @param overwrite Logical. Re-download files that already exist. Default FALSE.
#' @param quiet Logical. Suppress progress messages. Default FALSE.
#' @param timeout Seconds allowed per file. R's default of 60 s aborts the
#'   Bayesian ridge models, which are 0.9-2.4 GB each; elastic net models are
#'   about 1 MB. Default 3600.
#' @return The \code{clocks} data frame with an added \code{path} column.
#' @details A partial or failed download is removed rather than left on disk,
#'   and every file is checked to be a pickle (Zenodo answers some errors with
#'   an HTML page, which would otherwise be saved under the model's name).
#' @export
#' @examples
#' \dontrun{
#' clocks <- list_clocks(type = "EN", outcome = "Mortality",
#'                       species = "Multispecies", tissue = "Multi-Tissue")
#' clocks <- download_clocks(clocks, dest_dir = "clocks")
#' model_paths <- list(
#'   scaled_diff = clocks$path[clocks$scaling == "Scaled"],
#'   yugene_diff = clocks$path[clocks$scaling == "YuGene"]
#' )
#' }
download_clocks <- function(clocks, dest_dir = "clocks",
                            record = .TAGE_ZENODO_RECORD,
                            overwrite = FALSE, quiet = FALSE,
                            timeout = 3600) {
  if (is.data.frame(clocks)) {
    if (!"filename" %in% colnames(clocks)) {
      stop("`clocks` data frame must have a 'filename' column (use list_clocks()).")
    }
    files <- as.character(clocks$filename)
  } else {
    files <- as.character(clocks)
    clocks <- data.frame(filename = files, stringsAsFactors = FALSE)
  }
  if (length(files) == 0) stop("No clocks to download.")

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- character(length(files))

  old_timeout <- getOption("timeout", 60)
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout, timeout))

  for (i in seq_along(files)) {
    fn   <- files[i]
    dest <- file.path(dest_dir, fn)
    if (file.exists(dest) && !overwrite) {
      if (!quiet) message(sprintf("\u2713 Already present: %s", fn))
    } else {
      url <- sprintf("https://zenodo.org/records/%s/files/%s?download=1",
                     record, utils::URLencode(fn, reserved = TRUE))
      if (!quiet) message(sprintf("Downloading %s ...", fn))
      .tage_download_file(url, dest, quiet = quiet)
    }
    paths[i] <- dest
  }

  clocks$path <- paths
  clocks
}

# download.file() leaves a truncated file behind when it fails or times out,
# and a later call would then report it as "already present".
.tage_download_file <- function(url, dest, quiet = FALSE) {
  status <- tryCatch(
    utils::download.file(url, dest, mode = "wb", quiet = quiet),
    error = function(e) e, warning = function(w) w
  )
  if (inherits(status, "condition") || !identical(as.integer(status), 0L)) {
    unlink(dest)
    msg <- if (inherits(status, "condition")) conditionMessage(status) else "non-zero exit status"
    stop(sprintf("Download of %s failed (%s). The partial file was removed.", basename(dest), msg),
         call. = FALSE)
  }
  .tage_check_download(dest)
  invisible(dest)
}

# joblib pickles start with the pickle PROTO opcode (0x80); an HTML error page
# starts with "<". Anything else that is empty is an aborted transfer.
.tage_check_download <- function(dest) {
  size <- file.info(dest)$size
  first <- if (isTRUE(size > 0)) readBin(dest, "raw", n = 1L) else raw(0)
  ok <- length(first) == 1L && first == as.raw(0x80)
  if (!ok) {
    head <- if (length(first)) rawToChar(readBin(dest, "raw", n = min(size, 200L))) else ""
    unlink(dest)
    stop(sprintf("%s is not a model file (%s bytes%s); it was removed. Check the record id and the file name.",
                 basename(dest), format(size),
                 if (grepl("<", head, fixed = TRUE)) ", looks like an HTML page" else ""),
         call. = FALSE)
  }
  invisible(TRUE)
}
