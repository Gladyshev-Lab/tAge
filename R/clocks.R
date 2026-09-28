# Registry of the published transcriptomic clock models and helpers to list
# and download them from Zenodo. The registry mirrors the Zenodo record
# referenced below: the composite clocks as single files and the module
# clocks as two archives.

# Zenodo record holding the published clock models.
.TAGE_ZENODO_RECORD <- "22166800"

# Load a bundled registry table.
.tage_registry <- function(file) {
  path <- system.file("extdata", file, package = "tAge")
  if (path == "") stop(file, " not found in the installed package.")
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.clock_registry <- function() {
  df <- .tage_registry("clocks_metadata.csv")
  df$lifespan_scaled <- as.logical(df$lifespan_scaled)
  df
}

.module_clock_registry <- function() .tage_registry("module_clocks_metadata.csv")

.tage_filter <- function(df, ...) {
  filters <- list(...)
  for (col in names(filters)) {
    if (!is.null(filters[[col]])) df <- df[df[[col]] %in% filters[[col]], , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

# What a clock predicts, from the registries or -- for models outside them --
# from the file name: "Chronological", "Mortality", "Normalized age",
# "Lifespan" or NA. Module clocks predict the same quantities as the composite
# clocks and are scaled the same way.
.clock_outcome <- function(model_path) {
  fn  <- basename(as.character(model_path))
  for (reg in list(tryCatch(.clock_registry(), error = function(e) NULL),
                   tryCatch(.module_clock_registry(), error = function(e) NULL))) {
    hit <- reg[reg$filename == fn, , drop = FALSE]
    if (nrow(hit) == 1) return(as.character(hit$outcome))
  }
  key <- tolower(fn)
  if (grepl("chronoage|chronological", key)) return("Chronological")
  if (grepl("hazard|mortality", key))        return("Mortality")
  if (grepl("relage|normalizedage", key))    return("Normalized age")
  if (grepl("lifespan", key))                return("Lifespan")
  NA_character_
}

# Should a clock's output be multiplied by species maximum lifespan?
.clock_lifespan_scaled <- function(model_path) {
  outcome <- .clock_outcome(model_path)
  if (is.na(outcome)) return(NA)
  identical(outcome, "Chronological")
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
  .tage_filter(.clock_registry(), type = type, outcome = outcome, species = species,
               tissue = tissue, scaling = scaling)
}

#' List available module clock models
#'
#' Returns the registry of published module clocks: one elastic net clock per
#' co-expression module (plus \code{"allmodulegenes"}, trained on the genes of
#' all modules), for the rodent and the multispecies module sets. Pass a
#' returned (filtered) data frame to \code{\link{download_clocks}}, which
#' fetches the archive each set is published as.
#'
#' @param outcome Character. "Chronological" or "Mortality". Default NULL.
#' @param species Character. Module set: "Rodents" or "Multispecies". Default
#'   NULL.
#' @param color Character. Module colour(s), e.g. "blue". Default NULL.
#' @return A data frame with columns \code{filename}, \code{type},
#'   \code{outcome}, \code{species}, \code{tissue}, \code{scaling},
#'   \code{color}, \code{function} (the module's annotated biological process)
#'   and \code{archive} (the Zenodo archive holding the file).
#' @export
#' @examples
#' list_module_clocks(outcome = "Mortality", species = "Rodents")
list_module_clocks <- function(outcome = NULL, species = NULL, color = NULL) {
  .tage_filter(.module_clock_registry(), outcome = outcome, species = species, color = color)
}

#' Download clock models from Zenodo
#'
#' Downloads the given clock models from the Zenodo record and returns the
#' input augmented with a \code{path} column pointing to the local files.
#' Composite clocks (\code{\link{list_clocks}}) are single files; module
#' clocks (\code{\link{list_module_clocks}}) come in one archive per module
#' set, which is downloaded once and unpacked into \code{dest_dir}.
#'
#' @param clocks A data frame returned by \code{\link{list_clocks}} or
#'   \code{\link{list_module_clocks}}, or a character vector of composite clock
#'   file names.
#' @param dest_dir Directory to save the models into. Created if needed.
#'   Default "clocks".
#' @param record Character Zenodo record id. Defaults to the published record.
#' @param overwrite Logical. Re-download files that already exist. Default FALSE.
#' @param quiet Logical. Suppress progress messages. Default FALSE.
#' @param timeout Seconds allowed per file. R's default of 60 s aborts the
#'   Bayesian ridge models, which are 0.9-2.4 GB each; elastic net models are
#'   about 1 MB and the module clock archives under 1 MB. Default 3600.
#' @return The \code{clocks} data frame with an added \code{path} column.
#' @details The whole record is about 64 GB, nearly all of it the 30 Bayesian
#'   ridge models. A partial or failed download is removed rather than left on
#'   disk, and every file is checked to be a pickle or a zip archive (Zenodo
#'   answers some errors with an HTML page, which would otherwise be saved
#'   under the model's name).
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
#'
#' modules <- download_clocks(list_module_clocks(outcome = "Mortality", species = "Rodents"),
#'                            dest_dir = "clocks")
#' }
download_clocks <- function(clocks, dest_dir = "clocks",
                            record = .TAGE_ZENODO_RECORD,
                            overwrite = FALSE, quiet = FALSE,
                            timeout = 3600) {
  if (is.data.frame(clocks)) {
    if (!"filename" %in% colnames(clocks)) {
      stop("`clocks` data frame must have a 'filename' column (use list_clocks() or list_module_clocks()).")
    }
    files <- as.character(clocks$filename)
  } else {
    files <- as.character(clocks)
    clocks <- data.frame(filename = files, stringsAsFactors = FALSE)
  }
  if (length(files) == 0) stop("No clocks to download.")
  archives <- if ("archive" %in% colnames(clocks)) as.character(clocks$archive) else rep(NA_character_, length(files))

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- character(length(files))

  old_timeout <- getOption("timeout", 60)
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout, timeout))

  fetch <- function(name, dest) {
    if (file.exists(dest) && !overwrite) {
      if (!quiet) message(sprintf("\u2713 Already present: %s", name))
      return(invisible(FALSE))
    }
    url <- sprintf("https://zenodo.org/records/%s/files/%s?download=1",
                   record, utils::URLencode(name, reserved = TRUE))
    if (!quiet) message(sprintf("Downloading %s ...", name))
    .tage_download_file(url, dest, quiet = quiet)
    invisible(TRUE)
  }

  # Module clocks: one archive per set, unpacked into a folder of the same name.
  for (archive in unique(archives[!is.na(archives)])) {
    zip <- file.path(dest_dir, paste0(archive, ".zip"))
    fresh <- fetch(paste0(archive, ".zip"), zip)
    if (fresh || !dir.exists(file.path(dest_dir, archive))) {
      utils::unzip(zip, exdir = dest_dir, overwrite = TRUE)
    }
  }

  for (i in seq_along(files)) {
    if (is.na(archives[i])) {
      dest <- file.path(dest_dir, files[i])
      fetch(files[i], dest)
    } else {
      dest <- file.path(dest_dir, archives[i], files[i])
      if (!file.exists(dest)) {
        stop(sprintf("%s is not in the unpacked archive %s.", files[i], archives[i]), call. = FALSE)
      }
    }
    paths[i] <- dest
  }

  clocks$path <- paths
  clocks
}

# The transfer goes to `<dest>.part` and is renamed only once it is complete
# and checked: a failed or interrupted download never leaves a file under the
# model's name.
.tage_download_file <- function(url, dest, quiet = FALSE) {
  part <- paste0(dest, ".part")
  unlink(part)
  status <- tryCatch(
    utils::download.file(url, part, mode = "wb", quiet = quiet),
    error = function(e) e, warning = function(w) w
  )
  if (inherits(status, "condition") || !identical(as.integer(status), 0L)) {
    unlink(part)
    msg <- if (inherits(status, "condition")) conditionMessage(status) else "non-zero exit status"
    stop(sprintf("Download of %s failed (%s). The partial file was removed.", basename(dest), msg),
         call. = FALSE)
  }
  .tage_check_download(part, shown_as = basename(dest))
  if (!file.rename(part, dest)) {
    unlink(part)
    stop(sprintf("Could not move the downloaded file into place: %s", dest), call. = FALSE)
  }
  invisible(dest)
}

# joblib pickles start with the pickle PROTO opcode (0x80) and zip archives
# with "PK"; an HTML error page starts with "<". Anything else that is empty
# is an aborted transfer.
.tage_check_download <- function(dest, shown_as = basename(dest)) {
  size <- file.info(dest)$size
  first <- if (isTRUE(size > 0)) readBin(dest, "raw", n = 2L) else raw(0)
  ok <- length(first) >= 1L &&
    (first[1] == as.raw(0x80) || identical(first, charToRaw("PK")))
  if (!ok) {
    head <- if (length(first)) rawToChar(readBin(dest, "raw", n = min(size, 200L))) else ""
    unlink(dest)
    stop(sprintf("%s is not a model file (%s bytes%s); it was removed. Check the record id and the file name.",
                 shown_as, format(size),
                 if (grepl("<", head, fixed = TRUE)) ", looks like an HTML page" else ""),
         call. = FALSE)
  }
  invisible(TRUE)
}
