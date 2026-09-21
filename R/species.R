# Species maximum lifespans used to turn the clocks' normalised age back into
# age units, and the units each species is reported in. The values are those
# of the TACO reference application (AnAge maximum longevity); the paper text
# quotes 3.8 y for rat and 122 y for human, so check with the authors before
# relying on rat ages to better than ~10%.
.TAGE_SPECIES <- data.frame(
  species            = c("mouse", "rat", "human", "monkey"),
  max_lifespan_years = c(4, 4.2, 122.5, 39),
  default_units      = c("months", "months", "years", "years"),
  stringsAsFactors   = FALSE
)

#' Species supported by the clocks
#'
#' The species accepted by \code{\link{tAge_preprocessing}} and
#' \code{\link{predict_tAge}}, with the maximum lifespan used to rescale
#' chronological-age clocks and the age units each species is reported in by
#' default (months for rodents, years for primates).
#'
#' @return A data frame with columns \code{species}, \code{max_lifespan_years}
#'   and \code{default_units}.
#' @examples
#' tage_species()
#' @export
tage_species <- function() {
  .TAGE_SPECIES
}

.tage_species_row <- function(species) {
  if (is.null(species) || length(species) != 1L || is.na(species)) {
    stop("`species` must be one of: ", paste(.TAGE_SPECIES$species, collapse = ", "),
         call. = FALSE)
  }
  key <- tolower(as.character(species))
  hit <- .TAGE_SPECIES[.TAGE_SPECIES$species == key, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop("Unknown species '", species, "'. Supported: ",
         paste(.TAGE_SPECIES$species, collapse = ", "), " (see tage_species()).",
         call. = FALSE)
  }
  hit
}

# Multiplier turning a chronological clock's output (fraction of the species
# maximum lifespan) into `age_units`.
.tage_lifespan_factor <- function(species, age_units = c("auto", "months", "years")) {
  age_units <- match.arg(age_units)
  row <- .tage_species_row(species)
  if (age_units == "auto") age_units <- row$default_units
  factor <- row$max_lifespan_years
  if (age_units == "months") factor <- factor * 12
  list(factor = factor, units = age_units)
}

# tAge_preprocessing() records the species in the ExpressionSet so that
# predict_tAge() does not need to be told again.
.tage_set_species <- function(eset, species) {
  ed <- Biobase::experimentData(eset)
  ed@other$tage_species <- tolower(as.character(species))
  Biobase::experimentData(eset) <- ed
  eset
}

.tage_get_species <- function(eset) {
  other <- tryCatch(Biobase::notes(eset), error = function(e) NULL)
  sp <- if (is.list(other)) other$tage_species else NULL
  if (is.null(sp) || !nzchar(sp)) NULL else sp
}

.tage_resolve_species <- function(species, eset) {
  if (is.null(species)) species <- .tage_get_species(eset)
  if (is.null(species)) {
    stop("`species` is unknown: pass species = \"mouse\" / \"rat\" / \"human\" / \"monkey\", ",
         "or preprocess with tAge_preprocessing(), which records it in the ExpressionSet.",
         call. = FALSE)
  }
  .tage_species_row(species)$species
}
