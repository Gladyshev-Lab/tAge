# The tAge figure style: one palette, one type scale, one anatomy -- the same
# tokens as the Python package, so figures from either language read as one
# system. Colour is assigned by the job it does: categorical slots for groups
# (the reference group in neutral grey), one colour per clock outcome, a
# diverging scale centred on zero for signed effects, a single hue for
# magnitudes. Text never wears a data colour.

TAGE_SURFACE        <- "#ffffff"
TAGE_TEXT_PRIMARY   <- "#0b0b0b"
TAGE_TEXT_SECONDARY <- "#52514e"
TAGE_TEXT_MUTED     <- "#8a8987"
TAGE_GRID_COLOR     <- "#e6e5e1"
TAGE_AXIS_COLOR     <- "#8a8987"

# Colour-vision-safe categorical slots, in the order they are assigned.
TAGE_CATEGORICAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                      "#e87ba4", "#008300", "#4a3aa7", "#e34948")
# The reference / control group recedes; the comparisons carry the colour.
TAGE_REFERENCE_COLOR <- "#8a8987"

# Signed effects: blue (negative) - neutral grey - red (positive).
TAGE_DIVERGING <- c("#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f2aba9", "#e34948", "#8c1d1d")
# Magnitude: one hue, light -> dark.
TAGE_SEQUENTIAL <- c("#e4eefb", "#b7d3f6", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b")

#' The tAge colour tokens
#'
#' The colours every figure in the package is drawn with, for callers who
#' build their own plots and want them to match: categorical slots for groups,
#' the reference-group grey, one colour per clock outcome, the diverging
#' ramp for signed effects and the sequential ramp for magnitudes, plus the
#' text and grid tokens.
#'
#' @return A named list of hex colours and colour vectors.
#' @examples
#' tage_colors()$categorical
#' @export
tage_colors <- function() {
  list(
    categorical    = TAGE_CATEGORICAL,
    reference      = TAGE_REFERENCE_COLOR,
    outcome        = TAGE_OUTCOME_COLORS,
    diverging      = TAGE_DIVERGING,
    sequential     = TAGE_SEQUENTIAL,
    text_primary   = TAGE_TEXT_PRIMARY,
    text_secondary = TAGE_TEXT_SECONDARY,
    text_muted     = TAGE_TEXT_MUTED,
    grid           = TAGE_GRID_COLOR,
    axis           = TAGE_AXIS_COLOR,
    surface        = TAGE_SURFACE
  )
}

#' Colour per group level
#'
#' Assigns the categorical slots in the order the levels are given, with the
#' reference level in the neutral reference grey, so a colour stays with its
#' group whatever else is on the figure. This is how every group-coloured
#' tAge figure picks its colours.
#'
#' @param levels Group levels, in display order.
#' @param reference Level drawn in the reference grey. Default \code{NULL}
#'   (none).
#' @param palette Override: a named vector \code{c(level = colour)}, or an
#'   unnamed vector of colours in level order. Levels missing from a named
#'   vector are filled in automatically.
#'
#' @return A named character vector of hex colours.
#' @examples
#' tage_series_colors(c("WT", "KO", "HET"), reference = "WT")
#' @export
tage_series_colors <- function(levels, reference = NULL, palette = NULL) {
  levels <- as.character(levels)
  if (!is.null(palette)) {
    if (!is.null(names(palette)) && any(nzchar(names(palette)))) {
      out <- stats::setNames(unname(palette[levels]), levels)
      missing <- levels[is.na(out)]
      if (length(missing)) out[missing] <- tage_series_colors(missing, reference = reference)
      return(out)
    }
    if (length(palette) < length(levels)) {
      stop(sprintf("palette has %d colours for %d groups", length(palette), length(levels)),
           call. = FALSE)
    }
    return(stats::setNames(as.character(palette[seq_along(levels)]), levels))
  }
  out <- character(length(levels))
  slot <- 0L
  for (i in seq_along(levels)) {
    if (!is.null(reference) && identical(levels[i], as.character(reference))) {
      out[i] <- TAGE_REFERENCE_COLOR
    } else {
      out[i] <- TAGE_CATEGORICAL[slot %% length(TAGE_CATEGORICAL) + 1L]
      slot <- slot + 1L
    }
  }
  stats::setNames(out, levels)
}

# Move a colour towards white: the fill under a coloured edge.
.tage_tint <- function(colour, amount = 0.65) {
  rgb <- grDevices::col2rgb(colour) / 255
  grDevices::rgb(rgb[1] + (1 - rgb[1]) * amount, rgb[2] + (1 - rgb[2]) * amount,
                 rgb[3] + (1 - rgb[3]) * amount)
}

# Text colour that clears a filled cell: ink on light fills, white on dark.
.tage_ink_on <- function(colour) {
  rgb <- grDevices::col2rgb(colour) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum > 0.55, TAGE_TEXT_PRIMARY, TAGE_SURFACE)
}

#' The tAge ggplot2 theme
#'
#' Recessive axes (left and bottom only), a light grid on the value axis,
#' a left-aligned title with the subtitle in secondary ink and the caption in
#' muted ink, facet strips as plain semibold text, legend on the right. Every
#' tAge figure uses it; add it to your own plots to match.
#'
#' @param base_size Base font size in points.
#' @param grid \code{"y"} (default), \code{"x"}, \code{"both"} or
#'   \code{"none"}: which axis gets the light grid.
#'
#' @return A ggplot2 theme.
#' @examples
#' \dontrun{
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'   ggplot2::geom_point() +
#'   theme_tage()
#' }
#' @export
theme_tage <- function(base_size = 10, grid = c("y", "x", "both", "none")) {
  .tage_gg_require()
  grid <- match.arg(grid)
  el <- ggplot2::element_line(colour = TAGE_GRID_COLOR, linewidth = 0.3)
  none <- ggplot2::element_blank()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = TAGE_TEXT_PRIMARY),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1,
                                         colour = TAGE_TEXT_PRIMARY, hjust = 0),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = TAGE_TEXT_SECONDARY,
                                            hjust = 0, margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(size = base_size - 2.5, colour = TAGE_TEXT_MUTED,
                                           hjust = 0, margin = ggplot2::margin(t = 8)),
      plot.caption.position = "plot",
      plot.background = ggplot2::element_rect(fill = TAGE_SURFACE, colour = NA),
      panel.background = ggplot2::element_rect(fill = TAGE_SURFACE, colour = NA),
      panel.grid.minor = none,
      panel.grid.major.x = if (grid %in% c("x", "both")) el else none,
      panel.grid.major.y = if (grid %in% c("y", "both")) el else none,
      axis.line = ggplot2::element_line(colour = TAGE_AXIS_COLOR, linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = TAGE_AXIS_COLOR, linewidth = 0.4),
      axis.ticks.length = ggplot2::unit(2.5, "pt"),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = TAGE_TEXT_SECONDARY),
      axis.title = ggplot2::element_text(size = base_size, colour = TAGE_TEXT_SECONDARY),
      strip.background = none,
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1,
                                         colour = TAGE_TEXT_PRIMARY, hjust = 0),
      strip.placement = "outside",
      legend.position = "right",
      legend.title = ggplot2::element_text(size = base_size - 1, colour = TAGE_TEXT_SECONDARY),
      legend.text = ggplot2::element_text(size = base_size - 1, colour = TAGE_TEXT_PRIMARY),
      legend.key = none,
      legend.background = none
    )
}

#' Diverging fill scale for signed effects
#'
#' Blue for negative, neutral grey at zero, red for positive - the scale every
#' tAge heatmap uses. The limits are made symmetric around zero.
#'
#' @param limit Absolute limit of the scale; clamp values to it beforehand.
#' @param name Legend title.
#' @param ... Passed to \code{ggplot2::scale_fill_gradientn()}.
#'
#' @return A ggplot2 scale.
#' @export
scale_fill_tage_diverging <- function(limit, name = ggplot2::waiver(), ...) {
  .tage_gg_require()
  ggplot2::scale_fill_gradientn(colours = TAGE_DIVERGING, limits = c(-limit, limit),
                                name = name, ...)
}
