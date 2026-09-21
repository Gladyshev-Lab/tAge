# Publication-style figures built on the tAge statistics.
#
# Both functions take per-sample predictions, run the corresponding test from
# stats.R, and draw the result, so the figure and the numbers behind it can
# never drift apart. The layout mirrors the Python package's
# tage.pl.plot_clock_forest() and tage.pl.plot_module_heatmap(), so the two
# languages produce comparable figures from the same data.

# One colour per outcome, kept identical to the Python package.
TAGE_OUTCOME_COLORS <- c(
  "Chronological" = "#2a78d6",
  "Mortality"     = "#e34948",
  "Lifespan"      = "#4a3aa7",
  "NormalizedAge" = "#4a3aa7"
)
TAGE_OUTCOME_ORDER <- c("Chronological", "NormalizedAge", "Lifespan", "Mortality")

TAGE_OUTCOME_UNITS <- c(
  "Chronological" = "months",
  "Mortality"     = "log10 hazard ratio",
  "Lifespan"      = "fraction of max lifespan",
  "NormalizedAge" = "fraction of max lifespan"
)

# The registry writes "Normalized age"; the figures key panels by the
# one-word form. Anything unknown is passed through unchanged.
.tage_canonical_outcome <- function(x) {
  x <- as.character(x)
  key <- gsub("[^a-z]", "", tolower(x))
  known <- stats::setNames(TAGE_OUTCOME_ORDER, gsub("[^a-z]", "", tolower(TAGE_OUTCOME_ORDER)))
  hit <- unname(known[key])
  ifelse(is.na(hit), x, hit)
}

# Map a clock table (list_clocks() output, or any table with a `filename`
# column) onto the prediction columns of `data`. predict_tAge() names its
# output <normalisation>_<mode>_tAge, never by the model file, so a registry
# row whose `filename` is not a column of `data` is matched through its
# `scaling` and `type` instead: Scaled + EN -> scaled_diff_EN_tAge.
.tage_resolve_clock_columns <- function(clocks_meta, data, value_columns,
                                        label_column, outcome_column) {
  clocks_meta <- as.data.frame(clocks_meta, stringsAsFactors = FALSE)
  if (!"filename" %in% names(clocks_meta)) {
    stop("`clocks_meta` needs a 'filename' column naming the prediction columns.",
         call. = FALSE)
  }
  cols   <- as.character(clocks_meta$filename)
  direct <- cols %in% colnames(data)

  if (!all(direct) && all(c("scaling", "type") %in% names(clocks_meta))) {
    norm <- c(Scaled = "scaled_diff", YuGene = "yugene_diff")[as.character(clocks_meta$scaling)]
    derived <- paste0(norm, "_", as.character(clocks_meta$type), "_tAge")
    use <- !direct & !is.na(norm) & derived %in% colnames(data)
    cols[use] <- derived[use]
  }

  keep <- cols %in% colnames(data)
  if (!is.null(value_columns)) keep <- keep & cols %in% value_columns
  clocks_meta <- clocks_meta[keep, , drop = FALSE]
  cols <- cols[keep]

  if (anyDuplicated(cols)) {
    stop("Several rows of `clocks_meta` map onto the same prediction column (",
         paste(unique(cols[duplicated(cols)]), collapse = ", "), "). predict_tAge() ",
         "names its columns by normalisation and model type only, so keep one clock ",
         "outcome per results table, or rename the columns and put the new names in ",
         "`filename`.", call. = FALSE)
  }

  labels <- if (label_column %in% names(clocks_meta)) {
    as.character(clocks_meta[[label_column]])
  } else {
    .tage_clock_labels(clocks_meta, cols)
  }
  outcomes <- if (outcome_column %in% names(clocks_meta)) {
    .tage_canonical_outcome(clocks_meta[[outcome_column]])
  } else {
    rep(NA_character_, length(cols))
  }
  list(cols = cols, labels = labels, outcomes = outcomes)
}

# Row label built from the registry fields when no `name` column exists.
.tage_clock_labels <- function(meta, fallback) {
  parts <- intersect(c("type", "species", "tissue", "scaling"), names(meta))
  if (length(parts) == 0L || nrow(meta) == 0L) return(fallback)
  lab <- do.call(paste, c(lapply(parts, function(p) as.character(meta[[p]])), sep = " \u00b7 "))
  ifelse(is.na(lab) | !nzchar(lab), fallback, lab)
}

.tage_gg_require <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for tAge figures.", call. = FALSE)
  }
}

.tage_outcome_unit <- function(outcome, units = NULL) {
  merged <- c(TAGE_OUTCOME_UNITS, units)
  out <- unname(merged[as.character(outcome)])
  out[is.na(out)] <- "effect"
  out
}

#' Module colour to biological function map
#'
#' Reads the bundled annotation that names each co-expression module, e.g.
#' \code{"blue"} -> \code{"Respiration/Mitochondrial translation"}.
#'
#' @param version Module set: \code{"4.6"} (default), \code{"5.4"} or
#'   \code{"human"}.
#'
#' @return A named character vector, module colour to function, or an empty
#'   vector when the annotation is unavailable.
#'
#' @examples
#' head(load_module_functions("4.6"))
#' @export
load_module_functions <- function(version = c("4.6", "5.4", "human")) {
  version <- match.arg(version)
  file <- switch(version,
                 "4.6"   = "Module_to_function_map.csv",
                 "5.4"   = "Module_to_function_map_54.csv",
                 "human" = "Module_to_function_annotation_human.csv")
  path <- system.file("extdata", "modules", file, package = "tAge")
  if (!nzchar(path) || !file.exists(path)) return(stats::setNames(character(0), character(0)))
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("Module", "Function") %in% names(tab))) {
    return(stats::setNames(character(0), character(0)))
  }
  stats::setNames(as.character(tab$Function), as.character(tab$Module))
}

# Attach the size a figure wants to be printed at, and -- inside a notebook --
# ask the renderer for it. A ggplot carries no size of its own, so without this
# the caller has to rediscover a sensible one for every clock count.
.tage_set_size <- function(p, width, height, apply_option = TRUE) {
  size <- c(width = as.numeric(width), height = as.numeric(height))
  attr(p, "tage_size") <- size
  if (isTRUE(apply_option)) {
    options(repr.plot.width = size[["width"]], repr.plot.height = size[["height"]])
  }
  p
}

#' Size a tAge figure was designed for
#'
#' The figure functions grow the canvas with the number of clocks, modules and
#' panels they were given, and record the result. \code{tage_fig_size()} reads
#' it back so the same size can be reused for \code{ggsave()}, a knitr chunk
#' (\code{fig.width} / \code{fig.height}) or a manual tweak.
#'
#' @param p A plot from \code{\link{tage_clock_forest}} or
#'   \code{\link{tage_module_heatmap}}.
#'
#' @return Named numeric vector with \code{width} and \code{height} in inches,
#'   or \code{NULL} when the plot carries no recorded size.
#'
#' @examples
#' \dontrun{
#' p <- tage_clock_forest(results, clocks_meta = clocks,
#'                        group_column = "Genotype", reference_group = "WT")
#' tage_fig_size(p)
#' }
#' @export
tage_fig_size <- function(p) attr(p, "tage_size")

#' Save a tAge figure at the size it was designed for
#'
#' Thin wrapper over \code{ggplot2::ggsave()} that defaults \code{width} and
#' \code{height} to the size recorded by the figure function, so a figure with
#' thirty modules is not squeezed into the same canvas as one with three.
#'
#' @param p A plot from \code{\link{tage_clock_forest}} or
#'   \code{\link{tage_module_heatmap}}.
#' @param filename Output path; the extension picks the device.
#' @param width,height Size in inches. \code{NULL} (default) uses the recorded
#'   size, falling back to ggplot2's own default when there is none.
#' @param dpi Resolution for raster devices.
#' @param ... Passed to \code{ggplot2::ggsave()}.
#'
#' @return The path, invisibly.
#'
#' @examples
#' \dontrun{
#' p <- tage_clock_forest(results, clocks_meta = clocks,
#'                        group_column = "Genotype", reference_group = "WT")
#' tage_save_plot(p, "forest.png")                 # recorded size
#' tage_save_plot(p, "forest.pdf", width = 7)      # fixed width, recorded height
#' }
#' @export
tage_save_plot <- function(p, filename, width = NULL, height = NULL,
                           dpi = 300, ...) {
  .tage_gg_require()
  size <- attr(p, "tage_size")
  if (is.null(width) && !is.null(size))  width  <- unname(size[["width"]])
  if (is.null(height) && !is.null(size)) height <- unname(size[["height"]])
  args <- list(filename = filename, plot = p, dpi = dpi, ...)
  if (!is.null(width))  args$width  <- width
  if (!is.null(height)) args$height <- height
  do.call(ggplot2::ggsave, args)
  invisible(filename)
}

# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------

#' Forest plot of clock effects
#'
#' One row per clock: the effect between two groups, whiskers spanning the
#' confidence interval, filled when the adjusted p-value clears
#' \code{sig_threshold} and hollow otherwise. Panels are laid out as outcome
#' (rows) by stratum (columns). This is the figure to reach for when several
#' clocks are applied to one comparison -- a heatmap hides the uncertainty and a
#' box plot per clock does not fit on a page.
#'
#' The statistics come from \code{\link{tage_compare_groups}}, so the figure and
#' the numbers behind it cannot drift apart.
#'
#' @param data Data frame of per-sample predictions, e.g. from
#'   \code{\link{predict_tAge}}.
#' @param value_columns Prediction columns to show. May be omitted when
#'   \code{clocks_meta} is given, in which case its \code{filename} column is used.
#' @param group_column Column holding the experimental groups.
#' @param reference_group Level everything is compared against.
#' @param compare_groups Levels compared with the reference. \code{NULL} uses
#'   every other level; more than one is drawn as separate coloured series.
#' @param clocks_meta Clock table from \code{\link{list_clocks}}. Supplies the
#'   row labels and the outcome each clock predicts, which sets the panel rows,
#'   the colours and the x-axis units.
#' @param covariates,split_by,se_columns,variance_strata,p_adjust,p_adjust_scope,conf_level
#'   Passed to \code{\link{tage_compare_groups}}. \code{p_adjust_scope} defaults
#'   to \code{"across_columns"}, i.e. correction across the clocks in one panel,
#'   which is what "filled = significant" implies.
#' @param sig_threshold Adjusted p-value below which a marker is filled.
#' @param sort_by_effect Order clocks by effect size within each panel.
#' @param label_column,outcome_column Columns of \code{clocks_meta} holding the
#'   display label and the outcome.
#' @param units Named character vector overriding the x-axis unit per outcome,
#'   e.g. \code{c(Chronological = "years")} for human data.
#' @param title,subtitle,caption Figure text. \code{subtitle} defaults to a
#'   description of the markers; pass \code{NULL} to drop it.
#' @param base_size Base font size.
#' @param row_height,panel_width Inches per clock row and per panel column, used
#'   to size the canvas so that rows stay legible as the clock set grows.
#' @param width,height Absolute figure size in inches, overriding
#'   \code{row_height} and \code{panel_width}.
#' @param stats Pre-computed \code{\link{tage_compare_groups}} output to plot
#'   instead of recomputing.
#'
#' @return A \code{ggplot} object. The statistics are attached as the
#'   \code{"tage_stats"} attribute and the intended size as
#'   \code{"tage_size"}; see \code{\link{tage_save_plot}}.
#'
#' @seealso \code{\link{tage_compare_groups}}, \code{\link{tage_module_heatmap}}
#'
#' @examples
#' \dontrun{
#' clocks <- list_clocks(type = "EN", tissue = "Multi-Tissue")
#' tage_clock_forest(results, clocks_meta = clocks,
#'                   group_column = "Genotype", reference_group = "WT",
#'                   split_by = "Tissue")
#' }
#' @export
tage_clock_forest <- function(data,
                              value_columns = NULL,
                              group_column,
                              reference_group = NULL,
                              compare_groups = NULL,
                              clocks_meta = NULL,
                              covariates = NULL,
                              split_by = NULL,
                              se_columns = NULL,
                              variance_strata = c("subset", "all_data"),
                              p_adjust = "BH",
                              p_adjust_scope = "across_columns",
                              conf_level = 0.95,
                              sig_threshold = 0.05,
                              sort_by_effect = TRUE,
                              label_column = "name",
                              outcome_column = "outcome",
                              units = NULL,
                              title = "Clock effects",
                              subtitle = NA,
                              caption = NULL,
                              base_size = 11,
                              row_height = 0.30,
                              panel_width = 4.9,
                              width = NULL,
                              height = NULL,
                              stats = NULL) {
  .tage_gg_require()
  variance_strata <- match.arg(variance_strata)

  # ---- resolve which columns to show, and how to label them ----
  if (!is.null(clocks_meta)) {
    resolved <- .tage_resolve_clock_columns(clocks_meta, data, value_columns,
                                            label_column, outcome_column)
    cols     <- resolved$cols
    labels   <- resolved$labels
    outcomes <- resolved$outcomes
  } else {
    if (is.null(value_columns)) stop("Pass `value_columns`, `clocks_meta`, or both.", call. = FALSE)
    cols <- intersect(as.character(value_columns), colnames(data))
    labels <- cols
    outcomes <- rep(NA_character_, length(cols))
  }
  if (length(cols) == 0L) {
    stop("None of the requested clock columns are present in `data`.", call. = FALSE)
  }
  # Named clock_label, not label: tage_compare_groups() already returns a
  # `label` column holding the significance stars, and merging both would
  # silently produce label.x / label.y.
  lookup <- data.frame(value_column = cols, clock_label = labels, outcome = outcomes,
                       stringsAsFactors = FALSE)

  # ---- statistics ----
  if (is.null(stats)) {
    if (is.null(reference_group)) {
      stop("`reference_group` is required unless `stats` is supplied.", call. = FALSE)
    }
    stats <- tage_compare_groups(
      data = data, value_columns = cols, group_column = group_column,
      reference_group = reference_group, compare_groups = compare_groups,
      covariates = covariates, split_by = split_by, se_columns = se_columns,
      method = "trt.vs.ctrl", variance_strata = variance_strata,
      p_adjust = p_adjust, p_adjust_scope = p_adjust_scope, conf_level = conf_level
    )
  }
  if (nrow(stats) == 0L) {
    stop("No contrasts to plot - check the group labels and sample sizes.", call. = FALSE)
  }

  d <- merge(stats, lookup, by = "value_column", all.x = TRUE, sort = FALSE)
  d$clock_label[is.na(d$clock_label)] <- d$value_column[is.na(d$clock_label)]
  d$significant <- !is.na(d$p_adjusted) & d$p_adjusted < sig_threshold

  # Outcome drives the panel rows, the colour and the units.
  known <- intersect(TAGE_OUTCOME_ORDER, unique(stats::na.omit(d$outcome)))
  extra <- setdiff(unique(stats::na.omit(d$outcome)), known)
  d$outcome_f <- factor(ifelse(is.na(d$outcome), "Effect", d$outcome),
                        levels = c(known, sort(extra), "Effect"))
  d$outcome_f <- droplevels(d$outcome_f)

  # One row order per panel; ordering by effect puts the strongest on top.
  if (isTRUE(sort_by_effect)) {
    ord <- stats::aggregate(estimate ~ clock_label, data = d, FUN = mean)
    d$clock_label <- factor(d$clock_label,
                            levels = ord$clock_label[order(ord$estimate)])
  } else {
    keep <- lookup$clock_label[lookup$clock_label %in% d$clock_label]
    d$clock_label <- factor(d$clock_label, levels = rev(keep))
  }

  # The unit goes in the strip, because the x axis is shared down a column and
  # a single axis title cannot describe months and log10 hazard ratios at once.
  unit_by_outcome <- .tage_outcome_unit(levels(d$outcome_f), units)
  pretty_outcome <- c(NormalizedAge = "Normalized age")[levels(d$outcome_f)]
  pretty_outcome[is.na(pretty_outcome)] <- levels(d$outcome_f)[is.na(pretty_outcome)]
  strip <- ifelse(levels(d$outcome_f) == "Effect",
                  "Effect",
                  paste0(pretty_outcome, "\n(", unit_by_outcome, ")"))
  d$outcome_panel <- factor(strip[as.integer(d$outcome_f)], levels = strip)

  colours <- TAGE_OUTCOME_COLORS[levels(d$outcome_f)]
  colours[is.na(colours)] <- TAGE_TEXT_SECONDARY
  names(colours) <- strip

  if (identical(subtitle, NA)) {
    subtitle <- sprintf(
      "point = effect, whiskers = %d%% CI \u00b7 filled = %s < %s within panel",
      round(conf_level * 100), p_adjust, format(sig_threshold)
    )
  }

  multi_group <- length(unique(d$group2)) > 1L
  aes_point <- if (multi_group) {
    ggplot2::aes(x = .data$estimate, y = .data$clock_label,
                 colour = .data$group2, fill = .data$group2, shape = .data$significant)
  } else {
    ggplot2::aes(x = .data$estimate, y = .data$clock_label,
                 colour = .data$outcome_panel, fill = .data$outcome_panel,
                 shape = .data$significant)
  }

  p <- ggplot2::ggplot(d) +
    ggplot2::geom_vline(xintercept = 0, colour = TAGE_AXIS_COLOR, linewidth = 0.4) +
    # geom_errorbar(orientation = "y") rather than the deprecated geom_errorbarh().
    ggplot2::geom_errorbar(
      ggplot2::aes(y = .data$clock_label, xmin = .data$ci_low, xmax = .data$ci_high,
                   group = .data$group2),
      orientation = "y", width = 0.22, colour = TAGE_TEXT_MUTED, linewidth = 0.45,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::geom_point(aes_point, size = 2.4, stroke = 0.8,
                        position = ggplot2::position_dodge(width = 0.5)) +
    # 21 is a filled circle, 1 a hollow one: significance is legible in print
    # and without colour.
    ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 21), guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = NULL)

  if (multi_group) {
    group_cols <- tage_series_colors(unique(as.character(d$group2)))
    p <- p + ggplot2::labs(colour = "vs reference", fill = "vs reference") +
      ggplot2::scale_colour_manual(values = group_cols) +
      ggplot2::scale_fill_manual(values = group_cols)
  } else {
    p <- p + ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::scale_fill_manual(values = colours, guide = "none")
  }

  # Outcome has to sit on the columns: facet_grid frees the x scale per column
  # only, and months cannot share an axis with log10 hazard ratios.
  if (!is.null(split_by)) {
    p <- p + ggplot2::facet_grid(split ~ outcome_panel, scales = "free_x", drop = TRUE)
  } else {
    p <- p + ggplot2::facet_wrap(~outcome_panel, scales = "free_x", nrow = 1)
  }

  ref_lab <- unique(d$group1)[1]
  p <- p + ggplot2::labs(x = sprintf("effect vs %s", ref_lab)) +
    theme_tage(base_size = base_size, grid = "x") +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = base_size - 2.5),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
      panel.spacing = ggplot2::unit(10, "pt")
    )

  attr(p, "tage_stats") <- stats

  # Height follows the tallest panel, width the number of strata.
  n_rows <- length(levels(d$clock_label))
  n_cols <- if (is.null(split_by)) 1L else length(unique(d$split))
  n_panels <- length(levels(d$outcome_panel))
  if (is.null(height)) {
    height <- max(2.1, row_height * max(n_rows, 3) + 1.1) * n_panels
  }
  if (is.null(width)) width <- panel_width * max(n_cols, 1)
  p <- .tage_set_size(p, width, height)
  p
}

# ---------------------------------------------------------------------------
# Module heatmap
# ---------------------------------------------------------------------------

#' Heatmap of module-clock effects
#'
#' Rows are co-expression modules, columns are strata (or the comparison when
#' there is no stratification), and each cell carries the effect with a star
#' when the adjusted p-value clears \code{sig_threshold}. Modules live on
#' different native scales, so the fill is the effect divided by a robust scale
#' of its column while the printed number stays in the original units.
#'
#' The statistics come from \code{\link{tage_module_stats}}.
#'
#' @param data Data frame with one column per module clock.
#' @param module_columns Module-clock columns.
#' @param group_column Column holding the experimental groups.
#' @param reference_group Reference level.
#' @param compare_groups Levels compared against the reference; each becomes a
#'   facet. \code{NULL} uses every other level.
#' @param covariates,split_by,variance_strata,standardize,p_adjust,p_adjust_scope
#'   Passed to \code{\link{tage_module_stats}}.
#' @param sig_threshold Adjusted p-value below which a cell is starred.
#' @param modules_version Module set used for the row labels, one of
#'   \code{"4.6"}, \code{"5.4"} or \code{"human"}.
#' @param module_functions Named character vector overriding the bundled
#'   module-to-function annotation.
#' @param color_scale \code{"robust"} divides each column by its
#'   \code{robust_pct} percentile of \code{|effect|} so differently scaled
#'   columns share one colour bar; \code{"absolute"} uses the raw effect.
#' @param robust_pct,limit Percentile used by the robust scale, and the fill
#'   limit.
#' @param annotate Whether to print the effect in each cell.
#' @param digits Digits used for the printed effect.
#' @param title,subtitle,caption Figure text.
#' @param base_size Base font size.
#' @param cell_height,cell_width Inches per module row and per column, used to
#'   size the canvas so that row labels stay legible as the module set grows.
#' @param width,height Absolute figure size in inches, overriding
#'   \code{cell_height} and \code{cell_width}.
#' @param stats Pre-computed \code{\link{tage_module_stats}} output.
#'
#' @return A \code{ggplot} object, with the statistics attached as the
#'   \code{"tage_stats"} attribute and the intended size as
#'   \code{"tage_size"}; see \code{\link{tage_save_plot}}.
#'
#' @seealso \code{\link{tage_module_stats}}, \code{\link{tage_clock_forest}}
#'
#' @examples
#' \dontrun{
#' tage_module_heatmap(results, module_columns = modules,
#'                     group_column = "Genotype", reference_group = "WT",
#'                     split_by = "Tissue")
#' }
#' @export
tage_module_heatmap <- function(data,
                                module_columns,
                                group_column,
                                reference_group,
                                compare_groups = NULL,
                                covariates = NULL,
                                split_by = NULL,
                                variance_strata = c("all_data", "subset"),
                                standardize = TRUE,
                                p_adjust = "BH",
                                p_adjust_scope = "across_columns",
                                sig_threshold = 0.05,
                                modules_version = c("4.6", "5.4", "human"),
                                module_functions = NULL,
                                color_scale = c("robust", "absolute"),
                                robust_pct = 95,
                                limit = NULL,
                                annotate = TRUE,
                                digits = 2,
                                title = "Module clocks",
                                subtitle = NA,
                                caption = NULL,
                                base_size = 11,
                                cell_height = 0.30,
                                cell_width = 1.3,
                                width = NULL,
                                height = NULL,
                                stats = NULL) {
  .tage_gg_require()
  variance_strata <- match.arg(variance_strata)
  modules_version <- match.arg(modules_version)
  color_scale <- match.arg(color_scale)

  module_columns <- intersect(as.character(module_columns), colnames(data))
  if (length(module_columns) == 0L) {
    stop("None of the requested module columns are present in `data`.", call. = FALSE)
  }

  if (is.null(stats)) {
    stats <- tage_module_stats(
      data = data, module_columns = module_columns, group_column = group_column,
      reference_group = reference_group, compare_groups = compare_groups,
      covariates = covariates, split_by = split_by,
      variance_strata = variance_strata, standardize = standardize,
      p_adjust = p_adjust, p_adjust_scope = p_adjust_scope
    )
  }
  if (nrow(stats) == 0L) {
    stop("No module contrasts to plot - check the group labels and sample sizes.",
         call. = FALSE)
  }

  d <- stats
  d$column <- if (is.null(split_by)) as.character(d$group2) else as.character(d$split)
  # Name the facet after the contrast it shows, not just the treated group.
  d$panel <- factor(paste(d$group2, "vs", d$group1),
                    levels = unique(paste(d$group2, "vs", d$group1)))

  # Fill: normalise per column so columns with different units stay comparable.
  if (color_scale == "robust") {
    denom <- stats::aggregate(list(denom = abs(d$estimate)), by = list(column = d$column),
                              FUN = function(x) max(stats::quantile(x, robust_pct / 100,
                                                                    na.rm = TRUE), 1e-9))
    d <- merge(d, denom, by = "column", all.x = TRUE, sort = FALSE)
    d$fill <- d$estimate / d$denom
    lim <- if (is.null(limit)) 1.2 else limit
    fill_lab <- sprintf("effect /\n%gth pct |effect|\nof the column", robust_pct)
  } else {
    d$fill <- d$estimate
    lim <- if (is.null(limit)) max(abs(d$estimate), na.rm = TRUE) else limit
    fill_lab <- "effect"
  }
  d$fill <- pmax(pmin(d$fill, lim), -lim)

  funcs <- if (is.null(module_functions)) load_module_functions(modules_version) else module_functions
  pretty <- vapply(module_columns, function(m) {
    if (!is.null(funcs) && m %in% names(funcs)) paste0(m, " \u2014 ", funcs[[m]]) else m
  }, character(1))
  d$module_label <- factor(unname(pretty[as.character(d$module)]),
                           levels = rev(unname(pretty[module_columns])))

  d$star <- ifelse(!is.na(d$p_adjusted) & d$p_adjusted < sig_threshold, "*", "")
  d$text <- paste0(formatC(d$estimate, format = "f", digits = digits), d$star)

  if (identical(subtitle, NA)) {
    unit <- if (isTRUE(standardize)) "standardised mean difference" else "effect in native units"
    subtitle <- sprintf("cell values are the %s; * = %s < %s across modules",
                        unit, p_adjust, format(sig_threshold))
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$column, y = .data$module_label,
                                       fill = .data$fill)) +
    ggplot2::geom_tile(colour = TAGE_SURFACE, linewidth = 1) +
    # Column labels wrap on spaces so narrow cells do not collide.
    ggplot2::scale_x_discrete(labels = function(x) gsub(" ", "\n", x))

  if (isTRUE(annotate)) {
    # Ink on light cells, white on saturated ones; the ramp is symmetric so
    # the threshold is on |fill|.
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$text,
                   colour = abs(.data$fill) > 0.55 * lim),
      size = base_size / 4.2, show.legend = FALSE
    ) + ggplot2::scale_colour_manual(values = c(`FALSE` = TAGE_TEXT_PRIMARY, `TRUE` = TAGE_SURFACE),
                                     guide = "none")
  }

  p <- p +
    scale_fill_tage_diverging(limit = lim, name = fill_lab) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                  x = NULL, y = NULL) +
    theme_tage(base_size = base_size, grid = "none") +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = base_size - 2),
      legend.key.height = ggplot2::unit(0.9, "cm"),
      legend.key.width = ggplot2::unit(0.3, "cm")
    )

  attr(p, "tage_stats") <- stats

  n_modules <- length(levels(d$module_label))
  n_cols <- length(unique(d$column))
  n_panels <- length(levels(d$panel))
  if (is.null(height)) height <- max(2.6, cell_height * n_modules + 1.9)
  if (is.null(width)) {
    width <- max(3.6, cell_width * n_cols + 4.6) * max(n_panels, 1)
  }
  p <- .tage_set_size(p, width, height)
  p
}
