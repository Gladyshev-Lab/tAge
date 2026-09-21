#' Plot density curves for ExpressionSet data
#'
#' This function creates density plots for expression data in an ExpressionSet object.
#' It can plot density curves for all samples, with optional log transformation and
#' customizable styling options.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param title Character string for the plot title. Default is "Density Plot".
#' @param subtitle,caption Subtitle (secondary ink) and a provenance caption
#'   (muted ink) under the plot.
#' @param log_transform Logical indicating whether to apply log2 transformation
#'   before plotting. Default is TRUE.
#' @param na_rm Logical indicating whether to remove NA values when computing densities.
#'   Default is TRUE.
#' @param width Numeric value for plot width in inches. Default is 8.
#' @param height Numeric value for plot height in inches. Default is 6.
#' @param error_message Character string to display if plotting fails. Default is
#'   "Error: No data available".
#' @param palette \code{NULL} (default) draws every sample in muted grey, the
#'   shape of the distributions being the point; a \code{hcl.colors} palette
#'   name colours the samples and adds a legend when there are at most 20.
#' @param legend_position Character string specifying legend position. Options include
#'   "topright", "topleft", "bottomright", "bottomleft", etc. Default is "topright".
#' @return Invisibly returns NULL. Creates a density plot.
#' @export
#' @examples
#' # Load example data and create ExpressionSet
#' expr_data <- load_example_expression_data()
#' meta_data <- load_example_metadata()
#' eset <- make_ExpressionSet(expr_data, meta_data)
#' 
#' # Plot density curves
#' plot_eset_density(eset, title = "Expression Density", log_transform = TRUE)
plot_eset_density <- function(
  eset,
  title          = "Density Plot",
  log_transform  = TRUE,
  na_rm          = TRUE,
  width          = 8,
  height         = 6,
  error_message  = "Error: No data available",
  palette        = NULL,
  legend_position     = "topright"
) {
  # Set plot dimensions in Jupyter
  options(repr.plot.width = width, repr.plot.height = height)

  tryCatch({
    # Check input type
    if (!inherits(eset, "ExpressionSet")) stop("Input must be an ExpressionSet")

    # Extract expression data
    expr_data <- Biobase::exprs(eset)
    if (log_transform) expr_data <- log2(expr_data + 1)

    # Basic checks for empty data
    if (nrow(expr_data) == 0 || ncol(expr_data) == 0) stop("No data to plot")

    n_samples <- ncol(expr_data)

    # Get sample labels: prefer sampleNames(), fallback to colnames()
    labs <- Biobase::sampleNames(eset)
    if (is.null(labs) || length(labs) != n_samples || any(is.na(labs) | labs == ""))
      labs <- colnames(expr_data)
    if (is.null(labs)) labs <- paste0("Sample ", seq_len(n_samples))

    # Muted lines: the shape of the distributions is the point, not sample
    # identity. A named palette colours the samples instead.
    cols <- if (is.null(palette)) {
      rep(grDevices::adjustcolor(TAGE_TEXT_MUTED, alpha.f = 0.6), n_samples)
    } else if (n_samples == 1) TAGE_TEXT_PRIMARY else grDevices::hcl.colors(n_samples, palette = palette)

    dens <- lapply(seq_len(n_samples), function(i) density(expr_data[, i], na.rm = na_rm))
    op <- par(family = "sans", bty = "l", col.axis = TAGE_TEXT_SECONDARY, col.lab = TAGE_TEXT_SECONDARY,
              fg = TAGE_AXIS_COLOR, mgp = c(2.2, 0.6, 0), tcl = -0.25, mar = c(4, 4, 3, 1))
    on.exit(par(op), add = TRUE)
    plot(dens[[1]], main = "", xlab = "Expression values", ylab = "Density",
         col = cols[1], lwd = 1, xlim = range(unlist(lapply(dens, `[[`, "x"))),
         ylim = c(0, max(unlist(lapply(dens, `[[`, "y")))))
    grid(nx = NA, ny = NULL, col = TAGE_GRID_COLOR, lty = 1, lwd = 0.6)
    for (i in seq_len(n_samples)) lines(dens[[i]], col = cols[i], lwd = 1)
    title(main = title, adj = 0, font.main = 2, col.main = TAGE_TEXT_PRIMARY, cex.main = 1.05)
    mtext(sprintf("%d samples", n_samples), side = 3, adj = 0, line = 0.2,
          col = TAGE_TEXT_SECONDARY, cex = 0.85)

    if (n_samples > 1 && !is.null(palette) && n_samples <= 20) {
      legend(legend_position, legend = labs, col = cols, lty = 1, cex = 0.75, bty = "n",
             text.col = TAGE_TEXT_PRIMARY)
    }

    invisible(NULL)
  }, error = function(e) {
    # Plot an empty figure with error message
    op <- par(mar = c(1, 1, 1, 1)); on.exit(par(op), add = TRUE)
    plot.new()
    text(0.5, 0.5, paste(error_message, "\n", e$message), cex = 1.1, col = "red")
    cat("Error in plot_eset_density:", e$message, "\n")
    invisible(NULL)
  })
}

#' Box plot of tAge predictions with pairwise significance annotation
#'
#' Draws a box plot with jittered points for one prediction column, split by a
#' grouping variable and optionally faceted by a subgroup. Pairwise comparisons
#' are annotated with brackets; comparisons involving groups with too few
#' observations, and -- when \code{p_threshold} is set -- non-significant ones, are
#' dropped before plotting so the panel stays readable.
#'
#' @param data Data frame of per-sample values, typically the prediction table
#'   returned by \code{\link{predict_tAge}}.
#' @param x_var Column name defining the groups on the x axis.
#' @param y_var Column name holding the values to plot, e.g. a \code{*_tAge}
#'   column.
#' @param subgroup_var Optional column name used to facet the plot. Default
#'   \code{NULL}.
#' @param colors Optional colours keyed by the levels of \code{x_var}, overriding
#'   the default in which the reference group is grey and the other groups take
#'   the categorical slots of \code{\link{tage_series_colors}} in order.
#' @param point_size,point_alpha Size and opacity of the jittered points.
#' @param box_width Width of the boxes.
#' @param stat_method Test used for the comparisons. \code{"emmeans"} (default)
#'   uses \code{\link{tage_compare_groups}}: estimated marginal-mean
#'   contrasts, which support covariates, stratum-wise models and Bayesian
#'   ridge weighting. Any other value is passed to
#'   \code{ggpubr::stat_compare_means}, e.g.
#'   \code{"t.test"} or \code{"wilcox.test"}.
#' @param reference_group Reference level for \code{stat_method = "emmeans"}.
#'   Default \code{NULL} uses the first level of \code{x_var}.
#' @param covariates Covariate columns adjusted for when
#'   \code{stat_method = "emmeans"}. Supplying them also switches the plotted
#'   values to covariate-adjusted ones (\code{\link{tage_adjust_covariates}}).
#' @param se_column Column of per-sample prediction standard deviations for a
#'   Bayesian ridge clock, enabling the meta-regression test. Default
#'   \code{NULL}.
#' @param variance_strata Passed to \code{\link{tage_compare_groups}}.
#' @param p_adjust,p_adjust_scope Multiplicity correction passed to
#'   \code{\link{tage_compare_groups}}.
#' @param comparisons List of length-2 character vectors giving the pairs to
#'   test. Default \code{NULL} compares every level against
#'   \code{reference_group} for \code{stat_method = "emmeans"}, and tests all
#'   pairs of \code{x_var} levels otherwise.
#' @param p_label Label style for the annotations, e.g. \code{"p.signif"} for
#'   stars or \code{"p.format"} for numeric p-values.
#' @param p_threshold If supplied, comparisons whose p-value is not below this
#'   threshold are removed from the plot. When faceting, a comparison is kept if
#'   it is significant in at least one facet.
#' @param min_group_n Minimum number of non-missing observations required in both
#'   groups for a comparison to be shown. When faceting, every facet must meet it.
#' @param font_size Base font size; also scales the annotation text.
#' @param theme_type Deprecated; figures follow \code{\link{theme_tage}}.
#' @param title,xlab,ylab Plot title and axis labels.
#' @param legend_position Legend placement passed to \code{ggplot2::theme}.
#' @param y_center If given, the y axis is made symmetric around this value --
#'   useful for relative predictions centred on zero.
#' @param y_min,y_max Explicit y-axis limits, overriding the automatic range.
#' @param facet_scales Scale behaviour across facets, passed to
#'   \code{ggplot2::facet_wrap}. Default \code{"free_y"}.
#' @param x_order Character vector giving the order of the \code{x_var} levels.
#' @param width,height Plot size in inches, used to set the inline display size.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{predict_tAge}} for producing the input table.
#'
#' @examples
#' \dontrun{
#' results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "EN")
#' tage_boxplot(
#'   results,
#'   x_var = "Genotype",
#'   y_var = "yugene_diff_EN_tAge",
#'   subgroup_var = "Tissue",
#'   x_order = c("WT", "Klotho KO"),
#'   ylab = "Relative tAge, months"
#' )
#' }
#' @export
tage_boxplot <- function(
  data,
  x_var,
  y_var,
  subgroup_var = NULL,
  colors = NULL,
  point_size = 2,
  point_alpha = 0.7,
  box_width = 0.5,
  stat_method = "emmeans",
  reference_group = NULL,
  covariates = NULL,
  se_column = NULL,
  variance_strata = c("subset", "all_data"),
  p_adjust = "BH",
  p_adjust_scope = "within_column",
  comparisons = NULL,
  p_label = "p.signif",
  p_threshold = NULL,
  min_group_n = 2,
  font_size = 10,
  theme_type = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  xlab = NULL,
  ylab = NULL,
  legend_position = "right",
  y_center = NULL,
  y_min = NULL,
  y_max = NULL,
  facet_scales = "free_y",
  x_order = NULL,
  width  = 10,
  height = 6
) {
  options(repr.plot.width = width, repr.plot.height = height)

  if (!x_var %in% colnames(data)) stop("x_var not found in data")
  if (!y_var %in% colnames(data)) stop("y_var not found in data")
  if (!is.null(subgroup_var) && !subgroup_var %in% colnames(data)) stop("subgroup_var not found")

  if (!is.null(x_order)) {
    data[[x_var]] <- factor(data[[x_var]], levels = x_order)
  } else {
    data[[x_var]] <- factor(data[[x_var]])
  }
  x_levels <- levels(data[[x_var]])

  variance_strata <- match.arg(variance_strata)
  use_emmeans <- identical(stat_method, "emmeans")

  stat_table <- NULL
  if (use_emmeans) {
    if (is.null(reference_group)) reference_group <- x_levels[1]
    if (!reference_group %in% x_levels) {
      stop("reference_group '", reference_group, "' is not a level of x_var")
    }

    stat_table <- tage_compare_groups(
      data            = data,
      value_columns   = y_var,
      group_column    = x_var,
      reference_group = reference_group,
      covariates      = covariates,
      split_by        = subgroup_var,
      se_columns      = se_column,
      method          = "trt.vs.ctrl",
      variance_strata = variance_strata,
      p_adjust        = p_adjust,
      p_adjust_scope  = p_adjust_scope
    )

    # Plot what the model tested: with covariates, partial residuals rather
    # than raw predictions.
    if (!is.null(covariates)) {
      data[[".tage_adjusted"]] <- tage_adjust_covariates(
        data, y_var, covariates, split_by = subgroup_var, se_column = se_column
      )
      y_var <- ".tage_adjusted"
    }
  }

  # Filter comparisons to groups with enough observations
  if (is.null(comparisons)) {
    comparisons <- if (use_emmeans) {
      lapply(setdiff(x_levels, reference_group), function(g) c(reference_group, g))
    } else {
      combn(x_levels, 2, simplify = FALSE)
    }
  }

  # If faceting, filter comparisons per facet; otherwise filter globally
  if (!is.null(subgroup_var)) {
    # Keep comparison only if EVERY facet group has >= min_group_n for both levels
    valid_comparisons <- Filter(function(comp) {
      all(sapply(unique(data[[subgroup_var]]), function(grp) {
        sub <- data[data[[subgroup_var]] == grp, ]
        n1 <- sum(!is.na(sub[[y_var]]) & sub[[x_var]] == comp[1])
        n2 <- sum(!is.na(sub[[y_var]]) & sub[[x_var]] == comp[2])
        n1 >= min_group_n && n2 >= min_group_n
      }))
    }, comparisons)
  } else {
    valid_comparisons <- Filter(function(comp) {
      n1 <- sum(!is.na(data[[y_var]]) & data[[x_var]] == comp[1])
      n2 <- sum(!is.na(data[[y_var]]) & data[[x_var]] == comp[2])
      n1 >= min_group_n && n2 >= min_group_n
    }, comparisons)
  }

  # If p_threshold set, pre-compute and keep only significant comparisons
  if (!is.null(p_threshold) && length(valid_comparisons) > 0 && use_emmeans) {
    # Keep a comparison if it clears the threshold in at least one stratum.
    valid_comparisons <- Filter(function(comp) {
      hit <- stat_table$group1 == comp[1] & stat_table$group2 == comp[2]
      any(hit) && any(stat_table$p_adjusted[hit] < p_threshold, na.rm = TRUE)
    }, valid_comparisons)
  } else if (!is.null(p_threshold) && length(valid_comparisons) > 0) {
    if (!is.null(subgroup_var)) {
      # Keep if significant in at least one facet
      valid_comparisons <- Filter(function(comp) {
        any(sapply(unique(data[[subgroup_var]]), function(grp) {
          sub <- data[data[[subgroup_var]] == grp, ]
          v1 <- sub[[y_var]][sub[[x_var]] == comp[1]]
          v2 <- sub[[y_var]][sub[[x_var]] == comp[2]]
          if (length(v1) < min_group_n || length(v2) < min_group_n) return(FALSE)
          tryCatch({
            test <- do.call(stat_method, list(x = v1, y = v2))
            test$p.value < p_threshold
          }, error = function(e) FALSE)
        }))
      }, valid_comparisons)
    } else {
      valid_comparisons <- Filter(function(comp) {
        v1 <- data[[y_var]][data[[x_var]] == comp[1]]
        v2 <- data[[y_var]][data[[x_var]] == comp[2]]
        tryCatch({
          test <- do.call(stat_method, list(x = v1, y = v2))
          test$p.value < p_threshold
        }, error = function(e) FALSE)
      }, valid_comparisons)
    }
  }

  if (!is.null(theme_type)) {
    warning("`theme_type` is deprecated and has no effect: figures follow theme_tage().",
            call. = FALSE)
  }
  # Reference grey, comparisons in the fixed categorical slots; `colors`
  # overrides per level.
  ref_for_colour <- if (!is.null(reference_group)) reference_group else x_levels[1]
  cols <- tage_series_colors(x_levels, reference = ref_for_colour, palette = colors)

  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]],
                        fill = .data[[x_var]], colour = .data[[x_var]])) +
    # A coloured edge over a tinted fill (alpha over white), not a solid block.
    geom_boxplot(width = box_width, outlier.shape = NA, alpha = 0.35,
                 linewidth = 0.5, fatten = 2) +
    geom_point(position = position_jitter(width = 0.15, seed = 1), shape = 21,
               size = point_size, alpha = point_alpha, colour = TAGE_SURFACE, stroke = 0.5) +
    scale_x_discrete(limits = x_levels, drop = FALSE) +
    scale_fill_manual(values = cols, breaks = x_levels, drop = FALSE, name = x_var) +
    scale_colour_manual(values = cols, breaks = x_levels, drop = FALSE, guide = "none")

  # Add stat comparisons only if there are valid ones
  if (length(valid_comparisons) > 0 && use_emmeans) {
    keep <- mapply(
      function(g1, g2) any(vapply(valid_comparisons,
                                  function(cp) cp[1] == g1 && cp[2] == g2, logical(1))),
      stat_table$group1, stat_table$group2
    )
    ann <- stat_table[keep, , drop = FALSE]

    y_range_vals <- range(data[[y_var]], na.rm = TRUE)
    y_range_size <- diff(y_range_vals)

    if (nrow(ann) > 0) {
      # Stack brackets per stratum, above that stratum's own data.
      step <- 0.10 * y_range_size
      ann$y.position <- NA_real_

      if (is.null(subgroup_var)) {
        ann$y.position <- y_range_vals[2] + step * seq_len(nrow(ann))
      } else {
        tops <- tapply(data[[y_var]], as.character(data[[subgroup_var]]),
                       max, na.rm = TRUE)
        for (k in unique(as.character(ann$split))) {
          idx <- which(as.character(ann$split) == k)
          base_top <- if (k %in% names(tops)) tops[[k]] else y_range_vals[2]
          ann$y.position[idx] <- base_top + step * seq_along(idx)
        }
      }

      ann$label <- switch(
        p_label,
        "p.signif" = ifelse(nzchar(ann$label), ann$label, "ns"),
        "p.format" = format.pval(ann$p_adjusted, digits = 2, eps = 1e-16),
        format.pval(ann$p_adjusted, digits = 2, eps = 1e-16)
      )
      if (!is.null(subgroup_var)) ann[[subgroup_var]] <- ann$split

      y_max_needed <- max(ann$y.position, na.rm = TRUE) + y_range_size * 0.08

      .tage_require("ggpubr")
      p <- p + ggpubr::stat_pvalue_manual(
        ann,
        label      = "label",
        y.position = "y.position",
        xmin       = "group1",
        xmax       = "group2",
        tip.length = 0.01,
        size       = font_size / 3,
        color      = TAGE_TEXT_SECONDARY,
        bracket.size = 0.4
      )
    } else {
      y_max_needed <- y_range_vals[2] + y_range_size * 0.05
    }
  } else if (length(valid_comparisons) > 0) {
    y_range_vals  <- range(data[[y_var]], na.rm = TRUE)
    y_range_size  <- diff(y_range_vals)
    y_start       <- y_range_vals[2] + y_range_size * 0.05
    n_comp        <- length(valid_comparisons)
    label_y_pos   <- y_start + y_range_size * 0.10 * seq(0, n_comp - 1)
    y_max_needed  <- max(label_y_pos) + y_range_size * 0.08

    .tage_require("ggpubr")
    p <- p + ggpubr::stat_compare_means(
      method      = stat_method,
      comparisons = valid_comparisons,
      label       = p_label,
      size        = font_size / 3,
      tip.length  = 0.01,
      label.y     = label_y_pos,
      color       = TAGE_TEXT_SECONDARY,
      bracket.size = 0.4
    )
  } else {
    y_range_vals <- range(data[[y_var]], na.rm = TRUE)
    y_range_size <- diff(y_range_vals)
    y_max_needed <- y_range_vals[2] + y_range_size * 0.05
  }

  # Y-axis limits
  force_fixed <- (!is.null(y_center) || !is.null(y_min) || !is.null(y_max))
  ylim_vals <- NULL

  if (!is.null(y_center) && is.null(y_min) && is.null(y_max)) {
    half_range <- max(abs(data[[y_var]] - y_center), na.rm = TRUE) * 1.05
    half_range <- max(half_range, abs(y_max_needed - y_center))
    ylim_vals  <- c(y_center - half_range, y_center + half_range)
  } else if (!is.null(y_min) || !is.null(y_max)) {
    ylim_vals <- c(
      ifelse(is.null(y_min), y_range_vals[1], y_min),
      ifelse(is.null(y_max), y_range_vals[2], y_max)
    )
  } else {
    ylim_vals <- c(y_range_vals[1] - y_range_size * 0.05, y_max_needed)
  }

  if (!is.null(ylim_vals)) p <- p + coord_cartesian(ylim = ylim_vals, clip = "off")
  if (!is.null(y_center)) p <- p + geom_hline(yintercept = y_center, colour = TAGE_AXIS_COLOR, linewidth = 0.4)

  if (!is.null(subgroup_var)) {
    scales_use <- if (force_fixed && identical(facet_scales, "free_y")) "fixed" else facet_scales
    p <- p + facet_wrap(vars(.data[[subgroup_var]]), scales = scales_use)
  }

  p <- p + theme_tage(base_size = font_size, grid = "y") +
    theme(legend.position = legend_position)

  p <- p + labs(title = title, subtitle = subtitle, caption = caption,
                x = if (!is.null(xlab)) xlab else x_var,
                y = if (!is.null(ylab)) ylab else y_var)

  return(p)
}
