#' Plot density curves for ExpressionSet data
#'
#' This function creates density plots for expression data in an ExpressionSet object.
#' It can plot density curves for all samples, with optional log transformation and
#' customizable styling options.
#'
#' @param eset An ExpressionSet object containing expression data.
#' @param title Character string for the plot title. Default is "Density Plot".
#' @param log_transform Logical indicating whether to apply log2 transformation
#'   before plotting. Default is TRUE.
#' @param na_rm Logical indicating whether to remove NA values when computing densities.
#'   Default is TRUE.
#' @param width Numeric value for plot width in inches. Default is 8.
#' @param height Numeric value for plot height in inches. Default is 6.
#' @param error_message Character string to display if plotting fails. Default is
#'   "Error: No data available".
#' @param palette Character string specifying the color palette. See ?hcl.colors for
#'   available options. Default is "viridis".
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
  palette        = "viridis",   # see ?hcl.colors for palette options
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

    # Colors: use hcl.colors for better distinction
    cols <- if (n_samples == 1) "black" else grDevices::hcl.colors(n_samples, palette = palette)

    # Plot first sample
    d1 <- density(expr_data[, 1], na.rm = na_rm)
    plot(d1,
         main = title,
         xlab = "Expression values",
         ylab = "Density",
         col  = cols[1])

    # Add remaining samples
    if (n_samples > 1) {
      for (i in 2:n_samples) {
        lines(density(expr_data[, i], na.rm = na_rm), col = cols[i])
      }
      # Add legend with sample labels
      legend(legend_position, legend = labs, col = cols, lty = 1, cex = 0.8, bty = "n")
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
#' observations, and — when \code{p_threshold} is set — non-significant ones, are
#' dropped before plotting so the panel stays readable.
#'
#' @param data Data frame of per-sample values, typically the prediction table
#'   returned by \code{\link{predict_tAge}}.
#' @param x_var Column name defining the groups on the x axis.
#' @param y_var Column name holding the values to plot, e.g. a \code{*_tAge}
#'   column.
#' @param subgroup_var Optional column name used to facet the plot. Default
#'   \code{NULL}.
#' @param colors Optional named vector of fill colours, keyed by the levels of
#'   \code{x_var}.
#' @param point_size,point_alpha Size and opacity of the jittered points.
#' @param box_width Width of the boxes.
#' @param stat_method Test used for the pairwise comparisons, passed to
#'   \code{ggpubr::stat_compare_means}. Default \code{"t.test"}.
#' @param comparisons List of length-2 character vectors giving the pairs to
#'   test. Default \code{NULL} tests all pairs of \code{x_var} levels.
#' @param p_label Label style for the annotations, e.g. \code{"p.signif"} for
#'   stars or \code{"p.format"} for numeric p-values.
#' @param p_threshold If supplied, comparisons whose p-value is not below this
#'   threshold are removed from the plot. When faceting, a comparison is kept if
#'   it is significant in at least one facet.
#' @param min_group_n Minimum number of non-missing observations required in both
#'   groups for a comparison to be shown. When faceting, every facet must meet it.
#' @param font_size Base font size; also scales the annotation text.
#' @param theme_type ggplot2 theme to apply, e.g. \code{"classic"} or
#'   \code{"minimal"}.
#' @param title,xlab,ylab Plot title and axis labels.
#' @param legend_position Legend placement passed to \code{ggplot2::theme}.
#' @param y_center If given, the y axis is made symmetric around this value —
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
  stat_method = "t.test",
  comparisons = NULL,
  p_label = "p.signif",
  p_threshold = NULL,
  min_group_n = 2,
  font_size = 12,
  theme_type = "classic",
  title = NULL,
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
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("ggpubr required")

  library(ggplot2)
  library(ggpubr)

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

  # Filter comparisons to groups with enough observations
  if (is.null(comparisons)) {
    comparisons <- combn(x_levels, 2, simplify = FALSE)
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
  if (!is.null(p_threshold) && length(valid_comparisons) > 0) {
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

  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_boxplot(width = box_width, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = point_size, alpha = point_alpha) +
    scale_x_discrete(limits = x_levels, drop = FALSE)

  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors[x_levels], breaks = x_levels, drop = FALSE)
  }

  # Add stat comparisons only if there are valid ones
  if (length(valid_comparisons) > 0) {
    y_range_vals  <- range(data[[y_var]], na.rm = TRUE)
    y_range_size  <- diff(y_range_vals)
    y_start       <- y_range_vals[2] + y_range_size * 0.05
    n_comp        <- length(valid_comparisons)
    label_y_pos   <- y_start + y_range_size * 0.10 * seq(0, n_comp - 1)
    y_max_needed  <- max(label_y_pos) + y_range_size * 0.08

    p <- p + stat_compare_means(
      method      = stat_method,
      comparisons = valid_comparisons,
      label       = p_label,
      size        = font_size / 3,
      tip.length  = 0.01,
      label.y     = label_y_pos
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
  if (!is.null(y_center)) p <- p + geom_hline(yintercept = y_center, linetype = "dashed", linewidth = 0.5)

  if (!is.null(subgroup_var)) {
    scales_use <- if (force_fixed && identical(facet_scales, "free_y")) "fixed" else facet_scales
    p <- p + facet_wrap(vars(.data[[subgroup_var]]), scales = scales_use)
  }

  theme_func <- switch(theme_type,
    "classic" = theme_classic,
    "bw"      = theme_bw,
    "minimal" = theme_minimal,
    stop("Invalid theme_type")
  )

  p <- p + theme_func(base_size = font_size) +
    theme(
      axis.text = element_text(size = font_size, color = "black"),
      axis.title = element_text(size = font_size + 2, face = "bold"),
      legend.position = legend_position,
      strip.text = element_text(size = font_size, face = "bold")
    )

  if (!is.null(title)) p <- p + ggtitle(title)
  p <- p + xlab(ifelse(!is.null(xlab), xlab, x_var))
  p <- p + ylab(ifelse(!is.null(ylab), ylab, y_var))

  return(p)
}
