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


#' Create advanced boxplots for tAge analysis
#'
#' This function creates publication-ready boxplots with statistical comparisons,
#' customizable styling, and faceting options. It's designed for visualizing
#' transcriptomic age analysis results and other comparative data.
#'
#' @param data A data frame containing the data to plot.
#' @param x_var Character string specifying the column name for the x-axis variable
#'   (categorical grouping variable).
#' @param y_var Character string specifying the column name for the y-axis variable
#'   (continuous variable to plot).
#' @param subgroup_var Character string specifying an optional column name for faceting
#'   (e.g., "Tissue"). Default is NULL.
#' @param colors Named character vector for manual color mapping of groups. Names should
#'   correspond to group levels. Default is NULL (automatic colors).
#' @param point_size Numeric value for the size of jittered points. Default is 2.
#' @param point_alpha Numeric value for the transparency of points (0-1). Default is 0.7.
#' @param box_width Numeric value for the width of boxplots. Default is 0.5.
#' @param stat_method Character string specifying the statistical test method.
#'   Options: "t.test", "wilcox.test", etc. Default is "t.test".
#' @param comparisons List of group comparisons for statistical testing. If NULL,
#'   all pairwise comparisons are performed. Default is NULL.
#' @param p_label Character string specifying how to display p-values. Options:
#'   "p.signif" (stars), "p.format" (numeric). Default is "p.signif".
#' @param font_size Numeric value for the base font size. Default is 12.
#' @param theme_type Character string specifying the ggplot2 theme. Options:
#'   "classic", "bw", "minimal". Default is "classic".
#' @param title Character string for the plot title. Default is NULL.
#' @param xlab Character string for the x-axis label. Default is NULL (uses x_var).
#' @param ylab Character string for the y-axis label. Default is NULL (uses y_var).
#' @param legend_position Character string specifying legend position. Options:
#'   "none", "right", "bottom", "top", "left". Default is "right".
#' @param y_center Numeric value for a central reference line. Can also be used for
#'   symmetric scaling around this value. Default is NULL.
#' @param y_min Numeric value for fixed minimum y-axis limit. Default is NULL.
#' @param y_max Numeric value for fixed maximum y-axis limit. Default is NULL.
#' @param facet_scales Character string for facet scaling. Options: "fixed", "free",
#'   "free_y", "free_x". Default is "free_y".
#' @param x_order Character vector specifying the order of x-axis groups. Default is NULL.
#' @return A ggplot2 object.
#' @export
#' @examples
#' # Create a simple boxplot
#' library(ggplot2)
#' data(mtcars)
#' p <- tage_boxplot(mtcars, x_var = "cyl", y_var = "mpg", 
#'                   title = "Miles per Gallon by Cylinders")
#' 
#' # Create a boxplot with custom colors and statistical comparisons
#' p <- tage_boxplot(mtcars, x_var = "cyl", y_var = "mpg",
#'                   colors = c("4" = "blue", "6" = "green", "8" = "red"),
#'                   stat_method = "wilcox.test")
tage_boxplot <- function(
  data,
  x_var,
  y_var,
  subgroup_var = NULL,   # Optional: variable for faceting (e.g., "Tissue")
  colors = NULL,         # Optional: manual color mapping for groups
  point_size = 2,        # Size of jittered points
  point_alpha = 0.7,     # Transparency of points
  box_width = 0.5,       # Width of the boxplot
  stat_method = "t.test",# Statistical test (e.g., "t.test", "wilcox.test")
  comparisons = NULL,    # Optional: list of group comparisons
  p_label = "p.signif",  # How to display p-values: "p.signif" (stars) or "p.format" (numeric)
  font_size = 12,        # Base font size
  theme_type = "classic",# ggplot2 theme: "classic", "bw", "minimal"
  title = NULL,          # Optional: plot title
  xlab = NULL,           # Optional: x-axis label
  ylab = NULL,           # Optional: y-axis label
  legend_position = "right", # Legend position ("none", "right", "bottom", etc.)
  y_center = NULL,       # Value for a central reference line, can also be used for symmetric scaling
  y_min = NULL,          # Fixed minimum for y-axis
  y_max = NULL,          # Fixed maximum for y-axis
  facet_scales = "free_y", # Facet scaling: "fixed", "free", "free_y", "free_x"
  x_order = NULL         # Optional: specific order for x-axis groups
) {
  library(ggplot2)
  library(ggpubr)

  
  # --- validation ---
  if (!x_var %in% colnames(data)) stop("x_var not found in data")
  if (!y_var %in% colnames(data)) stop("y_var not found in data")
  if (!is.null(subgroup_var) && !subgroup_var %in% colnames(data)) stop("subgroup_var not found in data")
  if (!is.null(y_min) && !is.null(y_max) && (y_min >= y_max)) stop("y_min must be < y_max")

  # --- enforce x order as factor levels ---
  if (!is.null(x_order)) {
    missing <- setdiff(x_order, unique(data[[x_var]]))
    if (length(missing)) stop("x_order has values not present in x_var: ", paste(missing, collapse = ", "))
    data[[x_var]] <- factor(data[[x_var]], levels = x_order)
  } else {
    data[[x_var]] <- factor(data[[x_var]])
  }
  x_levels <- levels(data[[x_var]])

  # --- base plot ---
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_boxplot(width = box_width, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = point_size, alpha = point_alpha) +
    scale_x_discrete(limits = x_levels, drop = FALSE)

  # --- colors (respect order) ---
  if (!is.null(colors)) {
    if (!all(x_levels %in% names(colors))) {
      stop("colors must be a named vector covering all x groups: ", paste(x_levels, collapse = ", "))
    }
    p <- p + scale_fill_manual(values = colors[x_levels], breaks = x_levels, drop = FALSE)
  }

  # --- stats (default comparisons follow x order) ---
  if (is.null(comparisons)) {
    comparisons <- combn(x_levels, 2, simplify = FALSE)
  }
  p <- p + stat_compare_means(
    method = stat_method,
    comparisons = comparisons,
    label = p_label,
    size = font_size / 3
  )
  # --- y-axis scaling logic ---
  force_fixed_scales <- (!is.null(y_center) || !is.null(y_min) || !is.null(y_max))
  ylim_vals <- NULL
  
  if (!is.null(y_center) && is.null(y_min) && is.null(y_max)) {
    y_vals <- data[[y_var]]
    half_range <- max(abs(y_vals - y_center), na.rm = TRUE) * 1.05
    ylim_vals <- c(y_center - half_range, y_center + half_range)
  }

  if (!is.null(y_min) || !is.null(y_max)) {
    if (is.null(y_min)) y_min <- min(data[[y_var]], na.rm = TRUE)
    if (is.null(y_max)) y_max <- max(data[[y_var]], na.rm = TRUE)
    ylim_vals <- c(y_min, y_max)
  }

  if (!is.null(ylim_vals)) {
    p <- p + coord_cartesian(ylim = ylim_vals, clip = "off")
  }

  if (!is.null(y_center)) {
    p <- p + geom_hline(yintercept = y_center, linetype = "dashed", linewidth = 0.5)
  }

  # Faceting
  if (!is.null(subgroup_var)) {
    # If y scaling is forced, make all facets share the same y-axis unless overridden
    scales_to_use <- if (force_fixed_scales && identical(facet_scales, "free_y")) "fixed" else facet_scales
    p <- p + facet_wrap(vars(.data[[subgroup_var]]), scales = scales_to_use)
  }

  # Apply theme
  theme_func <- switch(theme_type,
                       "classic" = theme_classic,
                       "bw" = theme_bw,
                       "minimal" = theme_minimal,
                       stop("Invalid theme_type")
  )
  p <- p + theme_func(base_size = font_size) +
    theme(
      axis.text = element_text(size = font_size, color = "black"),
      axis.title = element_text(size = font_size + 2, face = "bold"),
      legend.position = legend_position,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = font_size, face = "bold")
    )

  # Add labels
  if (!is.null(title)) p <- p + ggtitle(title)
  if (!is.null(xlab)) p <- p + xlab(xlab) else p <- p + xlab(x_var)
  if (!is.null(ylab)) p <- p + ylab(ylab) else p <- p + ylab(y_var)

  return(p)
}
