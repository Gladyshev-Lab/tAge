#' @keywords internal
"_PACKAGE"

# Base-package functions used unqualified across the package. ggplot2 is
# imported whole because tage_boxplot() and the figures build grammar layers
# from many of its functions.
#' @import ggplot2
#' @importFrom stats median quantile prcomp mahalanobis qchisq setNames lm density var cov
#' @importFrom graphics plot lines legend par plot.new text
#' @importFrom methods is new
#' @importFrom utils read.csv head combn
NULL
