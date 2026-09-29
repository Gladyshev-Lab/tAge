#' Flag outlier samples
#'
#' Flags samples within each stratum; nothing is removed. The rule is the one
#' of the Python package's \code{tage.pp.flag_outliers}, and both give the
#' same numbers.
#'
#' \code{method = "distance"} (default) scores every sample by \code{1 - r},
#' the Pearson correlation of its log2(CPM + 1) profile with the median
#' profile of its stratum, and flags it when the robust z-score of that
#' distance (median and MAD within the stratum) exceeds \code{threshold}. It is
#' the paper's "correlation with the median profile of the tissue" rule with a
#' threshold relative to the stratum instead of a fixed rho < 0.5, which passes
#' even a sample from another tissue. Only samples further from the median than
#' the rest are flagged.
#'
#' \code{method = "pca_iqr"} is the paper's rule: PC1 or PC2 more than
#' \code{iqr_factor} interquartile ranges beyond the quartiles. The paper
#' applied it to hundreds of samples; with a dozen, the first components pick
#' out the noisiest sample and 1.5 IQR flags a normal sample in most strata.
#'
#' Both work on log2(CPM + 1) of the genes with a mean CPM of at least
#' \code{min_cpm} in the stratum. Depth only moves a sample's log profile by a
#' constant, which the correlation ignores, but low-count genes add noise that
#' grows as depth falls, and PCA on unnormalised counts mostly finds depth.
#'
#' Remove the flagged samples explicitly, after looking at them:
#' \code{eset <- eset[, !eset$tage_outlier]}.
#'
#' @param eset ExpressionSet of raw counts (bulk or pseudobulk).
#' @param split_by Column of the phenoData defining the strata (tissue,
#'   dataset, cell type) that are scored separately, as in
#'   \code{\link{tAge_preprocessing}}. \code{NULL}: one stratum.
#' @param method \code{"distance"} or \code{"pca_iqr"}.
#' @param threshold Robust z-score above which \code{"distance"} flags a
#'   sample. Default 5.
#' @param iqr_factor Interquartile ranges beyond the quartiles for
#'   \code{"pca_iqr"}. Default 1.5.
#' @param min_cpm Genes with a lower mean CPM in the stratum are left out.
#'   Default 10.
#' @param min_samples Strata with fewer samples are not scored (a warning says
#'   which). Default 6.
#' @return The ExpressionSet with phenoData columns \code{tage_outlier}
#'   (logical) and, for \code{"distance"}, \code{tage_outlier_distance} and
#'   \code{tage_outlier_z}; for \code{"pca_iqr"}, \code{tage_outlier_pc1} and
#'   \code{tage_outlier_pc2}. Samples of unscored strata get \code{FALSE} and
#'   \code{NA}.
#' @examples
#' eset <- make_ExpressionSet(load_example_expression_data(), load_example_metadata(),
#'                            verbose = FALSE)
#' eset <- flag_outliers(eset, split_by = "Tissue")
#' Biobase::pData(eset)[eset$tage_outlier, c("Tissue", "tage_outlier_z")]
#' eset <- eset[, !eset$tage_outlier]
#' @export
flag_outliers <- function(eset,
                          split_by = NULL,
                          method = c("distance", "pca_iqr"),
                          threshold = 5,
                          iqr_factor = 1.5,
                          min_cpm = 10,
                          min_samples = 6) {
  method <- match.arg(method)
  pd <- Biobase::pData(eset)
  if (!is.null(split_by) && !split_by %in% colnames(pd)) {
    stop("`split_by` '", split_by, "' is not a column of the phenoData.", call. = FALSE)
  }
  counts <- Biobase::exprs(eset)
  if (anyNA(counts) || any(counts < 0)) {
    stop("Outlier flagging needs raw counts: no missing or negative values.", call. = FALSE)
  }

  strata <- if (is.null(split_by)) rep("all", ncol(eset)) else as.character(pd[[split_by]])
  flags <- rep(FALSE, ncol(eset))
  score_names <- if (method == "distance") c("tage_outlier_distance", "tage_outlier_z")
                 else c("tage_outlier_pc1", "tage_outlier_pc2")
  scores <- matrix(NA_real_, ncol(eset), 2, dimnames = list(NULL, score_names))

  for (level in unique(strata)) {
    idx <- which(strata == level)
    where <- if (is.null(split_by)) "the data" else sprintf("%s = %s", split_by, level)
    if (length(idx) < min_samples) {
      warning(sprintf("%s: %d samples, fewer than min_samples = %d; not scored.",
                      where, length(idx), min_samples), call. = FALSE)
      next
    }
    L <- .tage_log_cpm(counts[, idx, drop = FALSE], min_cpm)
    if (nrow(L) < 100) {
      warning(sprintf("%s: only %d genes with mean CPM >= %g; not scored.", where, nrow(L), min_cpm),
              call. = FALSE)
      next
    }
    if (method == "distance") {
      s <- .tage_distance_scores(L)
      scores[idx, ] <- cbind(s$distance, s$z)
      flags[idx] <- s$z > threshold
    } else {
      s <- .tage_pca_iqr(L, iqr_factor)
      scores[idx, ] <- s$pcs
      flags[idx] <- s$outlier
    }
  }

  eset$tage_outlier <- flags
  for (nm in score_names) Biobase::pData(eset)[[nm]] <- scores[, nm]
  eset
}

# log2(CPM + 1) of the genes with a mean CPM of at least `min_cpm`; genes x samples.
.tage_log_cpm <- function(counts, min_cpm) {
  lib <- colSums(counts)
  if (any(lib <= 0)) stop("A sample has no counts.", call. = FALSE)
  cpm <- sweep(counts, 2, lib, "/") * 1e6
  log2(cpm[rowMeans(cpm) >= min_cpm, , drop = FALSE] + 1)
}

# 1 - Pearson r with the stratum's median profile, and its robust z.
.tage_distance_scores <- function(L) {
  median_profile <- apply(L, 1, stats::median)
  distance <- 1 - as.numeric(stats::cor(L, median_profile))
  centre <- stats::median(distance)
  spread <- stats::mad(distance, center = centre)
  # half the samples or more sit at exactly the same distance
  if (spread == 0) spread <- 1.2533 * mean(abs(distance - centre))
  z <- if (spread > 0) (distance - centre) / spread else rep(0, length(distance))
  list(distance = distance, z = z)
}

# PC1 / PC2 scores of the gene-centred profiles, and the IQR rule on them.
.tage_pca_iqr <- function(L, iqr_factor) {
  pcs <- stats::prcomp(t(L), center = TRUE, scale. = FALSE, rank. = 2)$x[, 1:2, drop = FALSE]
  outlier <- rep(FALSE, nrow(pcs))
  for (j in 1:2) {
    q <- stats::quantile(pcs[, j], c(0.25, 0.75), names = FALSE)
    spread <- iqr_factor * (q[2] - q[1])
    outlier <- outlier | pcs[, j] < q[1] - spread | pcs[, j] > q[2] + spread
  }
  list(pcs = unname(pcs), outlier = outlier)
}

#' Not available: use flag_outliers()
#'
#' @param ... Ignored.
#' @return Does not return; raises an error pointing to \code{\link{flag_outliers}}.
#' @keywords internal
#' @export
remove_outliers <- function(...) {
  stop("remove_outliers() is not available; use eset <- flag_outliers(eset, split_by = ...), ",
       "which adds the phenoData column tage_outlier, and remove the flagged samples with ",
       "eset <- eset[, !eset$tage_outlier].", call. = FALSE)
}
