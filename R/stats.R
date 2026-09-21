# Statistical tests for tAge predictions.
#
# Group comparisons are estimated marginal-mean contrasts from a linear model
# (elastic net clocks) or from a meta-regression
# weighted by the per-sample prediction standard deviation (Bayesian ridge
# clocks). Continuous predictors are single coefficients from the same models.
# Nothing here rescales tAge -- predict_tAge() has already put predictions on
# their reporting scale.

TAGE_P_ADJUST_METHODS <- c("BH", "fdr", "bonferroni", "holm", "hochberg",
                           "hommel", "BY", "none")

#' Significance stars for tAge statistics
#'
#' Maps p-values onto the label set of the published figures: \code{***} below 0.001, \code{**} below 0.01,
#' \code{*} below 0.05 and \code{^} below 0.1.
#'
#' @param p Numeric vector of p-values.
#'
#' @return Character vector of labels; \code{""} where nothing is significant
#'   and \code{NA} stays \code{NA}.
#'
#' @examples
#' tage_significance_stars(c(1e-4, 0.02, 0.08, 0.5))
#' @export
tage_significance_stars <- function(p) {
  p <- as.numeric(p)
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.1]   <- "^"
  out[!is.na(p) & p < 0.05]  <- "*"
  out[!is.na(p) & p < 0.01]  <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out[is.na(p)] <- NA_character_
  out
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.tage_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for tAge statistics. Install it with install.packages('%s').",
                 pkg, pkg), call. = FALSE)
  }
}

.tage_clean_chr <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "None"]
  if (length(x) == 0) NULL else unique(x)
}

# Resolve the per-value-column standard-deviation columns. Accepts NULL, a
# single name, a vector parallel to value_columns, or a named vector keyed by
# value column. Returns a character vector the same length as value_columns
# with NA where a column has no SD (elastic net clocks).
.tage_resolve_se <- function(se_columns, value_columns, data) {
  out <- rep(NA_character_, length(value_columns))
  if (is.null(se_columns)) return(out)

  se_columns <- as.character(se_columns)
  nms <- names(se_columns)

  if (!is.null(nms) && all(nzchar(nms))) {
    idx <- match(value_columns, nms)
    out <- unname(se_columns[idx])
  } else if (length(se_columns) == 1L) {
    out[] <- se_columns
  } else if (length(se_columns) == length(value_columns)) {
    out <- se_columns
  } else {
    stop("`se_columns` must be NULL, a single name, a vector parallel to `value_columns`, or a named vector.",
         call. = FALSE)
  }

  out[!is.na(out) & !nzchar(out)] <- NA_character_
  missing <- setdiff(stats::na.omit(out), colnames(data))
  if (length(missing)) {
    stop("SD column(s) not found in data: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out
}

.tage_rhs <- function(terms) {
  paste(terms, collapse = " + ")
}

# Split a data frame by a stratifying column. Always returns a list of
# list(label=, data=) so that an absent stratification (label NA) cannot be
# confused with a stratum literally named "NA".
.tage_split <- function(data, split_by) {
  if (is.null(split_by)) {
    return(list(list(label = NA_character_, data = data)))
  }
  parts <- split(data, factor(data[[split_by]]), drop = TRUE)
  Map(function(label, part) list(label = label, data = part), names(parts), parts)
}

# Contrast list in emmeans' own coefficient form. Building the weights
# explicitly avoids parsing contrast labels, which breaks whenever a factor
# level happens to start with the name of its column.
.tage_contrast_weights <- function(levels_present, reference, method) {
  levels_present <- as.character(levels_present)
  if (method == "pairwise") {
    pairs <- utils::combn(seq_along(levels_present), 2, simplify = FALSE)
    # Weights are built so the estimate is group2 - group1, matching the
    # trt.vs.ctrl direction; the reported pair is (levels[i], levels[j]).
    con <- lapply(pairs, function(ij) {
      w <- rep(0, length(levels_present))
      w[ij[2]] <- 1
      w[ij[1]] <- -1
      w
    })
    names(con) <- vapply(pairs, function(ij) {
      paste(levels_present[ij[1]], levels_present[ij[2]], sep = "\r")
    }, character(1))
    return(con)
  }

  ref_i <- which(levels_present == reference)
  if (length(ref_i) != 1L) return(NULL)
  others <- which(levels_present != reference)
  if (length(others) == 0L) return(NULL)

  con <- lapply(others, function(i) {
    w <- rep(0, length(levels_present))
    w[i] <- 1
    w[ref_i] <- -1
    w
  })
  names(con) <- paste(reference, levels_present[others], sep = "\r")
  con
}

# Two-sided critical value: normal for a meta-regression (df = Inf), t otherwise.
.tage_crit <- function(df, conf_level) {
  q <- 0.5 + conf_level / 2
  if (any(!is.finite(df))) stats::qnorm(q) else stats::qt(q, df)
}

# emmeans renames the test statistic depending on the degrees of freedom, so
# pull whichever of the two is present.
.tage_statistic <- function(cont_df) {
  if ("t.ratio" %in% names(cont_df)) return(cont_df$t.ratio)
  if ("z.ratio" %in% names(cont_df)) return(cont_df$z.ratio)
  rep(NA_real_, nrow(cont_df))
}

.tage_contrast_df <- function(emm, group_column, reference, method,
                              conf_level = 0.95) {
  levels_present <- as.character(emm@grid[[group_column]])
  con <- .tage_contrast_weights(levels_present, reference, method)
  if (is.null(con)) return(NULL)

  cont <- emmeans::contrast(emm, method = con, adjust = "none")
  cont_df <- as.data.frame(cont)

  keys <- strsplit(as.character(cont_df$contrast), "\r", fixed = TRUE)
  df_col <- if ("df" %in% names(cont_df)) cont_df$df else rep(NA_real_, nrow(cont_df))
  crit <- vapply(df_col, .tage_crit, numeric(1), conf_level = conf_level)

  data.frame(
    group1    = vapply(keys, `[`, character(1), 1L),
    group2    = vapply(keys, `[`, character(1), 2L),
    estimate  = cont_df$estimate,
    se        = cont_df$SE,
    ci_low    = cont_df$estimate - crit * cont_df$SE,
    ci_high   = cont_df$estimate + crit * cont_df$SE,
    statistic = .tage_statistic(cont_df),
    df        = df_col,
    p_value   = cont_df$p.value,
    stringsAsFactors = FALSE
  )
}

# A stratum that cannot be tested is reported, never dropped in silence: the
# caller sees which clock and stratum are missing from the table and why.
.tage_skip <- function(value_column, split, reason) {
  where <- if (is.null(split) || is.na(split)) value_column else sprintf("%s [%s]", value_column, split)
  warning(sprintf("Skipping %s: %s", where, reason), call. = FALSE)
  invisible(NULL)
}

.tage_try_reason <- function(x) {
  conditionMessage(attr(x, "condition"))
}

# Fit the group model and return emmeans for the grouping factor. Elastic net
# clocks go through lm(); Bayesian ridge clocks go through a REML
# meta-regression weighted by 1 / (sd^2 + tau^2), read back into emmeans with
# qdrg() on z-tests.
# Returns the emmeans object, or a character string saying why the fit failed.
.tage_fit_emm <- function(df, response, group_column, covariates, se_column) {
  rhs <- .tage_rhs(c(group_column, covariates))

  if (is.na(se_column)) {
    fml <- stats::as.formula(paste(response, "~", rhs))
    model <- try(stats::lm(fml, data = df), silent = TRUE)
    if (inherits(model, "try-error")) return(paste("lm failed:", .tage_try_reason(model)))
    if (any(is.na(stats::coef(model)))) {
      return("redundant predictor(s) dropped from the linear model - check for collinear covariates")
    }
    emm <- try(emmeans::emmeans(model, specs = group_column), silent = TRUE)
    if (inherits(emm, "try-error")) return(paste("emmeans failed:", .tage_try_reason(emm)))
    return(emm)
  }

  .tage_require("metafor")
  mods <- stats::as.formula(paste("~", rhs))
  model <- try(
    suppressWarnings(
      metafor::rma.uni(yi = df[[response]], sei = df[[se_column]],
                       mods = mods, data = df, method = "REML")
    ),
    silent = TRUE
  )
  if (inherits(model, "try-error")) return(paste("rma.uni failed:", .tage_try_reason(model)))

  qrg <- try(
    emmeans::qdrg(formula = mods, data = df,
                  coef = stats::coef(model), vcov = stats::vcov(model),
                  df = Inf),
    silent = TRUE
  )
  if (inherits(qrg, "try-error")) return(paste("qdrg failed:", .tage_try_reason(qrg)))

  emm <- try(emmeans::emmeans(qrg, specs = group_column), silent = TRUE)
  if (inherits(emm, "try-error")) return(paste("emmeans failed:", .tage_try_reason(emm)))
  emm
}

# Apply p.adjust over the scope requested. `key_cols` identify one adjustment
# family; rows sharing a key are adjusted together.
.tage_adjust <- function(res, p_adjust, scope, column_key = "value_column",
                         cell_keys = c("split", "group1", "group2")) {
  res$p_adjusted <- NA_real_
  if (nrow(res) == 0L) return(res)

  if (identical(scope, "none") || identical(p_adjust, "none")) {
    res$p_adjusted <- res$p_value
    return(res)
  }

  key <- switch(
    scope,
    within_column  = res[[column_key]],
    across_columns = do.call(paste, c(lapply(cell_keys, function(k) res[[k]]), sep = "\r")),
    global         = rep("all", nrow(res)),
    stop("Unsupported `p_adjust_scope`: ", scope, call. = FALSE)
  )
  key <- as.character(key)
  key[is.na(key)] <- "\rNA\r"

  for (k in unique(key)) {
    idx <- which(key == k)
    res$p_adjusted[idx] <- stats::p.adjust(res$p_value[idx], method = p_adjust)
  }
  res
}

.tage_finalise <- function(res, p_adjust, scope, column_key = "value_column",
                           cell_keys = c("split", "group1", "group2")) {
  res <- .tage_adjust(res, p_adjust, scope, column_key, cell_keys)
  res$label <- tage_significance_stars(res$p_adjusted)
  rownames(res) <- NULL
  res
}

# ---------------------------------------------------------------------------
# Group comparisons
# ---------------------------------------------------------------------------

#' Compare tAge between experimental groups
#'
#' Estimated marginal-mean contrasts of predicted tAge between a reference
#' group and one or more comparison groups, as used for the clock analyses in
#' the paper.
#'
#' For elastic net clocks the model is \code{value ~ group + covariates} fitted
#' with \code{\link[stats]{lm}}. For Bayesian ridge clocks -- signalled by
#' supplying \code{se_columns} -- the model is a REML meta-regression
#' (\code{\link[metafor]{rma.uni}}) that weights each sample by its own
#' prediction uncertainty, and the contrasts are z-tests rather than t-tests.
#' In both cases contrasts come from \code{\link[emmeans]{emmeans}} with no
#' built-in multiplicity adjustment; correction is applied afterwards over the
#' scope given by \code{p_adjust_scope}.
#'
#' @param data Data frame of per-sample predictions, typically the table
#'   returned by \code{\link{predict_tAge}}.
#' @param value_columns Character vector of columns holding tAge predictions.
#'   Several columns are accepted so that \code{p_adjust_scope =
#'   "across_columns"} can correct across clocks or modules.
#' @param group_column Column defining the experimental groups.
#' @param reference_group Level of \code{group_column} used as the reference.
#' @param compare_groups Levels compared against the reference. Default
#'   \code{NULL} uses every other level present.
#' @param covariates Character vector of covariate columns entering the model
#'   as fixed effects. Default \code{NULL}.
#' @param split_by Column to stratify on, fitting a separate model per level
#'   (e.g. \code{"Tissue"}). Default \code{NULL}.
#' @param se_columns Columns holding per-sample prediction standard deviations
#'   from Bayesian ridge clocks. \code{NULL} (default) fits ordinary linear
#'   models. May be a single name, a vector parallel to \code{value_columns},
#'   or a vector named by value column; \code{NA} entries fall back to
#'   \code{lm}.
#' @param method \code{"trt.vs.ctrl"} (default) contrasts every comparison
#'   group against the reference; \code{"pairwise"} contrasts all pairs.
#' @param variance_strata \code{"subset"} (default) estimates the residual
#'   variance from the compared groups only; \code{"all_data"} fits the model
#'   on every group present in the stratum and reports the requested contrasts
#'   from it.
#' @param p_adjust Method passed to \code{\link[stats]{p.adjust}}; \code{"BH"}
#'   by default.
#' @param p_adjust_scope Family over which p-values are corrected.
#'   \code{"within_column"} (default) corrects across all comparisons and
#'   strata of one clock. \code{"across_columns"} corrects across clocks within
#'   each comparison and stratum (the family used for module clocks).
#'   \code{"global"} corrects everything together, \code{"none"} disables it.
#' @param conf_level Two-sided confidence level for \code{ci_low} /
#'   \code{ci_high}. The critical value follows the test: normal for the
#'   meta-regression, t otherwise.
#'
#' @return A data frame with one row per contrast and the columns
#'   \code{value_column}, \code{split}, \code{group1} (reference),
#'   \code{group2}, \code{n}, \code{estimate} (always \code{group2 - group1}),
#'   \code{se}, \code{statistic}, \code{df}, \code{p_value}, \code{p_adjusted}
#'   and \code{label}.
#'
#' @seealso \code{\link{tage_regress_continuous}} for numeric predictors,
#'   \code{\link{tage_module_stats}} for module clocks, and
#'   \code{\link{tage_adjust_covariates}} for the matching plotting values.
#'
#' @examples
#' \dontrun{
#' results <- predict_tAge(tAge_eset, model_paths, species = "mouse", mode = "EN")
#' tage_compare_groups(
#'   results,
#'   value_columns   = "yugene_diff_EN_tAge",
#'   group_column    = "Genotype",
#'   reference_group = "WT",
#'   split_by        = "Tissue"
#' )
#' }
#' @export
tage_compare_groups <- function(data,
                                value_columns,
                                group_column,
                                reference_group,
                                compare_groups = NULL,
                                covariates = NULL,
                                split_by = NULL,
                                se_columns = NULL,
                                method = c("trt.vs.ctrl", "pairwise"),
                                variance_strata = c("subset", "all_data"),
                                p_adjust = TAGE_P_ADJUST_METHODS,
                                p_adjust_scope = c("within_column", "across_columns",
                                                   "global", "none"),
                                conf_level = 0.95) {
  .tage_require("emmeans")

  method          <- match.arg(method)
  variance_strata <- match.arg(variance_strata)
  p_adjust        <- match.arg(p_adjust, TAGE_P_ADJUST_METHODS)
  p_adjust_scope  <- match.arg(p_adjust_scope)

  data          <- as.data.frame(data)
  value_columns <- .tage_clean_chr(value_columns)
  covariates    <- .tage_clean_chr(covariates)
  split_by      <- .tage_clean_chr(split_by)

  if (is.null(value_columns)) stop("`value_columns` must name at least one column.", call. = FALSE)
  missing <- setdiff(c(value_columns, group_column, covariates, split_by), colnames(data))
  if (length(missing)) stop("Column(s) not found in data: ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(split_by) > 1L) stop("`split_by` must name a single column.", call. = FALSE)

  se_map <- .tage_resolve_se(se_columns, value_columns, data)

  data[[group_column]] <- as.factor(data[[group_column]])
  present <- levels(droplevels(data[[group_column]]))
  reference_group <- as.character(reference_group)
  if (!reference_group %in% present) {
    stop("`reference_group` '", reference_group, "' is not present in ", group_column, ".", call. = FALSE)
  }

  if (is.null(compare_groups)) {
    compare_groups <- setdiff(present, reference_group)
  } else {
    compare_groups <- setdiff(as.character(compare_groups), reference_group)
    unknown <- setdiff(compare_groups, present)
    if (length(unknown)) {
      stop("`compare_groups` not present in ", group_column, ": ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
  }
  if (length(compare_groups) == 0L) stop("No comparison groups left to test.", call. = FALSE)

  # Reference first so the model's baseline matches the reported direction.
  data[[group_column]] <- factor(
    data[[group_column]],
    levels = c(reference_group, compare_groups,
               setdiff(present, c(reference_group, compare_groups)))
  )

  requested <- c(reference_group, compare_groups)
  rows <- list()

  for (i in seq_along(value_columns)) {
    vc <- value_columns[i]
    se_col <- se_map[i]

    keep_cols <- c(vc, group_column, covariates, split_by)
    if (!is.na(se_col)) keep_cols <- c(keep_cols, se_col)
    work <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]), , drop = FALSE]
    if (nrow(work) == 0L) {
      .tage_skip(vc, NA, "no complete cases (missing values in the value, group, covariate or split columns)")
      next
    }

    for (stratum in .tage_split(work, split_by)) {
      df <- stratum$data

      if (variance_strata == "subset") {
        df <- df[as.character(df[[group_column]]) %in% requested, , drop = FALSE]
      }
      df[[group_column]] <- droplevels(df[[group_column]])
      if (nlevels(df[[group_column]]) < 2L) {
        .tage_skip(vc, stratum$label, "fewer than two groups present")
        next
      }
      if (!reference_group %in% levels(df[[group_column]])) {
        .tage_skip(vc, stratum$label, sprintf("reference group '%s' absent", reference_group))
        next
      }

      emm <- .tage_fit_emm(df, vc, group_column, covariates, se_col)
      if (is.character(emm)) {
        .tage_skip(vc, stratum$label, emm)
        next
      }

      cont_df <- .tage_contrast_df(emm, group_column, reference_group, method,
                                   conf_level = conf_level)
      if (is.null(cont_df) || nrow(cont_df) == 0L) {
        .tage_skip(vc, stratum$label, "no contrast could be formed")
        next
      }

      # as.integer() strips the 1-d array structure that table() arithmetic
      # would otherwise carry into the result column.
      n_by_group <- table(as.character(df[[group_column]]))
      cont_df$n <- as.integer(n_by_group[cont_df$group1]) +
        as.integer(n_by_group[cont_df$group2])
      cont_df$value_column <- vc
      cont_df$split <- stratum$label
      rows[[length(rows) + 1L]] <- cont_df
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(value_column = character(0), split = character(0),
                      group1 = character(0), group2 = character(0), n = integer(0),
                      estimate = numeric(0), se = numeric(0), ci_low = numeric(0),
                      ci_high = numeric(0), statistic = numeric(0),
                      df = numeric(0), p_value = numeric(0), p_adjusted = numeric(0),
                      label = character(0), stringsAsFactors = FALSE))
  }

  res <- do.call(rbind, rows)

  # With variance_strata = "all_data" the model also carries groups that were
  # never requested; drop their contrasts before correcting p-values.
  if (method == "trt.vs.ctrl") {
    res <- res[res$group2 %in% compare_groups, , drop = FALSE]
  } else {
    res <- res[res$group1 %in% requested & res$group2 %in% requested, , drop = FALSE]
  }

  res <- res[, c("value_column", "split", "group1", "group2", "n", "estimate",
                 "se", "ci_low", "ci_high", "statistic", "df", "p_value")]
  .tage_finalise(res, p_adjust, p_adjust_scope)
}

# ---------------------------------------------------------------------------
# Continuous predictors
# ---------------------------------------------------------------------------

#' Regress tAge on a continuous predictor
#'
#' Slope of predicted tAge against a numeric predictor such as chronological
#' age, dose or time in culture. Elastic net clocks use \code{\link[stats]{lm}}; Bayesian ridge clocks
#' use a REML meta-regression weighted by the per-sample prediction standard
#' deviation and report z-tests.
#'
#' @param data Data frame of per-sample predictions.
#' @param value_columns Character vector of tAge columns.
#' @param predictor Numeric column regressed against.
#' @param covariates Additional fixed-effect columns. Default \code{NULL}.
#' @param split_by Column to stratify on. Default \code{NULL}.
#' @param se_columns Per-sample standard-deviation columns for Bayesian ridge
#'   clocks; see \code{\link{tage_compare_groups}}.
#' @param p_adjust Method passed to \code{\link[stats]{p.adjust}}.
#' @param p_adjust_scope Correction family; see \code{\link{tage_compare_groups}}.
#'   Because there is one test per clock and stratum, \code{"across_columns"}
#'   corrects across clocks within a stratum.
#'
#' @param conf_level Two-sided confidence level for \code{ci_low} / \code{ci_high}.
#'
#' @return A data frame with \code{value_column}, \code{split}, \code{term},
#'   \code{n}, \code{estimate} (slope per unit of \code{predictor}), \code{se},
#'   \code{ci_low}, \code{ci_high}, \code{statistic}, \code{df},
#'   \code{p_value}, \code{p_adjusted} and \code{label}.
#'
#' @examples
#' \dontrun{
#' tage_regress_continuous(
#'   results,
#'   value_columns = "scaled_diff_EN_tAge",
#'   predictor     = "age_months",
#'   split_by      = "Tissue"
#' )
#' }
#' @export
tage_regress_continuous <- function(data,
                                    value_columns,
                                    predictor,
                                    covariates = NULL,
                                    split_by = NULL,
                                    se_columns = NULL,
                                    p_adjust = TAGE_P_ADJUST_METHODS,
                                    p_adjust_scope = c("within_column", "across_columns",
                                                       "global", "none"),
                                    conf_level = 0.95) {
  p_adjust       <- match.arg(p_adjust, TAGE_P_ADJUST_METHODS)
  p_adjust_scope <- match.arg(p_adjust_scope)

  data          <- as.data.frame(data)
  value_columns <- .tage_clean_chr(value_columns)
  covariates    <- .tage_clean_chr(covariates)
  split_by      <- .tage_clean_chr(split_by)

  if (is.null(value_columns)) stop("`value_columns` must name at least one column.", call. = FALSE)
  missing <- setdiff(c(value_columns, predictor, covariates, split_by), colnames(data))
  if (length(missing)) stop("Column(s) not found in data: ", paste(missing, collapse = ", "), call. = FALSE)

  se_map <- .tage_resolve_se(se_columns, value_columns, data)

  # as.numeric() on a factor silently returns its integer level codes, so check
  # before coercing rather than fitting a slope through arbitrary codes.
  if (is.factor(data[[predictor]]) || is.character(data[[predictor]])) {
    stop("Predictor '", predictor, "' is not numeric. tage_regress_continuous() fits ",
         "a slope; use tage_compare_groups() for a categorical variable.", call. = FALSE)
  }
  data[[predictor]] <- as.numeric(data[[predictor]])

  rows <- list()

  for (i in seq_along(value_columns)) {
    vc <- value_columns[i]
    se_col <- se_map[i]

    keep_cols <- c(vc, predictor, covariates, split_by)
    if (!is.na(se_col)) keep_cols <- c(keep_cols, se_col)
    work <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]), , drop = FALSE]
    if (nrow(work) == 0L) {
      .tage_skip(vc, NA, "no complete cases (missing values in the value, predictor, covariate or split columns)")
      next
    }

    for (stratum in .tage_split(work, split_by)) {
      df <- stratum$data
      if (nrow(df) < 3L) {
        .tage_skip(vc, stratum$label, sprintf("only %d sample(s), need at least 3", nrow(df)))
        next
      }

      rhs <- .tage_rhs(c(predictor, covariates))

      if (is.na(se_col)) {
        fml <- stats::as.formula(paste(vc, "~", rhs))
        model <- try(stats::lm(fml, data = df), silent = TRUE)
        if (inherits(model, "try-error")) {
          .tage_skip(vc, stratum$label, paste("lm failed:", .tage_try_reason(model)))
          next
        }
        coefs <- stats::coef(summary(model))
        if (!predictor %in% rownames(coefs)) {
          .tage_skip(vc, stratum$label, sprintf("'%s' is not estimable (constant or collinear)", predictor))
          next
        }
        est <- coefs[predictor, "Estimate"]
        se  <- coefs[predictor, "Std. Error"]
        st  <- coefs[predictor, "t value"]
        pv  <- coefs[predictor, "Pr(>|t|)"]
        dfr <- stats::df.residual(model)
      } else {
        .tage_require("metafor")
        mods <- stats::as.formula(paste("~", rhs))
        model <- try(
          suppressWarnings(
            metafor::rma.uni(yi = df[[vc]], sei = df[[se_col]],
                             mods = mods, data = df, method = "REML")
          ),
          silent = TRUE
        )
        if (inherits(model, "try-error")) {
          .tage_skip(vc, stratum$label, paste("rma.uni failed:", .tage_try_reason(model)))
          next
        }
        nm <- rownames(model$beta)
        if (!predictor %in% nm) {
          .tage_skip(vc, stratum$label, sprintf("'%s' is not estimable (constant or collinear)", predictor))
          next
        }
        j   <- which(nm == predictor)
        est <- as.numeric(model$beta)[j]
        se  <- model$se[j]
        st  <- model$zval[j]
        pv  <- model$pval[j]
        dfr <- Inf
      }

      crit <- .tage_crit(dfr, conf_level)
      rows[[length(rows) + 1L]] <- data.frame(
        value_column = vc,
        split        = stratum$label,
        term         = predictor,
        n            = nrow(df),
        estimate     = est,
        se           = se,
        ci_low       = est - crit * se,
        ci_high      = est + crit * se,
        statistic    = st,
        df           = dfr,
        p_value      = pv,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(value_column = character(0), split = character(0),
                      term = character(0), n = integer(0), estimate = numeric(0),
                      se = numeric(0), ci_low = numeric(0), ci_high = numeric(0),
                      statistic = numeric(0), df = numeric(0),
                      p_value = numeric(0), p_adjusted = numeric(0),
                      label = character(0), stringsAsFactors = FALSE))
  }

  res <- do.call(rbind, rows)
  .tage_finalise(res, p_adjust, p_adjust_scope, cell_keys = c("split", "term"))
}

# ---------------------------------------------------------------------------
# Module clocks
# ---------------------------------------------------------------------------

#' Module-clock effect sizes and p-values
#'
#' Per-module statistics for the module-clock heatmaps. With
#' \code{standardize = TRUE} (the default) the two quantities come from
#' different models:
#' the p-value is the estimated marginal-mean contrast of the full model fitted
#' on the stratum, while the effect size is the coefficient of a separate
#' two-group model fitted on values standardised to unit variance, which makes
#' modules comparable despite their different native scales. With
#' \code{standardize = FALSE} the estimate, standard error and p-value all come
#' from the same contrast.
#'
#' @param data Data frame with one column per module clock.
#' @param module_columns Character vector of module-clock columns.
#' @param group_column Column defining the experimental groups.
#' @param reference_group Reference level of \code{group_column}.
#' @param compare_groups Levels compared against the reference. Default
#'   \code{NULL} uses every other level.
#' @param covariates Covariate columns. Default \code{NULL}.
#' @param split_by Column to stratify on. Default \code{NULL}.
#' @param variance_strata \code{"all_data"} (default, matching the
#'   application) fits the p-value model on every group in the stratum;
#'   \code{"subset"} restricts it to the compared groups.
#' @param standardize Whether to report standardised effect sizes. Default
#'   \code{TRUE}.
#' @param p_adjust Method passed to \code{\link[stats]{p.adjust}}.
#' @param p_adjust_scope Correction family. Default \code{"across_columns"}
#'   corrects across modules within each comparison and stratum, which is the
#'   application's "per group" setting; \code{"global"} matches its "globally"
#'   setting.
#'
#' @param conf_level Two-sided confidence level for \code{ci_low} / \code{ci_high}.
#'
#' @return A data frame with \code{module}, \code{split}, \code{group1},
#'   \code{group2}, \code{n}, \code{estimate}, \code{se}, \code{ci_low},
#'   \code{ci_high}, \code{statistic}, \code{df}, \code{p_value},
#'   \code{p_adjusted} and \code{label}.
#'
#' @examples
#' \dontrun{
#' tage_module_stats(
#'   results,
#'   module_columns  = grep("^module_", names(results), value = TRUE),
#'   group_column    = "Genotype",
#'   reference_group = "WT",
#'   split_by        = "Tissue"
#' )
#' }
#' @export
tage_module_stats <- function(data,
                              module_columns,
                              group_column,
                              reference_group,
                              compare_groups = NULL,
                              covariates = NULL,
                              split_by = NULL,
                              variance_strata = c("all_data", "subset"),
                              standardize = TRUE,
                              p_adjust = TAGE_P_ADJUST_METHODS,
                              p_adjust_scope = c("across_columns", "within_column",
                                                 "global", "none"),
                              conf_level = 0.95) {
  variance_strata <- match.arg(variance_strata)
  p_adjust        <- match.arg(p_adjust, TAGE_P_ADJUST_METHODS)
  p_adjust_scope  <- match.arg(p_adjust_scope)

  base <- tage_compare_groups(
    data            = data,
    value_columns   = module_columns,
    group_column    = group_column,
    reference_group = reference_group,
    compare_groups  = compare_groups,
    covariates      = covariates,
    split_by        = split_by,
    se_columns      = NULL,
    method          = "trt.vs.ctrl",
    variance_strata = variance_strata,
    p_adjust        = "none",
    p_adjust_scope  = "none",
    conf_level      = conf_level
  )
  names(base)[names(base) == "value_column"] <- "module"

  if (standardize && nrow(base) > 0L) {
    data       <- as.data.frame(data)
    covariates <- .tage_clean_chr(covariates)
    split_by   <- .tage_clean_chr(split_by)

    for (i in seq_len(nrow(base))) {
      mc <- base$module[i]
      g2 <- base$group2[i]

      keep_cols <- c(mc, group_column, covariates, split_by)
      df <- data[stats::complete.cases(data[, keep_cols, drop = FALSE]), , drop = FALSE]
      if (!is.null(split_by) && !is.na(base$split[i])) {
        df <- df[as.character(df[[split_by]]) == base$split[i], , drop = FALSE]
      }
      df <- df[as.character(df[[group_column]]) %in% c(reference_group, g2), , drop = FALSE]
      if (nrow(df) < 3L) {
        .tage_skip(mc, base$split[i], "standardised effect not computed: fewer than 3 samples in the two groups")
        next
      }

      df$.tage_group <- factor(as.character(df[[group_column]]),
                               levels = c(reference_group, g2))
      y <- df[[mc]]

      if (!is.null(covariates)) {
        cov_fml <- stats::as.formula(paste(mc, "~", .tage_rhs(covariates)))
        cov_model <- try(stats::lm(cov_fml, data = df), silent = TRUE)
        if (inherits(cov_model, "try-error")) {
          .tage_skip(mc, base$split[i], paste("standardised effect not computed:", .tage_try_reason(cov_model)))
          next
        }
        y <- stats::resid(cov_model)
      }

      sd_y <- stats::sd(y)
      if (!is.finite(sd_y) || sd_y <= 0) {
        .tage_skip(mc, base$split[i], "standardised effect not computed: values are constant")
        next
      }
      df$.tage_std <- (y - mean(y)) / sd_y

      slope_model <- try(stats::lm(.tage_std ~ .tage_group, data = df), silent = TRUE)
      if (inherits(slope_model, "try-error")) {
        .tage_skip(mc, base$split[i], paste("standardised effect not computed:", .tage_try_reason(slope_model)))
        next
      }
      coefs <- stats::coef(summary(slope_model))
      row <- paste0(".tage_group", g2)
      if (!row %in% rownames(coefs)) {
        .tage_skip(mc, base$split[i], "standardised effect not computed: group effect not estimable")
        next
      }

      base$estimate[i] <- coefs[row, "Estimate"]
      base$se[i]       <- coefs[row, "Std. Error"]
      # The interval has to follow the effect it describes, so rebuild it
      # from the standardised model rather than leaving the composite one.
      crit <- .tage_crit(base$df[i], conf_level)
      base$ci_low[i]  <- base$estimate[i] - crit * base$se[i]
      base$ci_high[i] <- base$estimate[i] + crit * base$se[i]
      # statistic and df keep describing the test that produced p_value.
    }
  }

  base <- base[, c("module", "split", "group1", "group2", "n", "estimate",
                   "se", "ci_low", "ci_high", "statistic", "df", "p_value")]
  .tage_finalise(base, p_adjust, p_adjust_scope, column_key = "module")
}

# ---------------------------------------------------------------------------
# Covariate-adjusted values for plotting
# ---------------------------------------------------------------------------

#' Covariate-adjusted tAge values for plotting
#'
#' Removes the fitted covariate effects from a tAge column while keeping the
#' group effect, so that box plots show what the model tested: the covariate
#' model is fitted without the grouping variable, and the mean covariate effect is added
#' back so the adjusted values stay on the original scale.
#'
#' @param data Data frame of per-sample predictions.
#' @param value_column Column to adjust.
#' @param covariates Covariate columns to regress out.
#' @param split_by Optional stratifying column, included as a fixed effect for
#'   elastic net clocks and used to fit one model per stratum for Bayesian
#'   ridge clocks.
#' @param se_column Per-sample standard deviations for a Bayesian ridge clock.
#'   Default \code{NULL} uses \code{\link[stats]{lm}}.
#'
#' @return Numeric vector of adjusted values, in the row order of \code{data};
#'   rows with missing values in the model return \code{NA}.
#'
#' @examples
#' \dontrun{
#' results$adjusted <- tage_adjust_covariates(
#'   results, "yugene_diff_EN_tAge", covariates = "Sex", split_by = "Tissue"
#' )
#' }
#' @export
tage_adjust_covariates <- function(data,
                                   value_column,
                                   covariates,
                                   split_by = NULL,
                                   se_column = NULL) {
  data       <- as.data.frame(data)
  covariates <- .tage_clean_chr(covariates)
  split_by   <- .tage_clean_chr(split_by)
  se_column  <- .tage_clean_chr(se_column)

  if (is.null(covariates)) return(as.numeric(data[[value_column]]))

  missing <- setdiff(c(value_column, covariates, split_by, se_column), colnames(data))
  if (length(missing)) stop("Column(s) not found in data: ", paste(missing, collapse = ", "), call. = FALSE)

  out <- rep(NA_real_, nrow(data))
  keep_cols <- c(value_column, covariates, split_by, se_column)
  ok <- stats::complete.cases(data[, keep_cols, drop = FALSE])
  if (!any(ok)) return(out)

  if (is.null(se_column)) {
    rhs <- .tage_rhs(c(split_by, covariates))
    fml <- stats::as.formula(paste(value_column, "~", rhs))
    model <- stats::lm(fml, data = data[ok, , drop = FALSE])

    mm <- stats::model.matrix(model)
    terms_all <- colnames(mm)
    cov_cols <- setdiff(terms_all, "(Intercept)")
    if (!is.null(split_by)) {
      cov_cols <- setdiff(cov_cols, grep(paste0("^", split_by), terms_all, value = TRUE))
    }
    if (length(cov_cols) == 0L) {
      out[ok] <- data[[value_column]][ok]
      return(out)
    }

    beta <- stats::coef(model)[cov_cols]
    beta[is.na(beta)] <- 0
    effects <- as.vector(mm[, cov_cols, drop = FALSE] %*% beta)
    mean_effect <- sum(colMeans(mm[, cov_cols, drop = FALSE], na.rm = TRUE) * beta, na.rm = TRUE)
    out[ok] <- data[[value_column]][ok] - effects + mean_effect
    return(out)
  }

  .tage_require("metafor")
  mods <- stats::as.formula(paste("~", .tage_rhs(covariates)))
  sub <- data[ok, , drop = FALSE]

  if (is.null(split_by)) {
    model <- metafor::rma.uni(yi = sub[[value_column]], sei = sub[[se_column]],
                              mods = mods, data = sub, method = "REML")
    # Residuals are centred on zero; add the mean back so the adjusted values
    # stay on the tAge scale, as the lm branch and the per-stratum branch do.
    out[ok] <- as.numeric(stats::resid(model)) + mean(sub[[value_column]], na.rm = TRUE)
    return(out)
  }

  idx <- which(ok)
  for (lev in unique(as.character(sub[[split_by]]))) {
    sel <- as.character(sub[[split_by]]) == lev
    part <- sub[sel, , drop = FALSE]
    model <- try(
      metafor::rma.uni(yi = part[[value_column]], sei = part[[se_column]],
                       mods = mods, data = part, method = "REML"),
      silent = TRUE
    )
    if (inherits(model, "try-error")) {
      .tage_skip(value_column, lev, paste("rma.uni failed:", .tage_try_reason(model)))
      next
    }
    out[idx[sel]] <- as.numeric(stats::resid(model)) + mean(part[[value_column]], na.rm = TRUE)
  }
  out
}
