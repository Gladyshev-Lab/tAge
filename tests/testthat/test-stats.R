# Statistics tests. The reference values are computed inline with the same
# engine the TACO / tClock application uses (lm + emmeans for elastic net,
# metafor::rma.uni + qdrg for Bayesian ridge), so these lock the package to the
# application rather than to hand-copied numbers.

.stats_fixture <- function(seed = 42, n = 90) {
  set.seed(seed)
  d <- data.frame(
    Genotype = factor(rep(c("WT", "KO", "HET"), each = n / 3),
                      levels = c("WT", "KO", "HET")),
    Tissue   = factor(rep(c("Kidney", "Muscle"), times = n / 2)),
    Sex      = factor(sample(c("M", "F"), n, replace = TRUE)),
    Age      = runif(n, 3, 24)
  )
  d$tAge <- 0.4 * (d$Genotype == "KO") + 0.15 * (d$Genotype == "HET") +
    0.3 * (d$Sex == "M") + 0.02 * d$Age + rnorm(n, sd = 0.5)
  d$tAge_sd <- runif(n, 0.1, 0.6)
  d
}

# Reference implementation: the application's contrast code, transcribed.
.app_contrasts <- function(df, group1, Group_col, covariates = NULL,
                           split_group = NULL, use_br = FALSE,
                           variance_strata = "subset") {
  group2 <- setdiff(levels(df[[Group_col]]), group1)
  all_groups <- c(group1, group2)
  df$Temp_tAge <- df$tAge
  df$Temp_se <- df$tAge_sd

  split_list <- if (is.null(split_group)) {
    list(full_data = df)
  } else {
    split(df, factor(df[[split_group]]))
  }

  out <- data.frame()
  for (split_name in names(split_list)) {
    split_df <- split_list[[split_name]]
    if (variance_strata == "subset") {
      split_df <- subset(split_df, split_df[[Group_col]] %in% all_groups)
    }
    split_df <- droplevels(split_df)
    rhs <- paste(c(Group_col, covariates), collapse = " + ")

    if (!use_br) {
      model <- stats::lm(stats::as.formula(paste("Temp_tAge ~", rhs)), data = split_df)
      emm <- emmeans::emmeans(model, specs = Group_col)
    } else {
      fml <- stats::as.formula(paste("~", rhs))
      m <- metafor::rma.uni(yi = Temp_tAge, sei = Temp_se, mods = fml, data = split_df)
      qrg <- emmeans::qdrg(formula = fml, data = split_df, coef = stats::coef(m),
                           vcov = stats::vcov(m), df = Inf)
      emm <- emmeans::emmeans(qrg, specs = Group_col)
    }

    lv <- as.character(emm@grid[[Group_col]])
    ref_position <- which(lv == group1)
    cont <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_position,
                              adjust = "none")
    cont_df <- as.data.frame(cont)
    parsed <- do.call(rbind, lapply(cont_df$contrast, function(lbl) {
      parts <- strsplit(as.character(lbl), " - ", fixed = TRUE)[[1]]
      data.frame(lhs = sub(paste0("^", Group_col), "", trimws(parts[1])),
                 rhs = sub(paste0("^", Group_col), "", trimws(parts[2])),
                 stringsAsFactors = FALSE)
    }))
    stat_col <- if ("t.ratio" %in% names(cont_df)) cont_df$t.ratio else cont_df$z.ratio
    out <- rbind(out, data.frame(
      group2    = parsed$lhs,
      estimate  = cont_df$estimate,
      se        = cont_df$SE,
      statistic = stat_col,
      p_value   = cont_df$p.value,
      split     = split_name,
      stringsAsFactors = FALSE
    ))
  }
  out[order(out$split, out$group2), ]
}

.sorted <- function(res) {
  key <- ifelse(is.na(res$split), "full_data", as.character(res$split))
  res[order(key, res$group2), ]
}

test_that("significance stars follow the reference thresholds", {
  expect_equal(tage_significance_stars(c(1e-4, 0.005, 0.03, 0.08, 0.5)),
               c("***", "**", "*", "^", ""))
  expect_true(is.na(tage_significance_stars(NA_real_)))
  # Thresholds are strict "<", as in the application: exactly 0.05 is only a
  # trend, and exactly 0.1 is nothing.
  expect_equal(tage_significance_stars(0.05), "^")
  expect_equal(tage_significance_stars(0.1), "")
})

test_that("elastic net contrasts match the application", {
  d <- .stats_fixture()

  res <- .sorted(tage_compare_groups(d, "tAge", "Genotype", "WT", p_adjust = "none"))
  ref <- .app_contrasts(d, "WT", "Genotype")
  expect_equal(res$estimate, ref$estimate)
  expect_equal(res$se, ref$se)
  expect_equal(res$statistic, ref$statistic)
  expect_equal(res$p_value, ref$p_value)
  expect_equal(as.character(res$group1), rep("WT", nrow(res)))
})

test_that("covariates and stratification match the application", {
  d <- .stats_fixture()

  res <- .sorted(tage_compare_groups(d, "tAge", "Genotype", "WT",
                                     covariates = c("Sex", "Age"),
                                     split_by = "Tissue", p_adjust = "none"))
  ref <- .app_contrasts(d, "WT", "Genotype", covariates = c("Sex", "Age"),
                        split_group = "Tissue")
  expect_equal(res$estimate, ref$estimate)
  expect_equal(res$p_value, ref$p_value)
  expect_setequal(unique(res$split), c("Kidney", "Muscle"))
})

test_that("Bayesian ridge clocks go through the weighted meta-regression", {
  d <- .stats_fixture()

  res <- .sorted(tage_compare_groups(d, "tAge", "Genotype", "WT",
                                     se_columns = "tAge_sd", p_adjust = "none"))
  ref <- .app_contrasts(d, "WT", "Genotype", use_br = TRUE)
  expect_equal(res$estimate, ref$estimate)
  expect_equal(res$se, ref$se)
  expect_equal(res$p_value, ref$p_value)

  # Meta-regression reports z-tests, so the degrees of freedom are infinite.
  expect_true(all(is.infinite(res$df)))

  # And it must differ from the unweighted fit; otherwise the weights are
  # silently being ignored.
  unweighted <- .sorted(tage_compare_groups(d, "tAge", "Genotype", "WT", p_adjust = "none"))
  expect_false(isTRUE(all.equal(res$estimate, unweighted$estimate)))
})

test_that("variance_strata controls which groups enter the model", {
  d <- .stats_fixture()

  subset_fit <- tage_compare_groups(d, "tAge", "Genotype", "WT",
                                    compare_groups = "KO",
                                    variance_strata = "subset", p_adjust = "none")
  all_fit <- tage_compare_groups(d, "tAge", "Genotype", "WT",
                                 compare_groups = "KO",
                                 variance_strata = "all_data", p_adjust = "none")

  expect_equal(nrow(subset_fit), 1L)
  expect_equal(nrow(all_fit), 1L)
  # Same point estimate, different residual variance and hence different SE.
  expect_equal(subset_fit$estimate, all_fit$estimate)
  expect_false(isTRUE(all.equal(subset_fit$se, all_fit$se)))
})

test_that("estimate is always group2 minus group1", {
  d <- .stats_fixture()
  res <- tage_compare_groups(d, "tAge", "Genotype", "WT", compare_groups = "KO",
                             p_adjust = "none")
  raw_diff <- mean(d$tAge[d$Genotype == "KO"]) - mean(d$tAge[d$Genotype == "WT"])
  expect_equal(res$estimate, raw_diff)

  # Same for the Bayesian ridge path, which the application reports with the
  # opposite sign to its own elastic net output.
  res_br <- tage_compare_groups(d, "tAge", "Genotype", "WT", compare_groups = "KO",
                                se_columns = "tAge_sd", p_adjust = "none")
  expect_gt(res_br$estimate, 0)
})

test_that("p-value adjustment scopes define the right families", {
  d <- .stats_fixture()
  d$tAge2 <- d$tAge + rnorm(nrow(d), sd = 0.4)

  raw <- tage_compare_groups(d, c("tAge", "tAge2"), "Genotype", "WT",
                             split_by = "Tissue", p_adjust = "none")
  expect_equal(raw$p_adjusted, raw$p_value)

  within <- tage_compare_groups(d, c("tAge", "tAge2"), "Genotype", "WT",
                                split_by = "Tissue", p_adjust = "BH",
                                p_adjust_scope = "within_column")
  for (vc in unique(within$value_column)) {
    idx <- within$value_column == vc
    expect_equal(within$p_adjusted[idx],
                 stats::p.adjust(within$p_value[idx], "BH"))
  }

  across <- tage_compare_groups(d, c("tAge", "tAge2"), "Genotype", "WT",
                                split_by = "Tissue", p_adjust = "BH",
                                p_adjust_scope = "across_columns")
  key <- paste(across$split, across$group2)
  for (k in unique(key)) {
    idx <- key == k
    expect_equal(across$p_adjusted[idx], stats::p.adjust(across$p_value[idx], "BH"))
    expect_equal(sum(idx), 2L)  # one per clock
  }

  global <- tage_compare_groups(d, c("tAge", "tAge2"), "Genotype", "WT",
                                split_by = "Tissue", p_adjust = "BH",
                                p_adjust_scope = "global")
  expect_equal(global$p_adjusted, stats::p.adjust(global$p_value, "BH"))
})

test_that("pairwise mode returns every pair once, in the documented direction", {
  d <- .stats_fixture()
  res <- tage_compare_groups(d, "tAge", "Genotype", "WT", method = "pairwise",
                             p_adjust = "none")
  expect_equal(nrow(res), 3L)
  pairs <- paste(res$group1, res$group2, sep = "-")
  expect_setequal(pairs, c("WT-KO", "WT-HET", "KO-HET"))

  # estimate is group2 - group1 for every method, not just trt.vs.ctrl.
  for (i in seq_len(nrow(res))) {
    expected <- mean(d$tAge[d$Genotype == res$group2[i]]) -
      mean(d$tAge[d$Genotype == res$group1[i]])
    expect_equal(res$estimate[i], expected)
  }

  # And the trt.vs.ctrl rows agree with the matching pairwise rows.
  trt <- tage_compare_groups(d, "tAge", "Genotype", "WT", p_adjust = "none")
  shared <- merge(trt, res, by = c("group1", "group2"), suffixes = c(".trt", ".pw"))
  expect_equal(nrow(shared), 2L)
  expect_equal(shared$estimate.trt, shared$estimate.pw)
})

test_that("continuous regression returns the predictor coefficient", {
  d <- .stats_fixture()

  res <- tage_regress_continuous(d, "tAge", "Age", covariates = "Sex",
                                 p_adjust = "none")
  ref <- stats::coef(summary(stats::lm(tAge ~ Age + Sex, data = d)))["Age", ]
  expect_equal(res$estimate, unname(ref["Estimate"]))
  expect_equal(res$se, unname(ref["Std. Error"]))
  expect_equal(res$statistic, unname(ref["t value"]))
  expect_equal(res$p_value, unname(ref["Pr(>|t|)"]))
  expect_equal(res$n, nrow(d))

  res_br <- tage_regress_continuous(d, "tAge", "Age", covariates = "Sex",
                                    se_columns = "tAge_sd", p_adjust = "none")
  m <- metafor::rma.uni(yi = tAge, sei = tAge_sd, mods = ~ Age + Sex, data = d)
  j <- which(rownames(m$beta) == "Age")
  expect_equal(res_br$estimate, as.numeric(m$beta)[j])
  expect_equal(res_br$p_value, m$pval[j])
})

test_that("module statistics separate the effect size from the test", {
  d <- .stats_fixture()
  d$mod_a <- d$tAge
  d$mod_b <- d$tAge * 10 + 5   # same signal, different scale

  res <- tage_module_stats(d, c("mod_a", "mod_b"), "Genotype", "WT",
                           compare_groups = "KO", p_adjust = "none")
  expect_equal(nrow(res), 2L)

  # Standardised effect sizes make differently scaled modules comparable.
  expect_equal(res$estimate[res$module == "mod_a"],
               res$estimate[res$module == "mod_b"])

  # The p-value comes from the unstandardised model, so it is scale-free too.
  expect_equal(res$p_value[res$module == "mod_a"],
               res$p_value[res$module == "mod_b"])

  # Without standardisation the estimates track the native scale instead.
  raw <- tage_module_stats(d, c("mod_a", "mod_b"), "Genotype", "WT",
                           compare_groups = "KO", standardize = FALSE,
                           p_adjust = "none")
  expect_equal(raw$estimate[raw$module == "mod_b"],
               10 * raw$estimate[raw$module == "mod_a"])
})

test_that("module statistics default to correcting across modules", {
  d <- .stats_fixture()
  d$mod_a <- d$tAge
  d$mod_b <- d$tAge + rnorm(nrow(d), sd = 0.5)
  d$mod_c <- rnorm(nrow(d))

  res <- tage_module_stats(d, c("mod_a", "mod_b", "mod_c"), "Genotype", "WT",
                           compare_groups = "KO", p_adjust = "BH")
  expect_equal(res$p_adjusted, stats::p.adjust(res$p_value, "BH"))
})

test_that("covariate adjustment keeps the group effect and the scale", {
  d <- .stats_fixture()
  adj <- tage_adjust_covariates(d, "tAge", covariates = "Age")

  expect_length(adj, nrow(d))
  # Mean is preserved: the mean covariate effect is added back.
  expect_equal(mean(adj), mean(d$tAge))
  # The covariate is gone from the adjusted values.
  expect_lt(abs(stats::cor(adj, d$Age)), 1e-8)
  # The group signal survives.
  expect_gt(mean(adj[d$Genotype == "KO"]) - mean(adj[d$Genotype == "WT"]), 0)

  # No covariates is a no-op rather than an error.
  expect_equal(tage_adjust_covariates(d, "tAge", covariates = NULL), d$tAge)
})

test_that("missing values and bad arguments are handled", {
  d <- .stats_fixture()
  d$Age[1:5] <- NA

  # Rows missing a covariate are dropped, matching complete.cases in the app.
  # n counts the two groups entering each contrast, not the whole table.
  res <- tage_compare_groups(d, "tAge", "Genotype", "WT", covariates = "Age",
                             p_adjust = "none")
  complete <- d[!is.na(d$Age), ]
  expected_n <- vapply(res$group2, function(g) {
    sum(complete$Genotype == "WT") + sum(complete$Genotype == g)
  }, numeric(1))
  expect_equal(res$n, unname(expected_n))

  expect_error(tage_compare_groups(d, "nope", "Genotype", "WT"),
               "not found in data")
  expect_error(tage_compare_groups(d, "tAge", "Genotype", "Nonexistent"),
               "not present")
  expect_error(tage_compare_groups(d, "tAge", "Genotype", "WT",
                                   compare_groups = "Nope"),
               "not present")
})

test_that("a factor level starting with the column name is parsed correctly", {
  # The application derives group labels by stripping the column name off the
  # emmeans contrast string, which mangles levels like "GenotypeKO".
  d <- .stats_fixture()
  levels(d$Genotype) <- c("WT", "GenotypeKO", "HET")

  res <- tage_compare_groups(d, "tAge", "Genotype", "WT",
                             compare_groups = "GenotypeKO", p_adjust = "none")
  expect_equal(res$group2, "GenotypeKO")
  expect_equal(res$estimate,
               mean(d$tAge[d$Genotype == "GenotypeKO"]) - mean(d$tAge[d$Genotype == "WT"]))
})

test_that("a categorical predictor is rejected rather than coerced to codes", {
  d <- .stats_fixture()
  # as.numeric() on a factor yields level codes, which would fit a meaningless
  # slope through 1, 2, 3.
  expect_error(tage_regress_continuous(d, "tAge", "Genotype"), "not numeric")
  expect_error(
    tage_regress_continuous(d, "tAge", "Sex"), "not numeric"
  )
})

test_that("strata that cannot be tested are reported, not dropped in silence", {
  d <- .stats_fixture()
  # No WT in Muscle: that stratum has no reference group.
  d <- d[!(d$Tissue == "Muscle" & d$Genotype == "WT"), ]

  expect_warning(
    res <- tage_compare_groups(d, "tAge", "Genotype", "WT", split_by = "Tissue"),
    "Skipping tAge \\[Muscle\\]: reference group 'WT' absent"
  )
  expect_setequal(unique(res$split), "Kidney")
})

test_that("a collinear covariate is reported with the stratum it broke", {
  d <- .stats_fixture()
  d$Age2 <- d$Age * 2
  expect_warning(
    res <- tage_compare_groups(d, "tAge", "Genotype", "WT", covariates = c("Age", "Age2")),
    "Skipping tAge: redundant predictor"
  )
  expect_equal(nrow(res), 0L)
})

test_that("regress_continuous reports strata that are too small", {
  d <- .stats_fixture()
  d <- d[c(which(d$Tissue == "Kidney"), which(d$Tissue == "Muscle")[1:2]), ]
  expect_warning(
    res <- tage_regress_continuous(d, "tAge", "Age", split_by = "Tissue"),
    "Skipping tAge \\[Muscle\\]: only 2 sample"
  )
  expect_setequal(unique(res$split), "Kidney")
})

test_that("Bayesian ridge covariate adjustment keeps the tAge scale without split_by", {
  skip_if_not_installed("metafor")
  d <- .stats_fixture()
  adj_all   <- tage_adjust_covariates(d, "tAge", covariates = "Age", se_column = "tAge_sd")
  adj_split <- tage_adjust_covariates(d, "tAge", covariates = "Age", se_column = "tAge_sd",
                                      split_by = "Tissue")
  adj_lm    <- tage_adjust_covariates(d, "tAge", covariates = "Age")
  # Residuals alone are centred on zero; the adjusted values must sit where the
  # data sit, as the lm and per-stratum branches already did.
  # rma.uni residuals are weighted, so they only average out to ~zero.
  expect_equal(mean(adj_all), mean(d$tAge), tolerance = 1e-3)
  expect_equal(mean(adj_split), mean(d$tAge), tolerance = 0.05)
  expect_equal(mean(adj_lm), mean(d$tAge), tolerance = 1e-8)
})
