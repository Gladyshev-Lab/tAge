Thanks — both points were fair, and both are addressed in **tAge 1.4.0**
(`devtools::install_github("Gladyshev-Lab/tAge")`).

**`species`** is the species of your *samples*, not of the model. Its only effect is the
rescaling of chronological-age clocks from normalised age to age units (mouse × 4 y,
rat × 3.8 y, macaque × 39 y, human × 122 y — reported in months for rodents, years for
primates); mortality and normalized-age clocks ignore it. "Mouse / Rodents / Multispecies"
in `list_clocks()` describes the training set and does not enter this step. You were right
that it duplicated the argument of `tAge_preprocessing()`: since 1.4.0 the preprocessed
object records the species, so

```r
predict_tAge(tAge_eset, model_paths, mode = "EN")
```

works without it (an explicit `species` still overrides). The documentation now says all
this explicitly; see also `tage_species()`.

**BR standard deviations** were indeed discarded — a bug. `predict_tAge(..., mode = "BR")`
now returns them by default as `<normalisation>_BR_tAge_sd` (rescaled together with the
prediction for chronological clocks). For the downstream mixed-effects meta-regression:

```r
tage_compare_groups(
  results,
  value_columns   = "yugene_diff_BR_tAge",
  se_columns      = "yugene_diff_BR_tAge_sd",
  group_column    = "Genotype",
  reference_group = "WT",
  covariates      = "Sex",        # optional
  split_by        = "Tissue"      # optional
)
```

fits `metafor::rma.uni(yi, sei = sd, mods = ~ group + covariates, method = "REML")` with
contrasts via `emmeans::qdrg` (z-tests), the same model used for BR clocks in the paper;
`tage_regress_continuous()` does the same for a continuous predictor. Note the value passed
as `sei` is the model's predictive SD.

Two related changes in 1.4.0 you may want: `tAge_preprocessing(split_by = "Tissue", ...)`
preprocesses and centres each tissue on its own controls (how the clocks were trained), and
the species maximum lifespans were corrected to the training values (rat 3.8 y, human
122 y), so rat chronological predictions are ~10 % lower than in 1.1.0.

Closing as fixed in 1.4.0 — please reopen if anything doesn't match.
