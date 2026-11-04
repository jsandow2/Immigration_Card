# =============================================================================
# 21_wage_iv.R  --  Section 13: Wage IV models by occupation group
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 14_wage_iv.R.
# Renumbered 2026-10-21 from 13_wage_iv.R (archive-and-renumber pass).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

# CHECKPOINT LOAD REMOVED 2026-09-30, for the same reasons as in
# 20_migration_iv.R. This loaded after_sec12.RData and then overrode its
# contents with a fresh read of analysis_data.RData below; the remaining
# objects it carried (bpld_lookup, top175_list_named, meta_summary_full,
# multi_boot_summary) are never referenced here. With Sections 10-11 retired
# the pooled objects no longer exist on a clean rebuild.
#
# The specification here is unaffected by the 2026-09-30 respecification of the
# education design. This script uses total_actual_inflow and
# total_predicted_inflow at the SURVEY year, which is correct for its question:
# the effect of contemporaneous immigrant inflows on contemporaneous wages. The
# timing fix in 07_exposure_iv.R applied to the ages-15-25 exposure design,
# where instrument and regressor described different periods. Here they do not,
# so avg_exposure, avg_actual_exposure and exposure_missing do not enter.

# Section 13: Wage IV models by occupation group
#
# IV estimates of the effect of total_actual_inflow on native hourly wages,
# disaggregated by occupation group and demographic group. Two shift-share
# instruments are used:
#   (a) MSA x YEAR level: total_predicted_inflow (from Section 7)
#   (b) MSA x Occup_Group x YEAR level: occupation-specific shift-share
#
# The occupation-level shares follow the same national share definition as
# the main IV: s_{m,k,occ,1990} = fraction of all U.S. immigrants from
# origin k in occupation occ in 1990 who were located in MSA m.
# -----------------------------------------------------------------------------

# reduced_data_iv is not carried through the checkpoint chain past Section 8
# (see 08_load_analysis_subsets.R) -- reload it fresh here directly from
# analysis_data.RData, via a throwaway environment so nothing else already
# in memory (e.g. a processed all_natives_175) gets overwritten.
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
reduced_data_iv      <- .tmp_env$reduced_data_iv
natives_unrestricted <- .tmp_env$all_natives_175
rm(.tmp_env)

# -----------------------------------------------------------------------------
# Wage sample (added 2026-09-21)
#
# SAMPLE: ages 16-68, built from the unrestricted all_natives_175 rather than
# the checkpoint copy. The checkpoint carries ESTIMATION_MIN_AGE <- 30L from
# 08_load_analysis_subsets.R, a floor justified by an argument specific to
# EDUCATION outcomes (the ages-15-25 exposure window must have closed and
# attainment settled before the outcome is read). The outcome here is an
# hourly wage, so that rationale does not transfer. 16-68 matches Card (2001)
# and Section 12's Spec 1. Usable wage observations roughly double, from
# 699,619 under the old floor to 1,130,823.
#
# The 68 cap is not cosmetic: the 69+ band is badly behaved as a wage sample
# (p99 of $741/hr, 9.4% above $90/hr, 1st percentile of exactly $0.00),
# reflecting partial retirement, consulting and board income rather than a
# wage-setting labour market.
#
# CLEANING: hrwage is built in 03_variable_construction.R as
# INCEARN / (UHRSWORK * WEAKWK), set NA only when hours or weeks are zero or
# missing -- there is no trimming at either tail and no positivity check.
# Unfiltered, the sample contains negative hourly wages (self-employment
# losses, meaningless as a wage) and a maximum of $52,615/hr. Because
# mean_hrwage below is a cell mean over MSA x occupation x year, a handful of
# such values can move an entire cell.
#
# Card sets wages below $2/hr or above $90/hr to missing, but those are 1990
# dollars against a panel running to 2022: a fixed nominal ceiling would trim
# far more aggressively in later years purely from wage growth. Trimming to
# the 1st and 99th percentile WITHIN YEAR is the scale-free analogue and
# needs no external deflator. Applied once here, before the sex splits, so
# all three subsamples share identical thresholds.
# -----------------------------------------------------------------------------
WAGE_AGE_MIN <- 16L
WAGE_AGE_MAX <- 68L

# Survey year 2023 (added 2026-09-22).
#
# 06_shift_share_iv.R drops arrival years observable by fewer than two waves,
# which removes YRIMMIG == 2023 from the inflow series entirely. The merge
# onto individuals zero-fills rather than leaving NA, so every 2023 respondent
# carries actual_inflow == 0 -- structurally absent, not measured as zero.
#
# In the education and migration models a YEAR fixed effect absorbs that
# block. These wage regressions have NO fixed effects
# (mean_hrwage ~ actual_inflow + mean_age + mean_yrs_sch + frac_black), so
# 2023 cells would enter as genuine observations sitting at x = 0 with real
# wages attached, exerting leverage on the slope rather than being absorbed.
# Roughly one cell in nineteen. They are therefore excluded here.
#
# RETAINED (FALSE), decided 2026-09-22. Excluding 2023 raises every estimate
# -- Low-Skill Labor & Service goes from 0.011 (n.s.) to 0.062 per SD at 1%,
# and six of seven groups turn significant. Keeping the zero-inflow cells
# therefore attenuates the wage coefficients toward zero, which is the
# conservative direction and appropriate for a result this chapter reports as
# a rough first pass rather than a central claim. Set TRUE to exclude them and
# obtain the larger estimates.
#
# Dropped BEFORE the wage cleaning below, when enabled, so the within-year
# percentile cuts and retention figures describe the sample estimated on.
WAGE_DROP_2023 <- FALSE

natives_wage <- natives_unrestricted %>%
  filter(AGE >= WAGE_AGE_MIN, AGE <= WAGE_AGE_MAX)

if (WAGE_DROP_2023) {
  n_pre_2023 <- nrow(natives_wage)
  natives_wage <- natives_wage %>%
    filter(as.integer(as.character(YEAR)) != 2023L)
  cat(sprintf("Dropped survey year 2023 (zero inflow by construction): %d -> %d obs\n",
              n_pre_2023, nrow(natives_wage)))
}

n_pre_wage_clean <- sum(!is.na(natives_wage$hrwage))

natives_wage <- natives_wage %>%
  mutate(hrwage = ifelse(!is.na(hrwage) & hrwage <= 0, NA_real_, hrwage)) %>%
  group_by(YEAR) %>%
  mutate(
    lo_cut = quantile(hrwage, 0.01, na.rm = TRUE),
    hi_cut = quantile(hrwage, 0.99, na.rm = TRUE),
    hrwage = ifelse(!is.na(hrwage) & (hrwage < lo_cut | hrwage > hi_cut),
                    NA_real_, hrwage)
  ) %>%
  ungroup() %>%
  select(-lo_cut, -hi_cut)

cat(sprintf("Wage sample (AGE %d-%d): %d obs, %d with hrwage before cleaning\n",
            WAGE_AGE_MIN, WAGE_AGE_MAX, nrow(natives_wage), n_pre_wage_clean))
cat(sprintf("  After dropping non-positive and trimming 1st/99th pct by year: %d (%.1f%% retained)\n",
            sum(!is.na(natives_wage$hrwage)),
            100 * sum(!is.na(natives_wage$hrwage)) / n_pre_wage_clean))

rm(natives_unrestricted)

# Occupation-level actual inflows -- aggregated to MSA x Occup_Group x YEAR
# level directly from immigrant records. Do not join back to all_natives_175
# as that would re-expand to individual level.
actual_inflows_occ_year <- reduced_data_iv %>%
  filter(immigr == 1, YRIMMIG != 0, !is.na(Occup_Group)) %>%
  group_by(MSA        = as.character(MSA),
           Occup_Group,
           YEAR       = as.character(YRIMMIG)) %>%
  summarize(actual_inflow = sum(PERWT, na.rm = TRUE), .groups = "drop")

# Occupation-level shift-share IV
immigrant_shares_occ <- reduced_data_iv %>%
  filter(immigr == 1, YEAR == 1990, !is.na(Occup_Group)) %>%
  group_by(MSA, BPLD, Occup_Group) %>%
  summarize(weighted_count = sum(PERWT, na.rm = TRUE), .groups = "drop") %>%
  group_by(BPLD, Occup_Group) %>%
  mutate(
    national_total = sum(weighted_count),
    share          = weighted_count / national_total
  ) %>%
  ungroup() %>%
  select(MSA, BPLD, Occup_Group, share)

national_inflows_occ <- reduced_data_iv %>%
  filter(immigr == 1, YRIMMIG != 0, !is.na(Occup_Group)) %>%
  group_by(YEAR = YRIMMIG, BPLD, Occup_Group) %>%
  summarize(national_inflow = sum(PERWT, na.rm = TRUE), .groups = "drop")

# Full MSA x BPLD x Occup_Group x YEAR grid
msa_list    <- distinct(reduced_data_iv %>% filter(YEAR == 1990), MSA)
bpld_list   <- distinct(immigrant_shares_occ, BPLD)
occup_list  <- distinct(immigrant_shares_occ, Occup_Group)
year_list   <- distinct(national_inflows_occ, YEAR)

msa_bpld_occ_year_grid <- crossing(
  MSA         = msa_list$MSA,
  BPLD        = bpld_list$BPLD,
  Occup_Group = occup_list$Occup_Group,
  YEAR        = year_list$YEAR
)

shift_share_occ_data <- msa_bpld_occ_year_grid %>%
  left_join(immigrant_shares_occ,
            by = c("MSA", "BPLD", "Occup_Group")) %>%
  mutate(share = replace_na(share, 0)) %>%
  left_join(national_inflows_occ,
            by = c("YEAR", "BPLD", "Occup_Group")) %>%
  mutate(
    national_inflow  = replace_na(national_inflow, 0),
    predicted_inflow = share * national_inflow
  )

predicted_inflows_occ <- shift_share_occ_data %>%
  group_by(MSA, Occup_Group, YEAR) %>%
  summarize(
    predicted_inflow = sum(predicted_inflow, na.rm = TRUE),
    .groups = "drop"
  )

# Summary controls: weighted mean wage, age, schooling, fraction Black
build_wage_controls <- function(data, sex_filter = NULL) {
  if (!is.null(sex_filter)) data <- data %>% filter(SEX == sex_filter)
  data %>%
    filter(!is.na(hrwage), !is.na(Occup_Group)) %>%
    group_by(MSA  = as.character(MSA),
             Occup_Group,
             YEAR = as.character(YEAR)) %>%
    summarize(
      mean_hrwage  = weighted.mean(hrwage,            PERWT, na.rm = TRUE),
      mean_age     = weighted.mean(AGE,               PERWT, na.rm = TRUE),
      mean_yrs_sch = weighted.mean(YearSch,           PERWT, na.rm = TRUE),
      frac_black   = weighted.mean(as.numeric(Black), PERWT, na.rm = TRUE),
      # Respondents behind this cell mean. Deliberately unweighted: sum(PERWT)
      # would report implied population, which says nothing about precision --
      # one respondent carrying a weight of 150 implies 150 people but is
      # still a single observation. n_obs is what governs how noisy
      # mean_hrwage is. (The exactly-right measure under weighting is Kish's
      # effective N, (sum w)^2 / sum(w^2); n_obs is a close enough proxy for
      # screening thin cells.)
      n_obs        = n(),
      .groups = "drop"
    ) %>%
    left_join(
      natives_wage %>%
        distinct(MSA, YEAR, total_actual_inflow, total_predicted_inflow) %>%
        mutate(MSA = as.character(MSA), YEAR = as.character(YEAR)),
      by = c("MSA", "YEAR")
    ) %>%
    left_join(
      predicted_inflows_occ %>%
        mutate(MSA = as.character(MSA), YEAR = as.character(YEAR)),
      by = c("MSA", "Occup_Group", "YEAR")
    ) %>%
    left_join(actual_inflows_occ_year %>%
                select(MSA, Occup_Group, YEAR, actual_inflow),
              by = c("MSA", "Occup_Group", "YEAR")) %>%
    mutate(
      actual_inflow    = replace_na(actual_inflow, 0),
      predicted_inflow = replace_na(predicted_inflow, 0)
    )
}

iv_wage_data_all   <- build_wage_controls(natives_wage)
iv_wage_data_men   <- build_wage_controls(natives_wage, sex_filter = 1)
iv_wage_data_women <- build_wage_controls(natives_wage, sex_filter = 2)

# ELIGIBLE-WAVES DIVISOR (added 2026-11-02; approved specification change).
#
# The occupation-level inflow above is a raw PERWT sum by (MSA, Occup_Group,
# arrival year) pooled across every survey wave, with no division by the number
# of waves that could have observed that arrival year. An arrival year early in
# the panel is counted by many waves and 2023 by one, so measured inflow falls
# mechanically across the panel while nominal wages rise, putting a negative
# trend into the between-year covariance of wage and inflow.
#
# 06_shift_share_iv.R divides by n_waves for precisely this reason and drops
# arrival years with fewer than two eligible waves rather than dividing them.
# The same treatment is applied here, which also brings the wage exercise onto
# the same inflow construction as the rest of the chapter. Measured effect: the
# pooled and year-fixed-effect estimates, which differed by factors of three to
# twenty-five on the undivided series, come within about 20 per cent of each
# other once divided -- year effects had been doing the divisor's work.
#
# The n_waves >= 2 rule drops 2023, so WAGE_DROP_2023 above is now redundant
# for the inflow series; it is left in place because it also governs the wage
# sample itself.
.waves <- sort(unique(as.integer(reduced_data_iv$YEAR)))
.eligible <- data.frame(
  YEAR    = as.character(.waves),
  n_waves = vapply(.waves, function(y) sum(.waves >= y), integer(1)),
  stringsAsFactors = FALSE
)
.eligible <- .eligible[.eligible$n_waves >= 2L, ]

apply_wave_divisor <- function(d, label) {
  n0 <- nrow(d)
  d <- d %>%
    inner_join(.eligible, by = "YEAR") %>%
    mutate(actual_inflow    = actual_inflow    / n_waves,
           predicted_inflow = predicted_inflow / n_waves)
  cat(sprintf("  %-5s %s cells kept of %s, %d arrival years\n",
              paste0(label, ":"), format(nrow(d), big.mark = ","),
              format(n0, big.mark = ","), length(unique(d$YEAR))))
  d
}

cat("\nEligible-waves divisor (matching 06_shift_share_iv.R):\n")
iv_wage_data_all   <- apply_wave_divisor(iv_wage_data_all,   "All")
iv_wage_data_men   <- apply_wave_divisor(iv_wage_data_men,   "Men")
iv_wage_data_women <- apply_wave_divisor(iv_wage_data_women, "Women")

# Minimum cell count per occupation group, applied uniformly (2026-09-21).
# Previously this ran on the women's subsample only, which looked like a
# patch for a convergence failure rather than a design choice and left the
# three subsamples estimated on different inclusion rules.
MIN_CELLS_PER_OCCUP <- 30L

drop_thin_groups <- function(d, label) {
  kept <- d %>% group_by(Occup_Group) %>% filter(n() >= MIN_CELLS_PER_OCCUP) %>% ungroup()
  dropped <- setdiff(unique(as.character(d$Occup_Group)),
                     unique(as.character(kept$Occup_Group)))
  if (length(dropped)) {
    cat(sprintf("  %s: dropped %d group(s) with < %d cells: %s\n",
                label, length(dropped), MIN_CELLS_PER_OCCUP,
                paste(dropped, collapse = ", ")))
  }
  kept
}

cat("\nCell-count filter:\n")
iv_wage_data_all   <- drop_thin_groups(iv_wage_data_all,   "All")
iv_wage_data_men   <- drop_thin_groups(iv_wage_data_men,   "Men")
iv_wage_data_women <- drop_thin_groups(iv_wage_data_women, "Women")

# NOTE ON THIN CELLS
#
# Cells are estimated unweighted, which is worth stating plainly because cell
# headcounts are severely uneven. The n_obs column carried through
# build_wage_controls() reports the respondents behind each cell mean, and the
# diagnostic printed below gives the per-group distribution: median cell
# counts run 13-54 in the pooled sample, 4-28 in the sex splits, and 1 for
# military occupations, with a minimum of 1 almost everywhere. Unweighted, a
# one-respondent cell in a small metro carries the same influence as a
# three-thousand-respondent cell in New York, and cell size is correlated with
# the regressor, since large MSAs have both more respondents per cell and
# larger immigrant inflows.
#
# The unevenness is disclosed in the chapter, and the per-cell counts are
# reported alongside the estimates so a reader can weigh them.
# FIXED EFFECTS: METROPOLITAN AREA AND SURVEY YEAR (2026-11-02).
#
# Previously estimated with no fixed effects at all, and the reason recorded
# for that was itself wrong. Both are corrected here.
#
# Card (2001) is a single 1990 cross-section, and his estimating equations
# carry occupation AND city fixed effects: "City fixed effects absorb any
# citywide variables that might otherwise influence the levels of wages or
# employment in the local labor market." Year effects do not arise for him
# because he has no time dimension. Pooling survey years is this chapter's
# extension, so the panel analogue of his design is city and year effects
# together -- which is also the fixed-effect structure of the chapter's own
# education specification. Estimating separately by occupation group handles
# the occupation dimension and precludes an occupation dummy; it does not
# preclude city or year dummies, which is what the earlier note assumed.
#
# WHAT THIS CHANGES. The previously reported positive, significant estimates
# for professional and technical workers and for managers were cross-metro
# comparisons -- high-inflow metros are high-nominal-wage metros, the exact
# correlation Card's city effects exist to remove. Under city and year effects
# on the divided series, five of seven groups are indistinguishable from zero,
# including low-skill labour and service, and the two that survive point in
# opposite directions.
run_wage_iv_by_occup <- function(data) {
  # Names are taken from the SAME filtered list the models come from. The
  # previous version re-split the unfiltered data to build names, so a group
  # dropped by the completeness filter would have shifted every label.
  grps <- data %>%
    group_split(Occup_Group) %>%
    keep(~ nrow(.) > 0 & sum(complete.cases(.)) > 0)
  nms <- map_chr(grps, ~ as.character(unique(.$Occup_Group)))
  set_names(map(grps, function(g) list(
    iv = feols(mean_hrwage ~ mean_age + mean_yrs_sch + frac_black | MSA + YEAR |
                 actual_inflow ~ predicted_inflow, data = g, vcov = ~MSA),
    fs = feols(actual_inflow ~ predicted_inflow + mean_age + mean_yrs_sch +
                 frac_black | MSA + YEAR, data = g, vcov = ~MSA)
  )), nms)
}

iv_wage_models_all   <- run_wage_iv_by_occup(iv_wage_data_all)
iv_wage_models_men   <- run_wage_iv_by_occup(iv_wage_data_men)
iv_wage_models_women <- run_wage_iv_by_occup(iv_wage_data_women)

cat("\n=== Wage IV: All Natives ===\n")
walk(iv_wage_models_all,   ~ print(summary(.$iv)))
cat("\n=== Wage IV: Men ===\n")
walk(iv_wage_models_men,   ~ print(summary(.$iv)))
cat("\n=== Wage IV: Women ===\n")
walk(iv_wage_models_women, ~ print(summary(.$iv)))

# Raw coefficients are reported alongside a per-standard-deviation conversion
# (added 2026-09-21). The two are not competing numbers: the raw estimate is
# the effect of a one-unit change in occupation-specific inflow, but inflow
# volume differs enormously across occupation groups -- military inflow
# averages ~140 with an SD of ~307, against 3,400-14,000 and SDs of
# 12,000-49,000 elsewhere, because non-citizens are largely barred from
# military service. A raw coefficient is therefore measured against a
# different-sized shock in every row, and comparing rows directly invites the
# conclusion that military wages respond hundreds of times more strongly than
# anything else, when most of that gap is the denominator.
#
# Scaling Estimate and SE by the same group-specific SD leaves the
# t-statistic and significance stars unchanged; only the units move. Per-SD
# and the percentage-of-own-mean-wage column are the comparable figures.
compile_wage_table <- function(models_all,   data_all,
                               models_men,   data_men,
                               models_women, data_women) {
  extract_row <- function(models, data, label) {
    # sd_inflow is deliberately UNWEIGHTED: it describes the spread of
    # inflow across cells, so the per-SD column keeps the same units
    # regardless of WAGE_WEIGHT_MODE and stays comparable across runs. It is
    # now the SD of the wave-divided inflow, matching the regressor.
    scales <- data %>%
      group_by(Occup_Group) %>%
      summarise(sd_inflow     = sd(actual_inflow, na.rm = TRUE),
                grp_mean_wage = mean(mean_hrwage,  na.rm = TRUE),
                grp_yrs_sch   = mean(mean_yrs_sch, na.rm = TRUE),
                cells         = n(),
                median_cell_n = median(n_obs, na.rm = TRUE),
                min_cell_n    = min(n_obs,    na.rm = TRUE),
                .groups = "drop") %>%
      mutate(Occup_Group = as.character(Occup_Group))

    map_dfr(names(models), function(occ) {
      # feols carries the MSA-clustered VCOV requested at estimation, so no
      # separate sandwich/lmtest step is needed. The endogenous coefficient is
      # named fit_<endovar>.
      ct    <- models[[occ]]$iv$coeftable["fit_actual_inflow", ]
      fsct  <- models[[occ]]$fs$coeftable["predicted_inflow", ]
      coef  <- ct[["Estimate"]]
      se    <- ct[["Std. Error"]]
      pval  <- ct[["Pr(>|t|)"]]
      stars <- ifelse(pval < 0.01, "***",
               ifelse(pval < 0.05, "**",
               ifelse(pval < 0.10, "*", "")))
      s <- scales[scales$Occup_Group == occ, ]
      data.frame(
        Occup_Group   = occ,
        Group         = label,
        Estimate      = coef,
        SE            = se,
        Stars         = stars,
        SD_inflow     = s$sd_inflow,
        Est_per_SD    = coef * s$sd_inflow,
        SE_per_SD     = se   * s$sd_inflow,
        Pct_wage      = 100 * coef * s$sd_inflow / s$grp_mean_wage,
        FS_t          = fsct[["t value"]],
        FS_F          = fsct[["t value"]]^2,
        Mean_yrs_sch  = s$grp_yrs_sch,
        Cells         = s$cells,
        Median_cell_n = s$median_cell_n,
        Min_cell_n    = s$min_cell_n
      )
    })
  }
  bind_rows(
    extract_row(models_all,   data_all,   "All"),
    extract_row(models_men,   data_men,   "Men"),
    extract_row(models_women, data_women, "Women")
  )
}

wage_results <- compile_wage_table(iv_wage_models_all,   iv_wage_data_all,
                                   iv_wage_models_men,   iv_wage_data_men,
                                   iv_wage_models_women, iv_wage_data_women)

# Reported table: raw estimate and SE (per one-unit inflow), then the same
# estimate rescaled to a one-SD inflow increase, then stars. The per-SD column
# exists because inflow volume differs by an order of magnitude across
# occupation groups -- military inflow averages ~140 against 3,400-14,000
# elsewhere, since non-citizens are largely barred from military service -- so
# raw coefficients are measured against a different-sized shock in every row
# and are not comparable down the column. Rescaling estimate and SE by the
# same group-specific SD leaves the t-statistic and stars unchanged.
cat("\n--- Wage IV by occupation group ---\n")
print(wage_results[, c("Occup_Group", "Group", "Estimate", "SE",
                       "Est_per_SD", "Stars", "FS_t")],
      row.names = FALSE)

# Disclosed as a limitation rather than corrected -- see the note above
# run_wage_iv_by_occup(). Retained so the magnitude of the problem can be
# reported rather than asserted.
# Mean years of schooling by occupation group, for the motivation in the text:
# the mechanism proposed for the education result runs through the return to
# schooling, so which occupation groups are the high-education ones is the fact
# that connects a wage estimate to that channel.
cat("\n--- Mean years of schooling by occupation group ---\n")
print(unique(wage_results[wage_results$Group == "All",
                          c("Occup_Group", "Mean_yrs_sch", "FS_t", "FS_F")]),
      row.names = FALSE, digits = 4)


cat("\n--- Cell composition (respondents behind each cell mean) ---\n")
print(unique(wage_results[, c("Occup_Group", "Group", "Cells",
                              "Median_cell_n", "Min_cell_n")]),
      row.names = FALSE)

# Written to CSV so 24_results_tables.R can render from outputs/tables/ rather
# than from a session object passed down a checkpoint chain. 14 is now a pure
# rendering step with no checkpoint dependency; this is what makes that
# possible for the wage table.
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(wage_results, "outputs/tables/results_wage_iv.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/results_wage_iv.csv\n")

# reduced_data_iv was only needed above to build the occupation-level shift-
# share objects; drop it so it doesn't propagate into Section 14's
# checkpoint (Section 15 needs reduced_data_iv again but reloads its own
# copy separately).
rm(reduced_data_iv)

# -----------------------------------------------------------------------------

# NO CHECKPOINT WRITTEN, as of 2026-10-04.
#
# This wrote after_sec13.RData, which nothing reads. The checkpoint chain was
# dismantled between 2026-09-30 and 2026-10-04 -- every surviving script
# loads analysis_data.RData directly -- and after that the save was writing
# large objects to disk on every run for no consumer. Worse, a stale
# checkpoint left on disk is precisely how 15_rotemberg.R came to decompose
# a superseded estimate in September. The only checkpoint still in the
# Chapter 2 pipeline is after_sec22, written by 21 and read by 24.

cat("\n21_wage_iv.R complete (no checkpoint).\n")
