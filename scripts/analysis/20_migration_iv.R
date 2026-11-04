# =============================================================================
# 20_migration_iv.R  --  Section 12: Individual-level migration models (Card 2001 + diff_iv)
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 13_migration_iv.R.
# Renumbered 2026-10-21 from 12_migration_iv.R (archive-and-renumber pass).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

# CHECKPOINT LOAD REMOVED 2026-09-30.
#
# This previously loaded after_sec11.RData, whose contents it then overrode --
# see the note above the analysis_data.RData reload below, which states the
# section deliberately does not use the checkpoint's all_natives_175. The only
# other objects the checkpoint carried (bpld_lookup, top175_list_named,
# meta_summary_full, multi_boot_summary) are never referenced here; they were
# passed straight through to the next save.
#
# With Sections 10-11 retired, after_sec11.RData is no longer regenerated on a
# clean rebuild, so the load would fail for objects this script does not use.
# It reads analysis_data.RData directly instead.
#
# The specification here is unaffected by the 2026-09-30 respecification of the
# education design. This script uses total_actual_inflow and
# total_predicted_inflow at the SURVEY year, which is correct for its question:
# whether contemporaneous immigrant inflows induce contemporaneous out-migration.
# The timing mismatch fixed in 07_exposure_iv.R was specific to the ages-15-25
# exposure design, where the instrument described adolescence while the
# regressor described the survey year. Here the survey year is the period of
# interest, so avg_exposure, avg_actual_exposure and exposure_missing do not
# enter and nothing this script uses changed.

# Section 12: Individual-level migration models
#
# Estimates the causal effect of immigrant inflows on the probability that a
# U.S.-born individual moves out of their MSA of residence. Two specifications:
#
# Spec 1 (Card 2001 replication):
#   Instrument: total_predicted_inflow (MSA x year shift-share)
#   Fixed effects: YEAR only; MSA included as direct control
#   Sample: ages 16-68 (Card's own range), then a 1M random draw (all) /
#   500k (men/women) -- MSA as a direct control requires an explicit dummy
#   matrix, which exceeds RAM on the full sample
#   Clustered: MSA (state clustering only if significant at 10%)
#
# Spec 2 (diff_iv -- novel instrument):
#   Instrument: diff_iv = cr - bl where
#     bl = mean predicted_inflow in individual's MSA during ages 14-19
#     cr = mean predicted_inflow in individual's MSA over last 5 years
#   Fixed effects: YEAR + MSA (instrument varies within MSA x year cells
#   by birth cohort, solving the collinearity problem of Spec 1)
#   Sample: full ages 25-68 sample. The 25 floor is required by the
#   instrument itself -- below AGE 24 the ages-14-19 baseline window and the
#   last-5-years current window overlap and diff_iv collapses toward zero.
#   See the age-range block above the diff_iv construction for the algebra.
#   Clustered: MSA and state
#
# Controls (both specs, following Card 2001):
#   AGE, AGE^2, YearSch, Black, Asian, Aboriginal, SEX, Occup_Group
#
# Key identification insight for diff_iv:
#   Although the individual-level sample is observed only in census/ACS years,
#   the national inflow series {Delta_M_{k,t}} is reconstructed for all years
#   from immigrants' self-reported year of arrival (YRIMMIG). This allows
#   predicted inflows to be computed for any year within the historical range
#   of shift_share_data (1910-2023), enabling valid baseline windows for all
#   U.S.-born individuals aged 25+ regardless of birth cohort. Note this is
#   a separate question from the window-overlap floor: the series range is
#   what makes old cohorts usable, the overlap algebra is what sets the
#   lower bound at 25.
# -----------------------------------------------------------------------------

# Helper: run individual-level migration IV using feols() from fixest
# `specs` selects which specifications to estimate: "card" (Spec 1) and/or
# "diff" (Spec 2). The wrapper below asks for exactly one of each, because
# Spec 1 is reported from the random draw and Spec 2 from the full sample --
# estimating both on both datasets meant a costly full-sample Spec 1 fit that
# was computed and then discarded, plus a guaranteed "diff_iv not found"
# failure on the Spec 1 draw, which does not carry diff_iv.
run_migration_iv <- function(data, group_name, specs = c("card", "diff")) {

  cat(sprintf("\n=== Migration IV: %s (n = %d) ===\n",
              group_name, nrow(data)))

  d <- data %>%
    # The outcome is built as MOVED (uppercase) in 03_variable_construction.R
    # section 3.9; this script was written against a lowercase `moved` that
    # nothing in the pipeline ever creates, so every model formula below
    # referenced a column that did not exist. Aliased here rather than
    # renamed throughout so the formulas stay as written. MOVED is integer
    # 0/1 with default = 0L, so the is.na() filter below is a no-op kept for
    # safety. (Pre-existing bug, found 2026-09-19.)
    mutate(moved = as.integer(MOVED)) %>%
    filter(!is.na(moved),
           !is.na(Occup_Group),
           !is.na(YearSch)) %>%
    mutate(
      YEAR        = as.integer(as.character(YEAR)),
      MSA         = as.character(MSA),
      Occup_Group = as.character(Occup_Group),
      Black       = as.integer(Black),
      Asian       = as.integer(Asian),
      Aboriginal  = as.integer(Aboriginal),
      SEX         = as.integer(SEX)
    )

  # Guard the MOVED -> moved coercion above. If MOVED ever becomes a factor
  # upstream, as.integer() would silently return level codes (1/2) instead of
  # values (0/1) and every LPM below would be estimated on a corrupted
  # outcome. Fail loudly here rather than let that through.
  stopifnot(all(d$moved %in% c(0L, 1L)))

  cat(sprintf("  Observations after filtering: %d\n", nrow(d)))
  cat(sprintf("  Mean P(moved): %.4f\n", mean(d$moved)))

  # --- Specification 1: Card (2001) replication ---
  mod_card <- NULL
  if ("card" %in% specs) {
    cat("  Running Specification 1: Card replication (MSA as control)...\n")

    mod_card <- tryCatch(
      fixest::feols(
        moved ~ AGE + I(AGE^2) + YearSch + Black + Asian +
          Aboriginal + SEX + Occup_Group + MSA |
          YEAR |
          total_actual_inflow ~ total_predicted_inflow,
        data  = d,
        vcov  = ~MSA
      ),
      error = function(e) {
        cat(sprintf("  Spec 1 FAILED: %s\n", conditionMessage(e)))
        NULL
      }
    )
  }

  if (!is.null(mod_card)) {
    est_card  <- coef(mod_card)["fit_total_actual_inflow"]
    se_card   <- se(mod_card)["fit_total_actual_inflow"]
    t_card    <- tstat(mod_card)["fit_total_actual_inflow"]
    p_card    <- pvalue(mod_card)["fit_total_actual_inflow"]
    # CHANGED 2026-10-21: "ivf" -> "ivwald". fitstat's "ivf" is the
    # HOMOSKEDASTIC first-stage F. These models cluster on MSA, and the
    # instrument varies at the metro-year level, so a statistic computed as
    # though errors were iid across a million individual records is close to a
    # restatement of the sample size -- it was returning 2.37e6. "ivwald" uses
    # the model's own VCOV, so the statistic is the clustered Wald test on the
    # excluded instrument and is comparable to the standard errors reported
    # beside it. The tryCatch is gone with it: it mapped any failure to NA and
    # a silent NA in a strength diagnostic is worse than an error.
    fs_F_card <- fitstat(mod_card, "ivwald")[[1]]$stat

    cat(sprintf("  Spec 1 (Card, MSA-clustered):\n"))
    cat(sprintf("    Estimate: %.6e | SE: %.6e | t: %.3f | p: %.4f\n",
                est_card, se_card, t_card, p_card))
    cat(sprintf("    First-stage Wald (MSA-clustered): %.1f\n", fs_F_card))

    mod_card_state <- NULL
    if (!is.na(p_card) && p_card < 0.10) {
      cat("  Significant at 10% -- running state-clustered robustness...\n")
      mod_card_state <- tryCatch(
        fixest::feols(
          moved ~ AGE + I(AGE^2) + YearSch + Black + Asian +
            Aboriginal + SEX + Occup_Group + MSA |
            YEAR |
            total_actual_inflow ~ total_predicted_inflow,
          data  = d,
          vcov  = ~STATEFIP
        ),
        error = function(e) NULL
      )
    } else {
      cat("  Not significant -- skipping state-clustered robustness.\n")
    }
  } else {
    est_card <- se_card <- t_card <- p_card <- fs_F_card <- NA
    mod_card_state <- NULL
  }

  # --- Specification 2: diff_iv ---
  if ("diff" %in% specs) {
    cat("\n  Running Specification 2: diff_iv instrument...\n")
  }

  if ("diff" %in% specs &&
      "diff_iv" %in% names(d) && sum(!is.na(d$diff_iv)) > 1000) {

    mod_diff_msa <- tryCatch(
      fixest::feols(
        moved ~ AGE + I(AGE^2) + YearSch + Black + Asian +
          Aboriginal + SEX + Occup_Group |
          YEAR + MSA |
          total_actual_inflow ~ diff_iv,
        data  = d,
        vcov  = ~MSA
      ),
      error = function(e) {
        cat(sprintf("  Spec 2 MSA FAILED: %s\n", conditionMessage(e)))
        NULL
      }
    )

    mod_diff_state <- tryCatch(
      fixest::feols(
        moved ~ AGE + I(AGE^2) + YearSch + Black + Asian +
          Aboriginal + SEX + Occup_Group |
          YEAR + MSA |
          total_actual_inflow ~ diff_iv,
        data  = d,
        vcov  = ~STATEFIP
      ),
      error = function(e) {
        cat(sprintf("  Spec 2 State FAILED: %s\n", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(mod_diff_msa)) {
      est_diff      <- coef(mod_diff_msa)["fit_total_actual_inflow"]
      se_diff       <- se(mod_diff_msa)["fit_total_actual_inflow"]
      t_diff        <- tstat(mod_diff_msa)["fit_total_actual_inflow"]
      p_diff        <- pvalue(mod_diff_msa)["fit_total_actual_inflow"]
      # Same change as for the Card specification above.
      fs_F_diff     <- fitstat(mod_diff_msa, "ivwald")[[1]]$stat
      se_diff_state <- if (!is.null(mod_diff_state))
        se(mod_diff_state)["fit_total_actual_inflow"] else NA

      cat(sprintf("  Spec 2 (diff_iv, MSA-clustered):\n"))
      cat(sprintf("    Estimate: %.6e | SE: %.6e | t: %.3f | p: %.4f\n",
                  est_diff, se_diff, t_diff, p_diff))
      cat(sprintf("    State SE: %.6e\n", se_diff_state))
      cat(sprintf("    First-stage Wald (MSA-clustered): %.1f\n", fs_F_diff))
    } else {
      est_diff <- se_diff <- t_diff <- p_diff <- fs_F_diff <- NA
      se_diff_state <- NA
      mod_diff_state <- NULL
    }

  } else {
    # Only a genuine problem if Spec 2 was actually requested; silent when
    # this call was asked for "card" only.
    if ("diff" %in% specs) {
      cat("  diff_iv not found -- run diff_iv construction block first.\n")
    }
    est_diff <- se_diff <- t_diff <- p_diff <- fs_F_diff <- NA
    se_diff_state <- NA
    mod_diff_msa <- mod_diff_state <- NULL
  }

  fs_F    <- fs_F_card
  est     <- est_card
  msa_se  <- se_card
  msa_t   <- t_card
  msa_p   <- p_card
  state_se <- if (!is.null(mod_card_state))
    se(mod_card_state)["fit_total_actual_inflow"] else NA
  state_t  <- if (!is.null(mod_card_state))
    tstat(mod_card_state)["fit_total_actual_inflow"] else NA
  state_p  <- if (!is.null(mod_card_state))
    pvalue(mod_card_state)["fit_total_actual_inflow"] else NA

  cat(sprintf("\n  --- Summary: %s ---\n", group_name))
  if (!is.na(est)) {
    cat(sprintf("  Spec 1 (Card): est = %.6e | SE = %.6e | t = %.3f | p = %.4f\n",
                est, msa_se, msa_t, msa_p))
  }
  if (!is.na(est_diff)) {
    cat(sprintf("  Spec 2 (diff_iv): est = %.6e | SE = %.6e | t = %.3f | p = %.4f\n",
                est_diff, se_diff, t_diff, p_diff))
  }

  wh <- tryCatch(fitstat(mod_card, "wh"), error = function(e) NULL)
  if (!is.null(wh)) {
    cat(sprintf("  Wu-Hausman p (Spec 1): %.4f\n", wh[[1]]$p))
  }

  list(
    group          = group_name,
    n_obs          = nrow(d),
    mean_moved     = mean(d$moved),
    estimate       = est,
    msa_se         = msa_se,
    msa_t          = msa_t,
    msa_p          = msa_p,
    state_se       = state_se,
    state_t        = state_t,
    state_p        = state_p,
    first_stage_F  = fs_F_card,
    mod_card       = mod_card,
    mod_card_state = mod_card_state,
    est_diff       = est_diff,
    se_diff        = se_diff,
    t_diff         = t_diff,
    p_diff         = p_diff,
    se_diff_state  = se_diff_state,
    fs_F_diff      = fs_F_diff,
    mod_diff_msa   = mod_diff_msa,
    mod_diff_state = mod_diff_state
  )
}

# --- Construct diff_iv from full shift_share_data annual series ---
cat("\n--- Constructing diff_iv instrument ---\n")

# shift_share_data is not carried through the checkpoint chain past Section 8
# (see 08_load_analysis_subsets.R) -- reload it fresh here directly from
# analysis_data.RData, via a throwaway environment so nothing else already
# in memory (e.g. a processed all_natives_175) gets overwritten.
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
shift_share_data     <- .tmp_env$shift_share_data
natives_unrestricted <- .tmp_env$all_natives_175
rm(.tmp_env)

# -----------------------------------------------------------------------------
# Migration-specific age ranges (added 2026-09-18)
#
# This section deliberately does NOT use the checkpoint's all_natives_175.
# That object carries the ESTIMATION_MIN_AGE <- 30L floor applied in
# 08_load_analysis_subsets.R, which was justified by an argument specific to
# EDUCATION outcomes: the ages-15-25 exposure window must have closed and
# attainment must have settled before the outcome is read. The outcome here
# is `moved`, so that rationale does not transfer, and inheriting the floor
# silently cost ~350-400k observations and pushed Spec 1 further from the
# Card (2001) sample it is replicating.
#
# Spec 1 uses Card's own 16-68 range. Its instrument is the contemporaneous
# MSA x YEAR shift-share (total_predicted_inflow), which involves no
# individual age window, so nothing constrains the lower bound.
#
# Spec 2 keeps a 25+ floor because diff_iv = cr - bl needs its two windows
# to stay disjoint. The baseline ends at birth_year + 19 = YEAR - AGE + 19
# and the current window starts at YEAR - 4, so they separate only while
#
#     YEAR - AGE + 19  <  YEAR - 4   =>   AGE > 23
#
# At AGE <= 23 the windows overlap and diff_iv collapses mechanically toward
# zero; by 18 they are nearly identical years. 25 clears that bound with a
# year of margin. The bl_n < 3 guard below does NOT catch this -- an
# 18-year-old in 2022 passes it and still gets a near-zero diff_iv.
#
# The 68 upper bound is a comparability choice rather than a technical one
# (the 14-19 window is computable for anyone born after ~1896 given the
# series starts in 1910); it keeps both specs on the same age support and
# matches Card. all_natives_175 itself is left untouched here so it reaches
# Section 13 as the education-floored estimation sample it is elsewhere.
# -----------------------------------------------------------------------------
MIGRATION_AGE_MIN_SPEC1 <- 16L
MIGRATION_AGE_MIN_SPEC2 <- 25L
MIGRATION_AGE_MAX       <- 68L

natives_spec1 <- natives_unrestricted %>%
  filter(AGE >= MIGRATION_AGE_MIN_SPEC1, AGE <= MIGRATION_AGE_MAX)
natives_spec2 <- natives_unrestricted %>%
  filter(AGE >= MIGRATION_AGE_MIN_SPEC2, AGE <= MIGRATION_AGE_MAX)

cat(sprintf("  Spec 1 sample (AGE %d-%d): %d obs\n",
            MIGRATION_AGE_MIN_SPEC1, MIGRATION_AGE_MAX, nrow(natives_spec1)))
cat(sprintf("  Spec 2 sample (AGE %d-%d): %d obs\n",
            MIGRATION_AGE_MIN_SPEC2, MIGRATION_AGE_MAX, nrow(natives_spec2)))
rm(natives_unrestricted)

msa_year_full <- shift_share_data %>%
  group_by(MSA = as.character(MSA), YEAR) %>%
  summarize(total_predicted_inflow = sum(predicted_inflow, na.rm = TRUE),
            .groups = "drop")

cat(sprintf("  Years available: %d to %d\n",
            min(msa_year_full$YEAR), max(msa_year_full$YEAR)))

cat("  Computing baseline and current windows...\n")
t0 <- proc.time()

diff_iv_dt <- natives_spec2 %>%
  mutate(
    MSA_chr  = as.character(MSA),
    YEAR_int = as.integer(as.character(YEAR)),
    bl_start = birth_year + 14L,
    bl_end   = birth_year + 19L,
    cr_start = YEAR_int - 4L,
    cr_end   = YEAR_int
  ) %>%
  select(PERSON_ID, MSA_chr, YEAR_int, bl_start, bl_end, cr_start, cr_end)

bl_inflows <- diff_iv_dt %>%
  left_join(msa_year_full %>% rename(MSA_chr = MSA),
            by = "MSA_chr", relationship = "many-to-many") %>%
  filter(YEAR >= bl_start & YEAR <= bl_end) %>%
  group_by(PERSON_ID) %>%
  summarize(bl = mean(total_predicted_inflow, na.rm = TRUE),
            bl_n = n(), .groups = "drop") %>%
  mutate(bl = ifelse(bl_n < 3, NA_real_, bl))

cr_inflows <- diff_iv_dt %>%
  left_join(msa_year_full %>% rename(MSA_chr = MSA),
            by = "MSA_chr", relationship = "many-to-many") %>%
  filter(YEAR >= cr_start & YEAR <= cr_end) %>%
  group_by(PERSON_ID) %>%
  summarize(cr = mean(total_predicted_inflow, na.rm = TRUE),
            cr_n = n(), .groups = "drop") %>%
  mutate(cr = ifelse(cr_n < 3, NA_real_, cr))

diff_iv_values <- diff_iv_dt %>%
  select(PERSON_ID) %>%
  left_join(bl_inflows %>% select(PERSON_ID, bl), by = "PERSON_ID") %>%
  left_join(cr_inflows %>% select(PERSON_ID, cr), by = "PERSON_ID") %>%
  mutate(diff_iv = cr - bl)

cat(sprintf("  diff_iv computed in %.1fs\n", (proc.time() - t0)["elapsed"]))
cat(sprintf("  Non-missing: %d of %d in Spec 2 sample (%.1f%%)\n",
            sum(!is.na(diff_iv_values$diff_iv)),
            nrow(diff_iv_values),
            100 * mean(!is.na(diff_iv_values$diff_iv))))

# No /100 rescale (removed 2026-09-18): that convention dated from the
# PERWT-era inflow scale, the same leftover stripped out of Section 8.
# total_predicted_inflow now comes from shift_share_data under the raw-count
# + eligible-waves-divisor construction in 06_shift_share_iv.R and is already
# on a small scale, so diff_iv stays in raw units.
if ("diff_iv" %in% names(natives_spec2)) {
  natives_spec2 <- natives_spec2 %>% select(-diff_iv)
}
natives_spec2 <- natives_spec2 %>%
  left_join(diff_iv_values %>% select(PERSON_ID, diff_iv), by = "PERSON_ID")

cat("  diff_iv merged (raw units, no rescale)\n\n")

# --- Draw samples for Spec 1 ---
cat("--- Drawing samples for Spec 1 ---\n")
set.seed(42)
samp_1m_all <- as.data.frame(
  natives_spec1[sample(nrow(natives_spec1), min(1000000, nrow(natives_spec1))), ]
) %>% droplevels()

set.seed(42)
samp_500k_men <- as.data.frame(
  natives_spec1 %>% filter(SEX == 1) %>%
    .[sample(nrow(.), min(500000, nrow(.))), ]
) %>% droplevels()

set.seed(42)
samp_500k_women <- as.data.frame(
  natives_spec1 %>% filter(SEX == 2) %>%
    .[sample(nrow(.), min(500000, nrow(.))), ]
) %>% droplevels()

cat(sprintf("  All sample: %d | Men: %d | Women: %d\n\n",
            nrow(samp_1m_all), nrow(samp_500k_men), nrow(samp_500k_women)))

# --- Wrapper: Spec 1 on sample, Spec 2 on full data ---
run_migration_both_specs <- function(samp, full_data, group_name,
                                     full_group_name) {
  # Spec 1 is reported from the random draw, Spec 2 from the full sample --
  # see the note above run_migration_iv(). Only the fields actually harvested
  # below are estimated, so nothing costly gets computed and thrown away.
  result_samp <- run_migration_iv(samp,      group_name,      specs = "card")
  result_full <- run_migration_iv(full_data, full_group_name, specs = "diff")
  list(
    group          = group_name,
    n_obs_spec1    = result_samp$n_obs,
    n_obs_spec2    = result_full$n_obs,
    mean_moved     = result_samp$mean_moved,
    estimate       = result_samp$estimate,
    msa_se         = result_samp$msa_se,
    msa_t          = result_samp$msa_t,
    msa_p          = result_samp$msa_p,
    state_se       = result_samp$state_se,
    first_stage_F  = result_samp$first_stage_F,
    mod_card       = result_samp$mod_card,
    est_diff       = result_full$est_diff,
    se_diff        = result_full$se_diff,
    t_diff         = result_full$t_diff,
    p_diff         = result_full$p_diff,
    se_diff_state  = result_full$se_diff_state,
    fs_F_diff      = result_full$fs_F_diff,
    mod_diff_msa   = result_full$mod_diff_msa,
    mod_diff_state = result_full$mod_diff_state
  )
}

# --- Run models ---
cat("--- Running migration IV models ---\n")
t0 <- proc.time()

migr_all <- run_migration_both_specs(
  samp_1m_all,    natives_spec2,
  "All U.S.-born (Spec1: 1M, age 16-68)", "All U.S.-born (Spec2: full, age 25-68)")

migr_men <- run_migration_both_specs(
  samp_500k_men,  natives_spec2 %>% filter(SEX == 1),
  "U.S.-born Men (Spec1: 500k, age 16-68)", "U.S.-born Men (Spec2: full, age 25-68)")

migr_women <- run_migration_both_specs(
  samp_500k_women, natives_spec2 %>% filter(SEX == 2),
  "U.S.-born Women (Spec1: 500k, age 16-68)", "U.S.-born Women (Spec2: full, age 25-68)")

cat(sprintf("\nSection 12 total elapsed: %.1fs\n",
            (proc.time() - t0)["elapsed"]))

# --- Summary and save ---
migr_summary <- bind_rows(lapply(
  list(migr_all, migr_men, migr_women),
  function(r) {
    if (is.null(r)) return(NULL)
    data.frame(
      Group         = r$group,
      N_Spec1       = r$n_obs_spec1,
      N_Spec2       = r$n_obs_spec2,
      Mean_P_moved  = round(r$mean_moved, 4),
      Est_Card      = r$estimate,
      SE_Card_MSA   = r$msa_se,
      p_Card        = r$msa_p,
      FS_F_Card     = r$first_stage_F,
      Est_diff      = r$est_diff,
      SE_diff_MSA   = r$se_diff,
      SE_diff_State = r$se_diff_state,
      p_diff        = r$p_diff,
      FS_F_diff     = r$fs_F_diff
    )
  }
))

print(migr_summary)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(migr_summary, "outputs/tables/results_migration_iv.csv", row.names = FALSE)
cat("Migration IV results saved to outputs/tables/results_migration_iv.csv\n")

# shift_share_data was only needed above to build msa_year_full / diff_iv;
# drop it so it doesn't propagate into Section 13's checkpoint (Section 13
# needs reduced_data_iv, not shift_share_data, and reloads that separately).
rm(shift_share_data)

# -----------------------------------------------------------------------------

# migr_summary is already written to outputs/tables/results_migration_iv.csv
# above and nothing downstream reads it back as an object, so it's not
# carried forward -- passthrough of Section 11's objects only.
# NO CHECKPOINT WRITTEN, as of 2026-10-04.
#
# This wrote after_sec12.RData, which nothing reads. The checkpoint chain was
# dismantled between 2026-09-30 and 2026-10-04 -- every surviving script
# loads analysis_data.RData directly -- and after that the save was writing
# large objects to disk on every run for no consumer. Worse, a stale
# checkpoint left on disk is precisely how 15_rotemberg.R came to decompose
# a superseded estimate in September. The only checkpoint still in the
# Chapter 2 pipeline is after_sec22, written by 21 and read by 24.

cat("\n20_migration_iv.R complete (no checkpoint).\n")
