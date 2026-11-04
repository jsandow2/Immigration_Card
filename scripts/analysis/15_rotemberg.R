# =============================================================================
# 15_rotemberg.R  --  Section 19: Rotemberg weight decomposition (GPS 2020)
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 18_rotemberg.R.
# Renumbered 2026-10-21 from 19_rotemberg.R (archive-and-renumber pass).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

# UNCHAINED 2026-10-04. This loaded after_sec18.RData, which is written by
# 18_birth_cohort_check.R -- a CHAPTER 3 script. A Chapter 2 headline result
# should not depend on a Chapter 3 checkpoint, and on a clean rebuild that file
# does not exist, so this script would have failed to source.
#
# The load was already vestigial. after_sec18 carries only all_natives_175 and
# bpld_lookup; the refresh below overwrote the first immediately, and the second
# comes from 08_load_analysis_subsets.R regardless. Nothing is lost by dropping
# it and the guard below makes the one real dependency explicit.
if (!exists("bpld_lookup"))
  stop("bpld_lookup not found. Run 08_load_analysis_subsets.R first -- it is\n",
       "built from the IPUMS DDI and is not in analysis_data.RData.")

# --- Pipeline objects from analysis_data.RData (added 2026-09-29) -----------
#
# Loaded fresh rather than inherited. This is a stale-checkpoint hazard worth
# being explicit about: because the final specification uses NEW variable
# names, a stale object errors loudly against the guards below. Had it reused
# total_actual_inflow, this script would have run to completion on superseded
# data and produced a correctly-validated decomposition of the wrong estimate.
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
all_natives_175  <- .tmp_env$all_natives_175
shift_share_data <- .tmp_env$shift_share_data
rm(.tmp_env); gc()

.need <- c("avg_exposure_nodiv", "avg_actual_exposure_nodiv",
           "age_15_year", "age_25_year", "exposure_missing")
.miss <- setdiff(.need, names(all_natives_175))
if (length(.miss))
  stop("analysis_data.RData is stale -- all_natives_175 is missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R first.")
.miss <- setdiff(c("raw_count", "predicted_inflow_nodiv", "share"),
                 names(shift_share_data))
if (length(.miss))
  stop("analysis_data.RData is stale -- shift_share_data is missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R first.")
rm(.need, .miss)
cat("Refreshed all_natives_175 and shift_share_data from analysis_data.RData\n")

# Section 19: Rotemberg weight decomposition
#
# Decomposes the shift-share IV estimate into weighted averages of
# just-identified IV estimators, one per origin group k, following
# Goldsmith-Pinkham, Sorkin, and Swift (2020, AER).
#
# Runs on a 500k random sample from all_natives_175. The allow.cartesian
# merge in Step 2 requires ~4 GB of RAM at peak. On a 32 GB machine,
# earlier sections of the pipeline consume ~20 GB, leaving ~12 GB free
# which is sufficient. The cleanup block below removes objects no longer
# needed to ensure adequate headroom.
# -----------------------------------------------------------------------------

# Memory cleanup -- remove large objects from earlier sections no longer
# needed. Under the current explicit-save() checkpoint chain this list will
# mostly find nothing to remove (those objects were already dropped from
# each section's own checkpoint rather than carried forward) -- kept as a
# defensive no-op in case this script is ever run in a continuous session
# where they're still resident in memory from an earlier section.
cat("Freeing memory before Section 19...\n")
objs_to_remove <- c(
  # Wage model objects
  "iv_wage_data_all", "iv_wage_data_men", "iv_wage_data_women",
  "iv_wage_models_all", "iv_wage_models_men", "iv_wage_models_women",
  "actual_inflows_occ_year", "shift_share_occ_data",
  "predicted_inflows_occ", "national_inflows_occ",
  "immigrant_shares_occ", "msa_bpld_occ_year_grid",
  "wage_results", "wage_wide",
  # Migration model objects
  "diff_iv_dt", "diff_iv_values", "bl_inflows", "cr_inflows",
  "msa_year_full", "migr_all", "migr_men", "migr_women",
  "samp_1m_all", "samp_500k_men", "samp_500k_women",
  # Figure objects
  "plot_data", "immig_educ_prop", "native_educ_prop",
  "immig_college_grad", "native_college_grad", "educ_prop_plot",
  "plot_data_age25plus", "educ_prop_plot_age25plus",
  "plot_data_age18_24", "educ_prop_plot_age18_24",
  "immig_age25plus", "native_age25plus", "immig_age18_24",
  "build_educ_prop_data", "plot_educ_prop", "educ_colors",
  # Meta-analysis intermediate objects
  "comparison_results", "meta_table",
  # Section 17 objects (large subsample, no longer needed)
  "samp_em0", "iv_em0", "fs_em0", "fs_em0_red", "em_counts",
  # Demographic subsets no longer needed
  "native_men", "native_women", "sample_dt"
)
objs_to_remove <- objs_to_remove[sapply(objs_to_remove, exists)]
if (length(objs_to_remove) > 0) {
  rm(list = objs_to_remove)
  cat(sprintf("  Removed %d objects\n", length(objs_to_remove)))
}
gc()
cat("  Memory freed. Proceeding with Section 19.\n\n")

cat(sprintf("\n%s\n", strrep("=", 60)))
cat("Section 19: Rotemberg Weight Decomposition\n")
cat(sprintf("%s\n\n", strrep("=", 60)))

# Outcome for the decomposition (set 2026-09-24).
#
# The Rotemberg weights are built from residualised quantities, so they are
# outcome-specific, and Step 5 reports the sensitivity as a percent change
# against the baseline estimate. Run on Educ_Binary -- whose pooled estimate was
# indistinguishable from zero -- those percentages divide by approximately zero
# and are not interpretable. y_HSplus carries the detectable effect and gives a
# stable denominator. Set to "Educ_Binary" to reproduce the earlier
# decomposition.
ROTEMBERG_OUTCOME <- "y_HSplus"

# MUST be constructed BEFORE the sample is drawn (moved 2026-09-29). samp is a
# copy, so adding the outcome to all_natives_175 afterwards does not reach it.
# This block previously sat after the draw and never fired, because the
# all_natives_175 carried by after_sec18.RData already had y_HSplus from an
# earlier section. Now that the object is refreshed from analysis_data.RData --
# where 07_exposure_iv.R does not build outcome variables -- it does fire, and
# the ordering mattered.
setDT(all_natives_175)
if (!ROTEMBERG_OUTCOME %in% names(all_natives_175)) {
  if (ROTEMBERG_OUTCOME == "y_HSplus") {
    all_natives_175[, y_HSplus := as.integer(Educ_Group != "non_HS")]
  } else if (ROTEMBERG_OUTCOME == "Educ_Binary") {
    all_natives_175[, Educ_Binary := as.integer(Educ_Group %in% c("College", "Grad"))]
  } else {
    stop("Outcome not found on all_natives_175 and no rule to build it: ",
         ROTEMBERG_OUTCOME)
  }
}
cat(sprintf("Rotemberg decomposition outcome: %s (mean %.4f)\n",
            ROTEMBERG_OUTCOME, all_natives_175[[ROTEMBERG_OUTCOME]] |> mean()))

# Draw 500k sample (also used by Section 20)
set.seed(42)
n_sample   <- 500000
sample_idx <- sample(nrow(all_natives_175), n_sample)
samp       <- as.data.frame(all_natives_175[sample_idx, ]) %>% droplevels()
stopifnot(ROTEMBERG_OUTCOME %in% names(samp),
          "avg_exposure_nodiv"        %in% names(samp),
          "avg_actual_exposure_nodiv" %in% names(samp))

cat(sprintf("Sample size: %d observations\n", nrow(samp)))
cat(sprintf("MSAs in sample: %d\n", n_distinct(samp$MSA)))
cat(sprintf("Years in sample: %d\n\n", n_distinct(samp$YEAR)))

# Step 1: Estimate main IV on sample
cat("Step 1: Estimating main IV on sample...\n")

# ROTEMBERG_OUTCOME and its defensive construction now sit above the sample
# draw -- see the note there. Moved 2026-09-29.

# FINAL SPECIFICATION (updated 2026-09-29)
#
# Endogenous variable is avg_actual_exposure_nodiv -- actual inflow averaged
# over the individual's ages-15-25 window, undivided. It was previously
# total_actual_inflow, the count of arrivals in the SURVEY year, which is the
# timing mismatch documented at Step 4b of 07_exposure_iv.R and which inverts
# the first stage.
#
# This matters more here than almost anywhere else in the project. The
# decomposition check sum_k alpha_k * beta_k = beta_IV is an ALGEBRAIC identity
# in whatever endogenous variable is supplied -- it passes regardless. So a
# clean validation is not evidence the decomposition is of the right estimate.
# Both alpha_k and beta_k are built from x-perp, so feeding the survey-year
# measure produces a correctly-validated decomposition of the wrong quantity.
# CONVERTED TO fixest::feols 2026-09-29.
#
# Previously beta_IV came from AER::ivreg while alpha_k and beta_k came from
# residuals produced by SEPARATE lm() calls. The decomposition identity
#
#     sum_k alpha_k * beta_k = (sum_k z_k)' y-perp / (sum_k z_k)' x-perp
#
# only reproduces beta_IV if those two routes residualise on exactly the same
# column space. Nothing enforced that, and any difference -- a dropped factor
# level, a differing NA set after droplevels(), collinearity resolved
# differently -- surfaced as a discrepancy in the check. On 2026-09-29 it came
# to 11.87% even though sum_k z_ik reproduced the instrument exactly (ratio 1,
# sd 3e-16), which is what isolated the residualisation as the cause.
#
# Residualising ONCE with feols and deriving beta_IV from the same vectors
# makes the identity hold by construction. The check then tests the one thing
# it should -- that sum_k z_ik equals the instrument -- rather than failing for
# reasons that have nothing to do with the decomposition.
#
# feols also absorbs the 175 MSA and 19 year fixed effects rather than
# estimating 193 dummy coefficients, which is considerably faster on 500k rows.
CTRL_ROT <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
             Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

resid_on_ctrl <- function(v, dat)
  residuals(feols(as.formula(paste(v, "~", CTRL_ROT)), data = dat))

samp <- as.data.table(samp)
samp[, `:=`(
  resid_z = resid_on_ctrl("avg_exposure_nodiv",        samp),
  resid_x = resid_on_ctrl("avg_actual_exposure_nodiv", samp),
  resid_y = resid_on_ctrl(ROTEMBERG_OUTCOME,           samp)
)]

# Just-identified 2SLS from the residualised vectors. Identical to what a
# feols IV call returns, and it is the quantity the decomposition sums to.
beta_iv <- samp[, sum(resid_z * resid_y) / sum(resid_z * resid_x)]
cat(sprintf("IV estimate (beta_IV): %.6e\n", beta_iv))

# Cross-check against the packaged IV estimator. These should agree to
# floating-point; a gap means feols and the residual route disagree about the
# control space, which is exactly the failure this conversion removes.
iv_mod <- feols(as.formula(sprintf(
  "%s ~ %s | avg_actual_exposure_nodiv ~ avg_exposure_nodiv",
  ROTEMBERG_OUTCOME, CTRL_ROT)), data = samp, vcov = ~MSA)
beta_iv_pkg <- iv_mod$coeftable["fit_avg_actual_exposure_nodiv", "Estimate"]
cat(sprintf("  feols IV cross-check:  %.6e  (relative gap %.3e)\n",
            beta_iv_pkg, abs(beta_iv_pkg - beta_iv) / abs(beta_iv)))
if (abs(beta_iv_pkg - beta_iv) / abs(beta_iv) > 1e-6)
  warning("Residual-route and feols IV estimates disagree -- the control ",
          "space is not identical between them.")

# First stage, on the same control space.
fs_mod <- feols(as.formula(paste("avg_actual_exposure_nodiv ~ avg_exposure_nodiv +",
                                 CTRL_ROT)), data = samp, vcov = ~MSA)
cat(sprintf("First stage: %+.6e  (t %.1f)\n",
            coef(fs_mod)["avg_exposure_nodiv"],
            fs_mod$coeftable["avg_exposure_nodiv", "t value"]))
cat(sprintf("Correlation(instrument, endogenous) after residualising: %.4f\n",
            samp[, cor(resid_z, resid_x)]))
cat("  A high value here means a STRONG first stage, which is what a shift-share\n")
cat("  design wants. It does not imply 2SLS sits near OLS: the OLS-2SLS gap is\n")
cat("  the OLS bias, cov(x,e)/var(x), which does not depend on first-stage\n")
cat("  strength. Instrument and endogenous variable are averaged over the same\n")
cat("  window, so they differ only by the 1990 share structure -- that\n")
cat("  difference is the identifying content, not a defect.\n\n")

# Step 2: Construct z_{i,k} via MSA x cohort aggregation
cat("Step 2: Constructing z_{i,k}...\n")

# shift_share_data and all_natives_175 are both refreshed from
# analysis_data.RData at the top of this script, before the 500k sample is
# drawn -- see the note there. The reload that used to sit here has been
# removed to avoid loading a large file twice.
setDT(shift_share_data)
setDT(all_natives_175)

# UNDIVIDED national inflow (updated 2026-09-29). raw_count is the pooled
# respondent count before the eligible-waves division, and 06_shift_share_iv.R
# now carries it through shift_share_data. Previously this recovered the
# DIVIDED series via predicted_inflow / share; z_ik must be on the same scale
# as the instrument that produced beta_iv, which is now avg_exposure_nodiv, so
# the shocks have to be undivided too. Taking raw_count directly is exact and
# avoids the division (predicted_inflow_nodiv / share == raw_count).
national_inflow_dt <- shift_share_data[share > 0][, .(
  national_inflow = mean(raw_count, na.rm = TRUE)
), by = .(BPLD, YEAR)]
setkey(national_inflow_dt, BPLD, YEAR)

shares_dt <- shift_share_data[YEAR == min(shift_share_data$YEAR),
                               .(MSA, BPLD, share)]
shares_dt[, MSA := as.character(MSA)]
setkey(shares_dt, MSA, BPLD)

bpld_groups <- unique(shift_share_data$BPLD)
cat(sprintf("Number of origin groups (BPLD): %d\n\n", length(bpld_groups)))

samp_key <- as.data.table(samp)[, .(
  PERSON_ID,
  MSA           = as.character(MSA),
  birth_year,
  age_15_year,
  age_25_year   = pmin(age_25_year, 2022L),
  exposure_missing,
  resid_z, resid_x, resid_y
)]

msa_birth <- unique(samp_key[, .(MSA, birth_year, age_15_year, age_25_year)])
msa_birth[, cohort_id := .I]

cat("  Expanding cohort x BPLD combinations...\n")
t0 <- proc.time()

cohort_bpld <- msa_birth[
  rep(seq_len(.N), each = length(bpld_groups))
][, BPLD := rep(bpld_groups, times = nrow(msa_birth))]

cohort_bpld <- merge(cohort_bpld, national_inflow_dt,
                     by = "BPLD", allow.cartesian = TRUE)
cohort_bpld <- cohort_bpld[YEAR >= age_15_year & YEAR <= age_25_year]
cohort_bpld <- cohort_bpld[, .(
  mean_national_inflow = mean(national_inflow, na.rm = TRUE)
), by = .(cohort_id, MSA, birth_year, BPLD)]

cohort_bpld <- merge(cohort_bpld, shares_dt,
                     by = c("MSA", "BPLD"), all.x = TRUE)
cohort_bpld[is.na(share), share := 0]
# No /100 (removed 2026-09-24), matching the removal of the Section 8 rescale.
# z_ik must be on the same scale as the avg_exposure that produced beta_iv,
# since alpha_k divides one by the other. Scaling z_ik alone rescaled every
# weight by a constant -- harmless for the rankings and the top-N shares, which
# are ratios, but it put the printed alpha_k levels and the "Sum of weights"
# diagnostic off by two orders of magnitude. Rotemberg weights should sum to
# approximately 1; treat that printed line as the check that this is right.
cohort_bpld[, z_ik := share * mean_national_inflow]

cat(sprintf("  Cohort x BPLD done (%.1fs)\n\n", (proc.time() - t0)["elapsed"]))

samp_key <- merge(samp_key,
                  msa_birth[, .(MSA, birth_year, cohort_id)],
                  by = c("MSA", "birth_year"))

indiv_bpld <- merge(
  samp_key[, .(PERSON_ID, cohort_id, exposure_missing, resid_x, resid_y)],
  cohort_bpld[, .(cohort_id, BPLD, z_ik)],
  by = "cohort_id", allow.cartesian = TRUE)

# Zero the shock components for exposure_missing == 1 (added 2026-09-24).
#
# The decomposition identity beta = sum_k alpha_k * beta_k holds only if the
# per-origin components sum to the instrument that actually entered the
# regression: sum_k z_ik must equal avg_exposure_nodiv for every individual.
#
# Step 2 above rebuilds an exposure window for everyone from age_15_year to
# min(age_25_year, 2022), but 07_exposure_iv.R assigns exposure = 0 to anyone
# who fails the exposed_post_2000 condition -- roughly half the sample without
# the age restriction, whose ages-15-25 window begins before the 1990 share
# year. Without this correction those individuals carry positive z_ik against
# a zero instrument.
#
# Verified 2026-09-24 on the 500k sample under the then-current construction:
# for exposure_missing == 0, sum_k z_ik reproduced the instrument exactly
# (correlation 1, means both 266.6136); for exposure_missing == 1, mean
# instrument was 0 against mean sum_k z_ik of 112.78, which is why the
# decomposition check returned -1.42e-04 against a beta_IV of +9.52e-05.
#
# RE-VERIFY AFTER THE 2026-09-29 RESPECIFICATION. The instrument is now
# avg_exposure_nodiv and the shocks are built from raw_count, so the levels
# above no longer apply. The check below must be re-run and the correlation
# confirmed at 1 before any of this output is used -- the decomposition
# identity will validate whether or not the scales match, so it cannot catch
# a mismatch here on its own.
indiv_bpld[exposure_missing == 1L, z_ik := 0]

# Step 3: Compute Rotemberg weights
#
# REWRITTEN 2026-09-24. The previous implementation computed
#
#     alpha_k = sum(z_ik * resid_second) / (n_indiv * beta_iv)
#
# using SECOND-STAGE residuals. That is not the Goldsmith-Pinkham, Sorkin and
# Swift construction, and it failed both of that construction's defining
# properties: the weights did not sum to one (they summed to 6,301 under
# y_HSplus and 31,387 under Educ_Binary), and they changed completely with the
# outcome -- Mexico ranked first under one outcome and outside the top twenty
# under the other. Rotemberg weights describe the INSTRUMENT, so they cannot
# depend on the dependent variable.
#
# The correct decomposition is beta_Bartik = sum_k alpha_k * beta_k, with
#
#     alpha_k = (z_k' x-perp) / sum_k (z_k' x-perp)
#     beta_k  = (z_k' y-perp) / (z_k' x-perp)
#
# where x-perp and y-perp are the endogenous regressor and the outcome
# residualised on the controls and fixed effects. alpha_k sums to one by
# construction and involves no outcome; beta_k is the just-identified estimate
# using origin k's shock alone. Their weighted sum reproduces beta_IV, which is
# checked explicitly below and is a far stronger diagnostic than the old
# sum-of-weights line.
cat("Step 3: Computing Rotemberg weights...\n")
t0 <- proc.time()

n_indiv <- nrow(samp)

rotemberg <- indiv_bpld[, .(
  zx       = sum(z_ik * resid_x, na.rm = TRUE),
  zy       = sum(z_ik * resid_y, na.rm = TRUE),
  sum_z_ik = sum(z_ik, na.rm = TRUE),
  n_msas   = n_distinct(cohort_id)
), by = BPLD]

rotemberg[, `:=`(
  alpha_k = zx / sum(zx),
  beta_k  = fifelse(zx != 0, zy / zx, NA_real_)
)]
rotemberg[, alpha_k_norm := alpha_k / sum(abs(alpha_k))]
rotemberg <- rotemberg[order(-abs(alpha_k))]

# Validation: sum_k alpha_k * beta_k must reproduce beta_IV.
beta_recomposed <- rotemberg[!is.na(beta_k), sum(alpha_k * beta_k)]
cat(sprintf("  Decomposition check: sum(alpha_k * beta_k) = %.6e vs beta_IV = %.6e\n",
            beta_recomposed, beta_iv))
cat(sprintf("  Relative discrepancy: %.4f%%\n",
            100 * (beta_recomposed - beta_iv) / abs(beta_iv)))
cat("  (A discrepancy of more than a fraction of a percent means the\n")
cat("   decomposition does not reproduce the estimate and nothing below\n")
cat("   should be reported.)\n")

cat(sprintf("Rotemberg weights computed (%.1fs)\n\n", (proc.time() - t0)["elapsed"]))

# Step 4: Diagnostics
cat(sprintf("%s\n", strrep("=", 60)))
cat("ROTEMBERG WEIGHT DECOMPOSITION RESULTS\n")
cat(sprintf("%s\n\n", strrep("=", 60)))

cat(sprintf("IV estimate (beta_IV):       %.6e\n", beta_iv))
cat(sprintf("First stage:                 %+.6e (t %.1f)\n",
            coef(fs_mod)["avg_exposure_nodiv"],
            fs_mod$coeftable["avg_exposure_nodiv", "t value"]))
cat(sprintf("Number of origin groups (k): %d\n", nrow(rotemberg)))
cat(sprintf("Sum of weights:              %.4f  (should be 1.0000)\n",
            sum(rotemberg$alpha_k)))
cat(sprintf("Negative-weight groups:      %d of %d (%.1f%% of |alpha_k|)\n\n",
            rotemberg[alpha_k < 0, .N], nrow(rotemberg),
            100 * rotemberg[alpha_k < 0, sum(abs(alpha_k))] /
              rotemberg[, sum(abs(alpha_k))]))

for (top_n in c(1, 5, 10, 20)) {
  top_share <- sum(abs(rotemberg$alpha_k[1:min(top_n, nrow(rotemberg))]))
  cat(sprintf("Top %2d groups account for %.1f%% of |alpha_k|\n",
              top_n, top_share * 100 / sum(abs(rotemberg$alpha_k))))
}

# BPLD labels now sourced from the authoritative IPUMS codebook (bpld_lookup,
# loaded in Section 8 from usa_00007.xml) rather than a manual mapping. The
# manual vector this replaces mislabeled several codes (e.g. "51100" was
# tagged "India" but is actually Cambodia; India is 52100) and used a
# different regional code scheme than the actual BPLD values in the data.
# Kept as a named vector so all downstream indexing (bpld_labels[bpld_code])
# below is unchanged.
bpld_labels <- setNames(as.character(bpld_lookup$label),
                        as.character(bpld_lookup$BPLD))

cat(sprintf("\n--- Top 20 origin groups by |Rotemberg weight| ---\n"))
cat(sprintf("%-10s  %-34s  %-11s  %-11s  %-12s\n",
            "BPLD", "Description", "alpha_k", "alpha_k_norm", "beta_k"))
cat(strrep("-", 80), "\n")

top20 <- head(rotemberg, 20)
for (i in seq_len(nrow(top20))) {
  bpld_code <- as.character(top20$BPLD[i])
  label     <- ifelse(bpld_code %in% names(bpld_labels),
                      bpld_labels[bpld_code], paste("BPLD", bpld_code))
  cat(sprintf("%-10s  %-34s  %-11.4f  %-11.4f  %-12.4e\n",
              bpld_code, label, top20$alpha_k[i], top20$alpha_k_norm[i],
              top20$beta_k[i]))
}

# Step 5: Sensitivity -- drop top-weighted groups
cat(sprintf("\n--- Sensitivity: IV estimate dropping top origin groups ---\n"))
cat(sprintf("%-30s  %-12s  %-10s\n", "Specification", "IV estimate", "% change"))
cat(strrep("-", 55), "\n")

top_bplds <- head(rotemberg$BPLD, 10)
# Saved as well as printed, added 2026-10-05. These figures were console-only,
# which meant the chapter had to assert the sensitivity rather than show it.
sens_rows <- list()

for (drop_n in c(1, 3, 5, 10)) {
  drop_bplds   <- head(top_bplds, drop_n)
  ss_reduced   <- shift_share_data[!BPLD %in% drop_bplds]
  # predicted_inflow_nodiv, not predicted_inflow (updated 2026-09-29). The
  # reduced instrument has to sit on the same scale as avg_exposure_nodiv,
  # which produced beta_iv; otherwise the percent changes below compare two
  # differently-normalised instruments rather than measuring the effect of
  # dropping origin groups.
  msa_yr_red   <- ss_reduced[, .(
    total_predicted_inflow_reduced = sum(predicted_inflow_nodiv, na.rm = TRUE)
  ), by = .(MSA, YEAR)]
  msa_yr_red[, MSA := as.character(MSA)]
  setkey(msa_yr_red, MSA, YEAR)

  samp_red <- as.data.table(samp)[, .(
    PERSON_ID, MSA = as.character(MSA),
    birth_year, age_15_year,
    age_25_year = pmin(age_25_year, 2022L)
  )]

  samp_exp <- samp_red[, .(YEAR = seq(age_15_year, age_25_year)),
                       by = .(PERSON_ID, MSA)]
  samp_exp <- merge(samp_exp, msa_yr_red, by = c("MSA", "YEAR"), all.x = TRUE)
  samp_exp[is.na(total_predicted_inflow_reduced),
           total_predicted_inflow_reduced := 0]

  # No /100 (removed 2026-09-24). This mirrored the old Section 8 rescale,
  # which no longer exists, so dividing here put avg_exposure_reduced on a
  # different scale from the avg_exposure used for beta_iv -- making the
  # percent changes below a comparison between two differently-scaled
  # instruments rather than a sensitivity check.
  avg_exp_red <- samp_exp[, .(
    avg_exposure_reduced = mean(total_predicted_inflow_reduced, na.rm = TRUE)
  ), by = PERSON_ID]

  samp_sens <- merge(as.data.table(samp), avg_exp_red,
                     by = "PERSON_ID", all.x = TRUE)
  samp_sens[is.na(avg_exposure_reduced), avg_exposure_reduced := 0]

  # Zero the reduced instrument for exposure_missing == 1 (added 2026-09-24),
  # for the same reason as in Step 2. This loop rebuilds a window for everyone
  # from age_15_year to min(age_25_year, 2022), but 07_exposure_iv.R zeroes
  # avg_exposure for the 82.5% failing the exposed_post_2000 condition. Without
  # this, beta_sens is estimated on an instrument that is non-zero for those
  # individuals while beta_iv used one that is zero for them -- so the
  # comparison conflated dropping origin groups with un-zeroing most of the
  # sample, and reported a sign flip that the validated decomposition (19 of
  # the top 20 origins carrying positive beta_k) contradicts.
  samp_sens[exposure_missing == 1L, avg_exposure_reduced := 0]

  iv_sens <- tryCatch({
    feols(as.formula(sprintf(
      "%s ~ %s | avg_actual_exposure_nodiv ~ avg_exposure_reduced",
      ROTEMBERG_OUTCOME, CTRL_ROT)), data = samp_sens)
  }, error = function(e) NULL)

  if (!is.null(iv_sens)) {
    beta_sens  <- iv_sens$coeftable["fit_avg_actual_exposure_nodiv", "Estimate"]
    se_sens    <- iv_sens$coeftable["fit_avg_actual_exposure_nodiv", "Std. Error"]
    pct_change <- (beta_sens - beta_iv) / abs(beta_iv) * 100
    cat(sprintf("Drop top %-2d groups (BPLD: %-20s)  %-12.6e  %+.1f%%\n",
                drop_n, paste(head(drop_bplds, 3), collapse = ","),
                beta_sens, pct_change))
    sens_rows[[length(sens_rows) + 1L]] <- data.frame(
      drop_n        = drop_n,
      dropped_bplds = paste(drop_bplds, collapse = ";"),
      weight_dropped = sum(rotemberg$alpha_k[seq_len(drop_n)]),
      beta_full     = beta_iv,
      beta_reduced  = beta_sens,
      se_reduced    = se_sens,
      pct_change    = pct_change,
      stringsAsFactors = FALSE)
  } else {
    cat(sprintf("Drop top %-2d groups: FAILED\n", drop_n))
  }
}

# Reference row first, so the table carries the full-sample estimate it is
# measured against rather than requiring the reader to look it up.
rotemberg_sensitivity <- rbind(
  data.frame(drop_n = 0L, dropped_bplds = "", weight_dropped = 0,
             beta_full = beta_iv, beta_reduced = beta_iv, se_reduced = NA_real_,
             pct_change = 0, stringsAsFactors = FALSE),
  if (length(sens_rows)) rbindlist(sens_rows) else NULL)
write.csv(rotemberg_sensitivity, "outputs/tables/rotemberg_sensitivity.csv",
          row.names = FALSE)
cat("\nSaved: outputs/tables/rotemberg_sensitivity.csv\n")

# Step 6: F-statistic contribution by origin group
cat(sprintf("\n--- First-stage strength by origin group: top 10 ---\n"))
cat(sprintf("%-10s  %-40s  %-12s\n", "BPLD", "Description", "t on z_ik"))
cat(strrep("-", 65), "\n")

for (i in 1:min(10, nrow(rotemberg))) {
  k         <- rotemberg$BPLD[i]
  bpld_code <- as.character(k)
  label     <- ifelse(bpld_code %in% names(bpld_labels),
                      bpld_labels[bpld_code], paste("BPLD", bpld_code))

  single_z <- indiv_bpld[BPLD == k, .(PERSON_ID, z_ik)]
  samp_k   <- merge(as.data.table(samp), single_z,
                    by = "PERSON_ID", all.x = TRUE)
  samp_k[is.na(z_ik), z_ik := 0]

  # Converted to feols 2026-09-29, on the same control space as the main
  # estimate. Reported as the t-statistic on the single shock rather than an
  # F: with one restriction F = t^2, and the squared form produced figures in
  # the millions that were hard to compare and easy to misread as a joint
  # first-stage F for the full instrument.
  t_k <- tryCatch({
    m_k <- feols(as.formula(paste("avg_actual_exposure_nodiv ~ z_ik +", CTRL_ROT)),
                 data = samp_k, vcov = ~MSA)
    m_k$coeftable["z_ik", "t value"]
  }, error = function(e) NA_real_)

  cat(sprintf("%-10s  %-40s  %-12.1f\n", bpld_code, label, t_k))
}

# Step 7: Save results
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(rotemberg), "outputs/tables/rotemberg_results.csv", row.names = FALSE)

top20_out <- as.data.frame(top20) %>%
  mutate(label = ifelse(as.character(BPLD) %in% names(bpld_labels),
                        bpld_labels[as.character(BPLD)],
                        paste("BPLD", BPLD))) %>%
  select(BPLD, label, alpha_k, alpha_k_norm, sum_z_ik, n_msas)

write.csv(top20_out, "outputs/tables/rotemberg_top_groups.csv", row.names = FALSE)

cat(sprintf("\n%s\n", strrep("=", 60)))
cat("Section 19 complete.\n")
cat(sprintf("%s\n", strrep("=", 60)))

# shift_share_data was only needed above for the Rotemberg decomposition;
# drop it so it doesn't propagate into Section 20's checkpoint (Section 21
# reloads its own copy of the objects it needs separately).
rm(shift_share_data)

# -----------------------------------------------------------------------------

# NO CHECKPOINT WRITTEN, as of 2026-10-04.
#
# This used to save all_natives_175, bpld_lookup and samp (the 500k draw) to
# after_sec19.RData, justified by a comment saying "samp is explicitly reused
# by Section 20". That is no longer true. 16_adao_placebo.R was rebuilt
# self-sufficient on 2026-09-29 and now draws its own sample from
# analysis_data.RData; both it and 20b carry comments explaining that they
# stopped reading this checkpoint precisely because it carried
# pre-respecification objects. Nothing loads it.
#
# Writing a 500k-row object to disk for no consumer is not free, and a stale
# checkpoint that nothing reads is exactly the hazard this script's own header
# warns about -- so it goes. 19 now neither loads nor saves a checkpoint.
cat("\nNo checkpoint written -- nothing downstream reads one from this script.\n")
cat("15_rotemberg.R complete.\n")
