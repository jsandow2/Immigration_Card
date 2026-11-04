# =============================================================================
# 10_exposure_missing_check.R  --  Should the estimation sample include
#                                  observations that cannot identify anything?
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 17b_exposure_missing_check.R.
# Renumbered 2026-10-21 from 17_exposure_missing_check.R (archive-and-renumber pass).
# REBUILT 2026-09-30 on the final specification.
#
# THIS IS NO LONGER A MINOR ROBUSTNESS CHECK
# It was written as one, when the question was whether the main result survived
# restricting to the identifying subsample and the answer was "inconclusive."
# 11_birth_cohort_final_spec.R changed that. Restricting to observations that
# can contribute identifying variation moves the reduced form from 7.25e-06 to
# 1.94e-05 -- from 4.74 to 12.72 percentage points per SD, a factor of 2.7.
# This script is therefore choosing between two candidate headlines, not
# confirming one.
#
# WHAT THE TWO SAMPLES ARE
# `exposure_missing == 1` means avg_exposure is exactly zero by construction:
# the individual's ages-15-25 window opens before 1990, the year the shift-share
# shares are measured, so no predicted inflow can be assigned without using
# shares that postdate the treatment. See the WHY 1990 block in
# 07_exposure_iv.R. Without an age restriction that is 59.2% of the sample.
#
# WHY INCLUDING THEM CHANGES THE ESTIMATE AT ALL
# It is not that they contribute identifying variation -- they cannot. Their
# instrument is zero, and because `exposure_missing` is itself a control, the
# residualised instrument is approximately zero for them too, so they add
# almost nothing to either the numerator or the denominator of
# (z~'y~)/(z~'x~).
#
# They change the estimate by changing what "residualised" MEANS. The YEAR and
# MSA fixed effects and every control coefficient are fitted on the pooled
# sample, so 4.3 million non-identifying observations help determine the
# baseline that the identifying observations are then measured against.
#
# AND THAT IS THE SUBSTANTIVE QUESTION, because the two groups differ
# systematically on exactly those controls. The non-identifying group is older,
# has completed its schooling, and -- since Mom_Educ and Dad_Educ come from
# household parental linkage (03_variable_construction.R lines 268-282) -- is
# almost entirely "Missing" on parental education. So the control coefficients
# and fixed effects are being estimated substantially on a population unlike
# the one the estimate is about.
#
#   ARGUMENT FOR POOLING: the YEAR and MSA effects are more precisely estimated
#   with the full native population, and those effects are properties of the
#   metro-year, not of the cohort.
#
#   ARGUMENT FOR RESTRICTING: the control coefficients should be fitted on the
#   population the estimate describes. Returns to parental education, marital
#   status and poverty plainly differ between a 55-year-old with completed
#   schooling and a 25-year-old, and pooling imposes common slopes.
#
# The script reports the covariate divergence between the two groups so the
# choice can be made on evidence rather than on preference.
#
# Inference: MSA-clustered standard errors. No bootstrap -- that is the later
# chapter's material; the Adao placebo (script 20) is Chapter 2's inference test.
#
# Writes outputs/tables/exposure_missing_check.csv.
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(fixest)

cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
all_natives_175 <- as.data.table(.tmp_env$all_natives_175)
rm(.tmp_env); gc()

.need <- c("avg_exposure_nodiv", "avg_actual_exposure_nodiv",
           "exposure_missing", "Educ_Group", "AGE", "YEAR")
.miss <- setdiff(.need, names(all_natives_175))
if (length(.miss))
  stop("analysis_data.RData is stale or incomplete -- missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R first.")
rm(.need, .miss)

d <- all_natives_175[MSA != 0]
d[, `:=`(
  birth_year = YEAR - AGE,
  y_HSplus   = as.integer(Educ_Group != "non_HS"),
  y_CollPlus = as.integer(Educ_Group %in% c("College", "Grad")),
  z = avg_exposure_nodiv,
  x = avg_actual_exposure_nodiv
)]
d <- droplevels(d)

# One SD for both samples, taken from the identifying observations, so the
# per-SD figures are directly comparable across rows rather than each being
# scaled by its own sample's SD.
sd_z <- d[exposure_missing == 0, sd(z)]

cat(sprintf("\n%s\n", strrep("=", 68)))
cat("Section 17: should non-identifying observations be in the sample?\n")
cat(sprintf("%s\n\n", strrep("=", 68)))

cat(sprintf("  full sample            %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  exposure_missing == 0  %s  (%.1f%%)\n",
            format(d[exposure_missing == 0, .N], big.mark = ","),
            100 * d[exposure_missing == 0, .N] / nrow(d)))
cat(sprintf("  exposure_missing == 1  %s  (%.1f%%)\n",
            format(d[exposure_missing == 1, .N], big.mark = ","),
            100 * d[exposure_missing == 1, .N] / nrow(d)))
cat(sprintf("\n  SD of instrument among identifying observations: %.1f\n", sd_z))
cat(sprintf("  SD among non-identifying: %.1f  (zero by construction)\n",
            d[exposure_missing == 1, sd(z)]))

# --- How different are the two groups on the controls? ----------------------
cat("\n", strrep("-", 68), "\n", sep = "")
cat("COVARIATE DIVERGENCE\n")
cat("If the two groups are similar, pooling is innocuous. If they differ on\n")
cat("the controls, pooling imposes common slopes across populations that do\n")
cat("not share them.\n\n")

div <- d[, .(
  n            = .N,
  mean_age     = mean(AGE),
  birth_min    = min(birth_year),
  birth_max    = max(birth_year),
  pct_non_HS   = 100 * mean(Educ_Group == "non_HS"),
  pct_collplus = 100 * mean(y_CollPlus),
  pct_mom_miss = if ("Mom_Educ" %in% names(d)) 100 * mean(Mom_Educ == "Missing") else NA_real_,
  pct_dad_miss = if ("Dad_Educ" %in% names(d)) 100 * mean(Dad_Educ == "Missing") else NA_real_,
  pct_black    = 100 * mean(Black == 1)
), by = .(exposure_missing)][order(exposure_missing)]
print(div, row.names = FALSE, digits = 4)

# --- Estimation ---------------------------------------------------------------
# exposure_missing is a control in the pooled sample and is dropped in the
# restricted one, where it is constant by construction. feols reports the
# collinearity removal; that message is expected, not a warning of a problem.
CTRL_POOLED <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
                Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"
CTRL_RESTRICTED <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
                    Mom_Educ + Dad_Educ | YEAR + MSA"

run <- function(dat, ctrl, outcome, label) {
  dat <- droplevels(dat)
  fs <- feols(as.formula(sprintf("x ~ z + %s", ctrl)), data = dat, vcov = ~MSA)
  rf <- feols(as.formula(sprintf("%s ~ z + %s", outcome, ctrl)),
              data = dat, vcov = ~MSA)
  iv <- feols(as.formula(sprintf("%s ~ %s | x ~ z", outcome, ctrl)),
              data = dat, vcov = ~MSA)
  ct <- iv$coeftable["fit_x", ]
  b  <- unname(coef(rf)["z"])

  cat(sprintf("\n  %-34s n %10s\n", label, format(nrow(dat), big.mark = ",")))
  cat(sprintf("    first stage   %+.6e  (t %6.2f)\n",
              coef(fs)["z"], fs$coeftable["z", "t value"]))
  cat(sprintf("    reduced form  %+.6e  SE %.4e  t %6.2f  p %.4f\n",
              b, rf$coeftable["z", "Std. Error"],
              rf$coeftable["z", "t value"], rf$coeftable["z", "Pr(>|t|)"]))
  cat(sprintf("    2SLS          %+.6e  SE %.4e  t %6.2f  p %.4f\n",
              ct[["Estimate"]], ct[["Std. Error"]],
              ct[["t value"]], ct[["Pr(>|t|)"]]))
  cat(sprintf("    %+.3f pp per identifying SD\n", 100 * b * sd_z))

  data.frame(sample = label, outcome = outcome, n = nrow(dat),
             # Base rate of the outcome in this sample. A percentage-point
             # effect is only interpretable against the share who have the
             # outcome at all, and the two samples differ sharply here.
             outcome_mean = 100 * mean(dat[[outcome]]),
             first_stage = unname(coef(fs)["z"]),
             reduced_form = b, rf_se = rf$coeftable["z", "Std. Error"],
             rf_p = rf$coeftable["z", "Pr(>|t|)"],
             pp_per_sd = 100 * b * sd_z,
             tsls = ct[["Estimate"]], tsls_se = ct[["Std. Error"]],
             tsls_p = ct[["Pr(>|t|)"]], stringsAsFactors = FALSE)
}

rows <- list()
for (oc in c("y_HSplus", "y_CollPlus")) {
  cat("\n", strrep("-", 68), "\n", sep = "")
  cat(sprintf("%s\n", oc))
  rows[[length(rows)+1]] <- run(d, CTRL_POOLED, oc,
                                "pooled (all observations)")
  rows[[length(rows)+1]] <- run(d[exposure_missing == 0], CTRL_RESTRICTED, oc,
                                "restricted (identifying only)")
}
res <- rbindlist(rows)

cat("\n", strrep("=", 68), "\n", sep = "")
cat("SUMMARY\n\n")
print(res[, .(sample, outcome, n, first_stage, reduced_form, rf_p, pp_per_sd)],
      row.names = FALSE, digits = 3)

cat("\nRatio of restricted to pooled, per outcome:\n")
for (oc in unique(res$outcome)) {
  p <- res[sample %like% "pooled"     & outcome == oc, pp_per_sd]
  r <- res[sample %like% "restricted" & outcome == oc, pp_per_sd]
  cat(sprintf("  %-12s pooled %+7.3f  restricted %+7.3f  ratio %.2fx\n",
              oc, p, r, r / p))
}

cat("\n", strrep("-", 68), "\n", sep = "")
cat("READING THIS\n\n")
cat("The two samples are not a robustness pair where agreement is reassuring\n")
cat("and disagreement is a problem. They are two defensible specifications,\n")
cat("and the covariate table above is the evidence for choosing between them.\n\n")
cat("If the groups are close on the controls, pooling is innocuous and the\n")
cat("larger sample is preferable for estimating the fixed effects.\n\n")
cat("If they diverge -- and in particular if parental education is almost\n")
cat("entirely missing in the non-identifying group while substantially\n")
cat("observed in the identifying one -- then the pooled specification is\n")
cat("fitting control slopes on one population and applying them to another.\n")
cat("In that case the restricted sample is the better specification, and the\n")
cat("chapter should report it with the pooled figure as the robustness row\n")
cat("rather than the reverse.\n\n")
cat("Whichever is chosen, say which and why. They differ by roughly a factor\n")
cat("of three, so this is not a detail that can be left implicit in a table\n")
cat("note.\n\n")
cat("WHY NO AGE RESTRICTION IS OFFERED HERE. The identifying requirement and\n")
cat("the sample's youth are the same constraint: birth_year >= 1975 follows\n")
cat("from the 1990 share year, and it is what makes the identifying group\n")
cat("young (mean age 21.4, non_HS 43.2% against 21.5% among the excluded).\n")
cat("An age cut would not resolve that -- it would trade one selected sample\n")
cat("for another, and age is not the quantity of interest anyway. Someone who\n")
cat("left school at sixteen has completed their schooling; someone in a\n")
cat("doctoral programme at twenty-eight has not. What is unobserved is the\n")
cat("FINAL state, and no threshold on age recovers it.\n\n")
cat("Card's filter is not portable either. His 16-68 sample with at least one\n")
cat("year of potential experience conditions on years of schooling -- fine\n")
cat("when the outcome is wages, but here schooling IS the outcome, so the\n")
cat("filter would select on the dependent variable.\n\n")
cat("So this is a property of the design to document, not to filter. The base\n")
cat("year that buys share exogeneity is the same choice that concentrates\n")
cat("identification in cohorts whose attainment is still being realised. That\n")
cat("belongs in the limitations and in the hand-off to the later chapter,\n")
cat("stated plainly, with the mechanism named.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/exposure_missing_check.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/exposure_missing_check.csv\n")

# NO CHECKPOINT WRITTEN, as of 2026-10-04.
#
# This used to save after_sec17.RData solely to keep a chain alive:
# 18_birth_cohort_check.R loaded it, and 15_rotemberg.R loaded 18's output in
# turn. The note that stood here said both could go once 18 was archived.
# That has now happened -- 19 was unchained on 2026-10-04 and 18 archived the
# same day -- so nothing reads this file and the save goes with it.
#
# 18 was archived because its question is already answered on the final
# specification by 14_gender_difference_test.R. Birth cohort is YEAR - AGE,
# so with YEAR fixed effects already in the model, age and cohort fixed effects
# span the same space: 18's cohort-FE test and 31b's age-FE test are the same
# test. 31b runs it without 18's exposure_missing == 0 restriction, which 18
# needed only to dodge a collinearity the age formulation does not have.
cat("\nNo checkpoint written -- nothing downstream reads one from this script.\n")
cat("10_exposure_missing_check.R complete.\n")
