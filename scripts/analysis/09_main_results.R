# =============================================================================
# 09_main_results.R  --  Chapter 2 headline estimates
# Renumbered 2026-10-21 from 29_main_results.R (archive-and-renumber pass).
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# THE SPECIFICATION
#   outcome    y_HSplus (non_HS vs HighSch/College/Grad) and y_CollPlus
#   instrument avg_exposure_nodiv        predicted inflow, ages 15-25, no divisor
#   endogenous avg_actual_exposure_nodiv actual inflow,    ages 15-25, no divisor
#   sample     natives, top-175 MSAs, NO age restriction
#   controls   Black, SEX, Marital, ChildAgeGroup, POVERTY, Mom_Educ, Dad_Educ,
#              exposure_missing | YEAR + MSA        (NO age fixed effects)
#   inference  MSA-clustered standard errors, with the Adao shock-randomisation
#              placebo (script 20) as the inference test. That placebo is the
#              canonical robustness check for this design (Adao, Kolesar &
#              Morales); alternative procedures are taken up in the later
#              chapter and are deliberately not reported here.
#
# BOTH REGRESSORS NOW COME FROM THE PIPELINE
# 07_exposure_iv.R averages the actual inflow series over the same ages-15-25
# window as the instrument, and 06_shift_share_iv.R emits the undivided series,
# so nothing is reconstructed here. That matters: scripts 24-28 and 18b each
# built the window and reversed the divisor in-script, six copies of logic that
# could drift. This script reads columns.
#
# It also matters that the divisor is reversed in 06 rather than here. avg_*
# variables are means over a window whose years carry DIFFERENT n_waves, so no
# single factor applied to the average can undo the division -- it has to happen
# at MSA x arrival year, before averaging.
#
# WHAT THE CHAPTER CLAIMS, AND WHAT IT DOES NOT
# The point estimate is large -- on the order of 8 percentage points per
# identifying SD -- which is not a plausible magnitude for an immigration effect
# on native attainment. The chapter says so, and points to 18b, where the
# estimate falls from 13.4 to 4.3 pp/SD as the sample is restricted to cohorts
# whose schooling is complete, and to 28, where adding age fixed effects shrinks
# it roughly thirty-fold. The language should therefore describe a direction and
# a mechanism to be resolved, not a precisely estimated magnitude.
#
# A STRUCTURAL FACT TO REPORT ALONGSIDE n
# avg_exposure_nodiv is a deterministic function of (MSA, birth_year). Birth
# years 1975-2007 across 175 metros is roughly 5,800 distinct instrument values
# behind 4.2 million identifying observations, nesting inside 175 clusters and,
# after unequal contribution, about 2.5 effective ones. The observation count is
# not the sample size for this design, and "n = 8.5 million" invites a
# confidence it does not support. The script prints all four figures.
#
# Writes outputs/tables/main_results.csv.
# =============================================================================

# Clears the session but KEEPS bpld_lookup. It is built from the IPUMS DDI in
# 08_load_analysis_subsets.R, is not carried in analysis_data.RData, and cannot
# be rebuilt without re-running 08. A blanket rm(list = ls()) here deletes it,
# so any later script that needs it -- 15, 19, 21, 40, 41 -- fails on its guard
# even though 08 was run at the start of the session.
rm(list = setdiff(ls(), "bpld_lookup")); gc()
source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(fixest)

cat("Loading analysis_data.RData...\n")
load("analysis_data.RData")
setDT(reduced_data_iv)

REQUIRED <- c("avg_exposure_nodiv", "avg_actual_exposure_nodiv",
              "avg_exposure", "avg_actual_exposure", "exposure_missing")
missing_cols <- setdiff(REQUIRED, names(reduced_data_iv))
if (length(missing_cols))
  stop("analysis_data.RData is stale -- missing ", paste(missing_cols, collapse = ", "),
       ".\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R before this script.")

# --- Estimation sample -------------------------------------------------------
d <- reduced_data_iv[immigr == 0 & MSA != 0 & MSA %in% selected_msas]
d[, `:=`(
  birth_year = YEAR - AGE,
  y_HSplus   = as.integer(Educ_Group != "non_HS"),
  y_CollPlus = as.integer(Educ_Group %in% c("College", "Grad")),
  z = avg_exposure_nodiv,
  x = avg_actual_exposure_nodiv
)]
d <- droplevels(d)

id <- d[exposure_missing == 0]
sd_z <- id[, sd(z)]

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SAMPLE\n\n")
cat(sprintf("  observations                 %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  identifying observations     %s\n", format(nrow(id), big.mark = ",")))
cat(sprintf("  distinct instrument values   %s   <- (MSA x birth_year) cells\n",
            format(id[, uniqueN(paste(MSA, birth_year))], big.mark = ",")))
cat(sprintf("  MSA clusters                 %d\n", uniqueN(d$MSA)))
cat(sprintf("  birth years (identifying)    %d-%d\n",
            id[, min(birth_year)], id[, max(birth_year)]))
cat(sprintf("\n  non_HS share                 %.1f%%\n",
            100 * d[, mean(Educ_Group == "non_HS")]))
cat(sprintf("  y_HSplus mean                %.4f\n", d[, mean(y_HSplus)]))
cat(sprintf("  y_CollPlus mean              %.4f\n", d[, mean(y_CollPlus)]))
cat(sprintf("  SD of instrument (identifying) %.1f\n", sd_z))
cat("\n  Note: the observation count is not the sample size for this design.\n")
cat("  Inference is governed by the cluster count, and effectively by fewer.\n")

# --- Machinery ---------------------------------------------------------------
CTRL <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
         Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

# Race-augmented control set, for the robustness row below.
CTRL_RACE <- "Black + Asian + Aboriginal + SEX + Marital + ChildAgeGroup +
              POVERTY + Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

# STATE-CLUSTERED STANDARD ERRORS ADDED 2026-10-05.
#
# This is a reporting addition, not a specification change: the model, sample
# and point estimates are untouched, and only a second variance estimator is
# computed alongside the first. The assertion below enforces that -- clustering
# cannot move a coefficient, so if it does, the two fits differ in some other
# way and the pair is not comparable.
#
# WHY. \citet{adao2019shift} caution that clustering at too fine a level can
# understate sampling variability in shift-share designs, and the chapter needs
# to report the state-to-MSA ratio. It cannot be borrowed from
# 16_adao_placebo.R: 20 runs CTRL_PLA, which excludes the race indicators, and
# more importantly the ratio is NOT one number. Measured 2026-10-05 it is
# 0.880 for the 2SLS on y_HSplus but 0.680 on y_CollPlus, and 0.747 for the
# REDUCED FORM on y_HSplus against 0.880 for the 2SLS on the same outcome. A
# single borrowed figure would describe one cell and misdescribe the other
# five, so each is computed where it is reported.
#
# All are below one, which is the uninformative direction: coarsening from 175
# MSA clusters to 45 state clusters biases cluster-robust variance downward on
# its own, so a smaller state-clustered SE is partly mechanical. A ratio ABOVE
# one would have been the finding that mattered.
estimate <- function(outcome, label, ctrl = CTRL) {
  fs <- feols(as.formula(sprintf("x ~ z + %s", ctrl)), data = d, vcov = ~MSA)
  rf <- feols(as.formula(sprintf("%s ~ z + %s", outcome, ctrl)), data = d, vcov = ~MSA)
  iv <- feols(as.formula(sprintf("%s ~ %s | x ~ z", outcome, ctrl)),
              data = d, vcov = ~MSA)
  ct <- iv$coeftable["fit_x", ]

  rf_st <- feols(as.formula(sprintf("%s ~ z + %s", outcome, ctrl)),
                 data = d, vcov = ~STATEFIP)
  iv_st <- feols(as.formula(sprintf("%s ~ %s | x ~ z", outcome, ctrl)),
                 data = d, vcov = ~STATEFIP)
  ct_st <- iv_st$coeftable["fit_x", ]
  stopifnot(
    isTRUE(all.equal(unname(coef(rf)["z"]), unname(coef(rf_st)["z"]))),
    isTRUE(all.equal(ct[["Estimate"]], ct_st[["Estimate"]])))

  b_rf <- unname(coef(rf)["z"])
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(sprintf("%s\n\n", label))
  cat(sprintf("  FIRST STAGE   %+.6e   SE %.4e   t %6.2f\n",
              coef(fs)["z"], fs$coeftable["z", "Std. Error"],
              fs$coeftable["z", "t value"]))
  cat(sprintf("  REDUCED FORM  %+.6e   SE %.4e   t %6.2f\n",
              b_rf, rf$coeftable["z", "Std. Error"], rf$coeftable["z", "t value"]))
  cat(sprintf("                p %.4f\n", rf$coeftable["z", "Pr(>|t|)"]))
  cat(sprintf("  2SLS          %+.6e   SE %.4e   t %6.2f   p %.4f\n",
              ct[["Estimate"]], ct[["Std. Error"]],
              ct[["t value"]], ct[["Pr(>|t|)"]]))
  cat(sprintf("\n  State-clustered SEs (45 clusters against %d MSA):\n",
              uniqueN(d$MSA)))
  cat(sprintf("    reduced form  SE %.4e   t %6.2f   ratio state/MSA %.4f\n",
              rf_st$coeftable["z", "Std. Error"], rf_st$coeftable["z", "t value"],
              rf_st$coeftable["z", "Std. Error"] / rf$coeftable["z", "Std. Error"]))
  cat(sprintf("    2SLS          SE %.4e   t %6.2f   ratio state/MSA %.4f\n",
              ct_st[["Std. Error"]], ct_st[["t value"]],
              ct_st[["Std. Error"]] / ct[["Std. Error"]]))
  cat(sprintf("\n  Effect size:  %+.4f pp per identifying SD of the instrument\n",
              100 * b_rf * sd_z))

  data.frame(spec = label, outcome = outcome, n = nobs(rf),
             first_stage = unname(coef(fs)["z"]),
             fs_se = fs$coeftable["z", "Std. Error"],
             reduced_form = b_rf, rf_se = rf$coeftable["z", "Std. Error"],
             rf_p = rf$coeftable["z", "Pr(>|t|)"],
             pp_per_sd = 100 * b_rf * sd_z,
             tsls = ct[["Estimate"]], tsls_se = ct[["Std. Error"]],
             tsls_p = ct[["Pr(>|t|)"]],
             rf_se_state   = rf_st$coeftable["z", "Std. Error"],
             rf_ratio      = rf_st$coeftable["z", "Std. Error"] /
                             rf$coeftable["z", "Std. Error"],
             tsls_se_state = ct_st[["Std. Error"]],
             tsls_ratio    = ct_st[["Std. Error"]] / ct[["Std. Error"]],
             stringsAsFactors = FALSE)
}

res <- rbind(
  estimate("y_HSplus",   "HEADLINE: non_HS vs high school or above"),
  estimate("y_CollPlus", "SECONDARY: college or above"),
  # --- Robustness: race controls ------------------------------------------
  # 16_adao_placebo.R finds the instrument predicts Asian identity (14.3% of
  # the reference effect, p = 0.0055) while Black and Aboriginal pass. That is
  # not a hidden confounder. MSA fixed effects absorb the cross-sectional
  # level, so what remains is within-metro variation across cohorts -- and
  # within a gateway metro, the cohorts that grew up during high Asian-
  # immigration years contain more Asian natives, because the US-born children
  # of Asian immigrants ARE natives in exactly those cohorts. The placebo is
  # detecting a definitional property of the sample.
  #
  # The channel through which it could still matter is composition: Asian
  # natives have higher average attainment, so a compositional shift toward
  # them could inflate the estimate. Asian and Aboriginal were absent from the
  # control set (Black was not). Adding them moves the estimate by 0.09%
  # (4.745 -> 4.740 pp/SD), so the composition channel is orthogonal to the
  # effect. Reported here so the obvious referee question is answered in the
  # table rather than in correspondence.
  estimate("y_HSplus", "ROBUSTNESS: headline with Asian and Aboriginal added",
           CTRL_RACE)
)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SUMMARY\n\n")
res$pct_of_headline <- 100 * res$pp_per_sd / res$pp_per_sd[1]
print(res[, c("outcome", "first_stage", "reduced_form", "rf_se", "rf_p",
              "pp_per_sd", "tsls", "tsls_p")],
      row.names = FALSE, digits = 3)

cat("\nSpecifications, with effect size relative to the headline row:\n\n")
print(res[, c("spec", "pp_per_sd", "pct_of_headline")],
      row.names = FALSE, digits = 4)
cat("\n(The pct_of_headline column is only meaningful for rows sharing the\n")
cat("headline's outcome -- y_CollPlus is a different margin, not a robustness\n")
cat("check on the same one.)\n")

cat("\n", strrep("-", 72), "\n", sep = "")
cat("FOR THE WRITE-UP\n\n")
cat("The first stage should be strongly POSITIVE. If it is negative, the\n")
cat("pipeline has regressed to the survey-year endogenous variable -- see the\n")
cat("Step 4b note in 07_exposure_iv.R.\n\n")
cat("Report the Adao shock-randomisation placebo (script 20) as the inference\n")
cat("test. It is the canonical robustness check for shift-share designs and\n")
cat("the one a referee in this literature will expect.\n\n")
cat("The magnitude remains large for an immigration effect on native\n")
cat("attainment and the chapter should say so, then point to 18b (the cohort\n")
cat("gradient) and 28 (the age-fixed-effects comparison) as the hand-off to\n")
cat("the age-structure chapter.\n\n")
cat("Both of those scripts must be RE-RUN before their figures are quoted.\n")
cat("Their previously recorded numbers were computed on the construction that\n")
cat("preceded the two 07_exposure_iv.R fixes of 2026-09-29 -- the seq()\n")
cat("direction bug, which admitted ~713,000 observations of children aged\n")
cat("0-15 carrying a single year's inflow as their formative exposure, and\n")
cat("the AGE_15_25 branch, which admitted pre-1975 cohorts whose windows\n")
cat("predated the 1990 share year. Together those moved this headline from\n")
cat("8.19 to 4.74 pp per SD.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/main_results.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/main_results.csv\n")
cat("09_main_results.R complete.\n")
