# =============================================================================
# 17_asian_placebo_concentration.R  --  Where does the Asian placebo failure
# Renumbered 2026-10-21 from 20b_asian_placebo_concentration.R (archive-and-renumber pass).
#                                        come from, and is it where the real
#                                        result comes from?
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# THE QUESTION
# 16_adao_placebo.R finds Asian identity significantly predicted by the
# instrument while Black and Aboriginal pass cleanly. On the final
# specification (2026-09-29) the figures are t = +2.814 MSA-clustered, +3.601
# state-clustered, with the effect equal to 14.3% of the reference y_HSplus
# effect in per-SD terms. Black comes to 4.3% and Aboriginal to 0.02%.
#
# WHAT REMAINS WORTH ASKING, and it is the better question anyway: is the
# contamination concentrated where the main estimate is identified? If the
# Asian failure sits in metros and cohorts that contribute little to non_HS vs
# higher, it is a localised nuisance and the limitation sentence can say so. If
# it sits exactly where the main result is identified, the same unmodelled
# correlation plausibly contaminates both, and that is a materially stronger
# caveat. Note that MSA fixed effects absorb the cross-sectional level, so
# whatever drives either profile is within-metro variation across cohorts.
#
# METHOD
# For a 2SLS estimate, beta = (z~'y~) / (z~'x~), where ~ denotes residualising
# on controls and fixed effects. Each individual contributes z~_i * y~_i to the
# numerator. Aggregating those contributions by MSA and by birth cohort shows
# where the instrument-outcome covariance actually comes from. Computed for
# both Asian and y_HSplus on an identical control set so the two are directly
# comparable, then the two concentration profiles are correlated.
#
# SELF-SUFFICIENT as of 2026-09-29. It previously read samp from
# after_sec19.RData, which carried pre-respecification objects whenever Section
# 19 had not just been re-run. It now rebuilds the sample from
# analysis_data.RData under the same seed and construction as Sections 19 and
# 20, so all three use an identical draw without depending on session state or
# checkpoint freshness. Nothing written to the pipeline.
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
           "exposure_missing", "Educ_Group", "Asian")
.miss <- setdiff(.need, names(all_natives_175))
if (length(.miss))
  stop("analysis_data.RData is stale or incomplete -- missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R first.")
rm(.need, .miss)

# Outcome built BEFORE the draw -- samp is a copy.
all_natives_175[, y_HSplus := as.integer(Educ_Group != "non_HS")]

set.seed(42)
samp <- droplevels(all_natives_175[sample(nrow(all_natives_175), 500000)])
samp[, birth_year := as.integer(as.character(YEAR)) - AGE]

cat(sprintf("Sample: %d observations, %d MSAs\n",
            nrow(samp), uniqueN(samp$MSA)))

# --- Residualise on a common control set ------------------------------------
# Black is omitted throughout. 16_adao_placebo.R omits it (it cannot control
# for an outcome that is itself a placebo), and holding the control set fixed
# across both outcomes is what makes the two concentration profiles comparable.
# This means the y_HSplus contributions here will not reproduce the main
# estimate exactly; the object of interest is WHERE the covariance sits, not
# its level.
CTRL <- "SEX + Marital + ChildAgeGroup + POVERTY + Mom_Educ + Dad_Educ +
         exposure_missing | YEAR + MSA"

resid_of <- function(v) {
  residuals(feols(as.formula(paste(v, "~", CTRL)), data = samp))
}

cat("Residualising...\n")
samp[, `:=`(
  z_t      = resid_of("avg_exposure_nodiv"),
  y_asian  = resid_of("Asian"),
  y_hsplus = resid_of("y_HSplus")
)]

# RAW instrument, not the residualised one (changed 2026-09-29).
#
# The 2SLS numerator is z'y-perp, and because the fitted part of z lies in the
# control space while y-perp is orthogonal to it,
#
#     sum_i z_i * y-perp_i  ==  sum_i z-tilde_i * y-perp_i
#
# so both give the same TOTAL. They attribute it differently, though, and the
# residualised version attributes badly here. For anyone with
# exposure_missing == 1 the raw instrument is exactly zero, but z-tilde is
# 0 - fitted, which is not zero -- so those individuals picked up nonzero
# contributions despite carrying no exposure variation at all. That is why the
# earlier cohort table showed substantial Asian contributions from the
# 1900-1974 cohorts (+0.0161 for 1947, +0.0180 for 1962), every one of which
# has z = 0 by construction through the exposed_post_2000 condition.
#
# Using the raw instrument assigns exactly zero to non-identifying individuals
# and sums to the same numerator, so the resulting profiles describe where the
# covariance actually comes from. The outcomes stay residualised -- that part
# is required for the decomposition to be of the estimate rather than of a raw
# correlation.
samp[, `:=`(
  contrib_asian  = avg_exposure_nodiv * y_asian,
  contrib_hsplus = avg_exposure_nodiv * y_hsplus
)]
stopifnot(samp[exposure_missing == 1L, all(contrib_asian == 0 & contrib_hsplus == 0)])

# --- Concentration by MSA ----------------------------------------------------
by_msa <- samp[, .(
  n              = .N,
  contrib_asian  = sum(contrib_asian),
  contrib_hsplus = sum(contrib_hsplus)
), by = MSA]

tot_a <- sum(abs(by_msa$contrib_asian))
tot_h <- sum(abs(by_msa$contrib_hsplus))
by_msa[, `:=`(
  share_asian  = contrib_asian  / tot_a,
  share_hsplus = contrib_hsplus / tot_h
)]

# order() inside [ rather than setorder(), which takes column names only and
# will not accept an expression like -abs(share_asian).
by_msa <- by_msa[order(-abs(share_asian))]

cat("\n--- Top 15 MSAs by contribution to the ASIAN placebo covariance ---\n")
cat("(35620 NY, 31080 LA, 16980 Chicago, 41860 SF, 41940 San Jose,\n")
cat(" 46520 Honolulu, 42660 Seattle, 47900 DC, 33100 Miami, 26420 Houston)\n\n")
cat(sprintf("%-8s  %-8s  %-13s  %-13s\n", "MSA", "n", "share_asian", "share_hsplus"))
cat(strrep("-", 50), "\n")
for (i in seq_len(min(15, nrow(by_msa)))) {
  cat(sprintf("%-8s  %-8d  %+12.4f  %+12.4f\n",
              by_msa$MSA[i], by_msa$n[i],
              by_msa$share_asian[i], by_msa$share_hsplus[i]))
}

cat(sprintf("\nTop 1 MSA:  %.1f%% of |Asian| covariance, %.1f%% of |y_HSplus|\n",
            100 * abs(by_msa$share_asian[1]),
            100 * abs(by_msa$share_hsplus[1])))
cat(sprintf("Top 5 MSAs: %.1f%% of |Asian| covariance, %.1f%% of |y_HSplus|\n",
            100 * sum(abs(by_msa$share_asian[1:5])),
            100 * sum(abs(by_msa$share_hsplus[1:5]))))

cat(sprintf("\nCorrelation of MSA contribution profiles: %.4f\n",
            cor(by_msa$share_asian, by_msa$share_hsplus)))
cat(sprintf("Rank correlation:                        %.4f\n",
            cor(by_msa$share_asian, by_msa$share_hsplus, method = "spearman")))

# --- Concentration by birth cohort ------------------------------------------
by_coh <- samp[, .(
  n              = .N,
  contrib_asian  = sum(contrib_asian),
  contrib_hsplus = sum(contrib_hsplus)
), by = birth_year][order(birth_year)]

by_coh[, `:=`(
  share_asian  = contrib_asian  / sum(abs(contrib_asian)),
  share_hsplus = contrib_hsplus / sum(abs(contrib_hsplus))
)]

cat("\n--- Contribution by birth cohort ---\n")
cat(sprintf("%-8s  %-8s  %-13s  %-13s\n",
            "cohort", "n", "share_asian", "share_hsplus"))
cat(strrep("-", 50), "\n")
for (i in seq_len(nrow(by_coh))) {
  cat(sprintf("%-8d  %-8d  %+12.4f  %+12.4f\n",
              by_coh$birth_year[i], by_coh$n[i],
              by_coh$share_asian[i], by_coh$share_hsplus[i]))
}

cat("\nREADING THIS\n\n")
cat("If the Asian covariance is concentrated in a handful of metros that\n")
cat("contribute little to y_HSplus -- and the profile correlation is low --\n")
cat("the placebo failure is localised and can be reported as such. If the two\n")
cat("profiles track each other, the same unmodelled correlation plausibly\n")
cat("affects both, and the placebo failure is a warning about the main\n")
cat("estimate rather than a curiosity beside it.\n\n")
cat("Compare the Pearson and rank correlations. A high Pearson with a low\n")
cat("rank correlation means the two profiles agree because the same one or two\n")
cat("large metros dominate both, while the ordering elsewhere is unrelated.\n")
cat("For this question the mass is what matters, so the Pearson figure is the\n")
cat("relevant one.\n")

cat("\n", strrep("-", 60), "\n", sep = "")
cat("NOTE ON THE YOUNGEST IDENTIFYING COHORTS\n\n")
cat("For the write-up, not for a further sample restriction. The exposure\n")
cat("window is capped at 2022, so birth years 2005, 2006 and 2007 have windows\n")
cat("of three, two and one covered years -- someone born in 2007 has their\n")
cat("entire formative-years exposure measured from 2022 alone, at age fifteen.\n\n")
cat("On the raw-instrument attribution their contribution is NEGATIVE and\n")
cat("roughly proportional to their share of the sample, not outsized. Indeed\n")
cat("the whole 1998-2007 stretch contributes negatively to the y_HSplus\n")
cat("covariance while 1975-1997 contributes positively, so the censored\n")
cat("youngest cohorts pull the estimate DOWN. Dropping them would raise it.\n")
cat("That is the opposite of what the residualised attribution suggested and\n")
cat("the opposite of the pre-fix 18b gradient, so 18b must be re-run before\n")
cat("any 'restrict to older cohorts and the estimate shrinks' claim is made.\n\n")
cat("This is the same mechanism as the 2008-and-later cohorts removed from\n")
cat("07_exposure_iv.R on 2026-09-29, differing only in falling inside the\n")
cat("age_15_year <= 2022 filter rather than outside it. It is left in\n")
cat("deliberately. The age range is already narrowed, and documenting the\n")
cat("mechanism is what this script does in place of filtering on it.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(by_msa, "outputs/tables/asian_placebo_by_msa.csv", row.names = FALSE)
write.csv(by_coh, "outputs/tables/asian_placebo_by_cohort.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/asian_placebo_by_msa.csv, _by_cohort.csv\n")
cat("17_asian_placebo_concentration.R complete.\n")
