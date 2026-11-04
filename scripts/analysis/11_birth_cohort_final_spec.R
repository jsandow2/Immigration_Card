# =============================================================================
# 11_birth_cohort_final_spec.R  --  Birth-cohort check on the final
# Renumbered 2026-10-21 from 18b_birth_cohort_final_spec.R (archive-and-renumber pass).
#                                    specification
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Supersedes 18_birth_cohort_check.R, which runs on the survey-year endogenous
# variable, the AGE >= 30 floor, and the wave-divided series -- all superseded.
#
# THE JOB THIS SCRIPT HAS TO DO
# The headline specification (levels, windowed, undivided pooled counts, no base
# population, no age restriction, no age fixed effects) returns a reduced-form
# effect of roughly 8.9 percentage points per one-SD increase in exposure. That
# is implausibly large for an immigration effect on native attainment, and the
# chapter says so. This script is the diagnostic the text then points to, and
# the hand-off to the chapter that takes up age structure.
#
# WHY COHORT IS THE RIGHT DIAGNOSTIC
# Two facts line up. Later exposure windows carry lower z_lev. And later
# windows belong to younger people, who are far more likely to be coded non_HS
# because their schooling is not finished -- 30.4% of this sample against 7.8%
# under the old age floor. Low z_lev therefore travels with low y_HSplus for
# reasons unrelated to immigration, and without age controls the regression
# reads that as effect.
#
# If that is what inflates the pooled estimate, it has a signature: the
# estimate should FALL as the sample is restricted to older cohorts, whose
# schooling is complete and who therefore cannot contribute censoring.
#
# WHAT WOULD BREAK THE NARRATIVE
# A roughly FLAT profile across cohorts. That would mean cohort composition is
# not what drives the magnitude, the pooled estimate's size is unexplained, and
# the cross-reference forward would be pointing at a diagnosis this check does
# not support. That outcome should be reported rather than worked around.
#
# STRUCTURE
#   Panel A  Five-year birth-cohort bands, each estimated separately.
#   Panel B  Cumulative restriction -- birth_year <= C for a descending set of
#            C. Because restricting to earlier birth years means older at
#            survey, this is a de facto age restriction and produces the
#            "narrow the sample and the estimate drops" evidence directly.
#
# Every row reports the non_HS rate and mean age alongside the estimate, so the
# censoring is visible rather than inferred.
#
# Cohorts born before 1975 have avg_exposure = 0 by construction (the
# exposed_post_2000 condition in 07_exposure_iv.R requires the ages-15-25
# window to open in 1990 or later, the share year). They carry no exposure
# variation and cannot be estimated; the script confirms this rather than
# assuming it.
#
# Writes outputs/tables/birth_cohort_final_spec.csv.
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

# --- Final-spec inflow series: undivided pooled counts ----------------------
waves <- sort(unique(reduced_data_iv$YEAR))
pred  <- as.data.table(msa_year_predicted_inflows)
nw <- data.table(YEAR = sort(unique(pred$YEAR)))
nw[, n_waves := vapply(YEAR, function(y) sum(waves >= y), integer(1))]
pred <- merge(pred, nw, by = "YEAR", all.x = TRUE)
pred[, I_pred := total_predicted_inflow * n_waves]

cat("Recovering actual_inflows_msa_year from after_sec06.RData...\n")
.e <- new.env()
load(file.path(local_checkpoint_dir, "after_sec06.RData"), envir = .e)
act <- as.data.table(get("actual_inflows_msa_year", envir = .e))
rm(.e); gc()
act <- merge(act, nw, by = "YEAR", all.x = TRUE)
act[, I_act := total_actual_inflow * n_waves]

setnames(pred, "YEAR", "window_year")
setnames(act,  "YEAR", "window_year")

# --- Sample and windowed exposure -------------------------------------------
d <- reduced_data_iv[immigr == 0 & MSA != 0 & MSA %in% selected_msas]
d[, birth_year := YEAR - AGE]

elig <- d[exposure_missing == 0, .(PERSON_ID, MSA, AGE, YEAR)]
elig[, birth_year := YEAR - AGE]
elig[, `:=`(start_year = birth_year + 15L,
            end_year   = pmin(birth_year + 25L, 2022L))]
elig <- elig[start_year <= end_year]
expanded <- elig[, .(window_year = seq(start_year, end_year)),
                 by = .(PERSON_ID, MSA)]

w <- merge(expanded, pred[, .(MSA, window_year, I_pred)],
           by = c("MSA", "window_year"), all.x = TRUE)
w <- merge(w, act[, .(MSA, window_year, I_act)],
           by = c("MSA", "window_year"), all.x = TRUE)
agg <- w[, .(z_lev = mean(I_pred, na.rm = TRUE),
             x_lev = mean(I_act,  na.rm = TRUE)), by = PERSON_ID]

d <- merge(d, agg, by = "PERSON_ID", all.x = TRUE)
d[is.na(z_lev), z_lev := 0]
d[is.na(x_lev), x_lev := 0]
d[, y_HSplus := as.integer(Educ_Group != "non_HS")]
d <- droplevels(d)

sd_z <- d[exposure_missing == 0, sd(z_lev)]

cat(sprintf("\nSample: %s obs | %d MSAs | %s identifying\n",
            format(nrow(d), big.mark = ","), uniqueN(d$MSA),
            format(d[exposure_missing == 0, .N], big.mark = ",")))
cat(sprintf("non_HS overall: %.1f%% | SD of z_lev among identifying: %.1f\n",
            100 * d[, mean(Educ_Group == "non_HS")], sd_z))

# --- Confirm the pre-1975 cohorts carry no exposure variation ---------------
cat("\nExposure variation by cohort era (the exposed_post_2000 condition):\n")
print(d[, .(n = .N, mean_z = mean(z_lev), sd_z = sd(z_lev),
            pct_exposure_missing = 100 * mean(exposure_missing == 1)),
        by = .(era = fifelse(birth_year < 1975, "born < 1975", "born >= 1975"))],
      row.names = FALSE)
cat("Bands below are estimated on born >= 1975 only -- earlier cohorts have no\n")
cat("within-cohort exposure variation and cannot identify anything.\n")

# --- Machinery ---------------------------------------------------------------
CTRL <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
         Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

fit_band <- function(dd, label, panel) {
  if (nrow(dd) < 5000 || dd[, uniqueN(MSA)] < 20 ||
      dd[exposure_missing == 0, .N] < 1000) {
    cat(sprintf("  %-22s SKIPPED (too thin: n %s, MSAs %d)\n",
                label, format(nrow(dd), big.mark = ","), dd[, uniqueN(MSA)]))
    return(NULL)
  }
  dd <- droplevels(dd)
  ok <- tryCatch({
    rf <- feols(as.formula(sprintf("y_HSplus ~ z_lev + %s", CTRL)),
                data = dd, vcov = ~MSA)
    fs <- feols(as.formula(sprintf("x_lev ~ z_lev + %s", CTRL)),
                data = dd, vcov = ~MSA)
    iv <- feols(as.formula(sprintf("y_HSplus ~ %s | x_lev ~ z_lev", CTRL)),
                data = dd, vcov = ~MSA)
    TRUE
  }, error = function(e) { cat(sprintf("  %-22s FAILED: %s\n", label,
                                       conditionMessage(e))); FALSE })
  if (!ok) return(NULL)

  ct <- iv$coeftable["fit_x_lev", ]
  b_rf <- unname(coef(rf)["z_lev"])

  cat(sprintf("  %-22s RF %+.4e  p %.4f | %+7.3f pp/SD | non_HS %5.1f%% | age %4.1f | n %9s\n",
              label, b_rf, rf$coeftable["z_lev", "Pr(>|t|)"], 100 * b_rf * sd_z,
              100 * dd[, mean(Educ_Group == "non_HS")], dd[, mean(AGE)],
              format(nrow(dd), big.mark = ",")))

  data.frame(panel = panel, band = label, n = nrow(dd),
             n_identifying = dd[exposure_missing == 0, .N],
             mean_age = dd[, mean(AGE)],
             pct_non_HS = 100 * dd[, mean(Educ_Group == "non_HS")],
             first_stage = unname(coef(fs)["z_lev"]),
             reduced_form = b_rf,
             rf_se = rf$coeftable["z_lev", "Std. Error"],
             rf_p = rf$coeftable["z_lev", "Pr(>|t|)"],
             pp_per_sd = 100 * b_rf * sd_z,
             tsls = ct[["Estimate"]], tsls_se = ct[["Std. Error"]],
             stringsAsFactors = FALSE)
}

res <- list()

# --- Reference: the headline -------------------------------------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("REFERENCE: pooled headline specification\n")
res[[length(res)+1]] <- fit_band(d, "all cohorts (headline)", "reference")

# --- Panel A: five-year cohort bands ----------------------------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PANEL A: five-year birth-cohort bands (born >= 1975)\n")
cat("Expect the estimate to FALL toward the older bands if censoring drives it.\n\n")
brk <- c(1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010)
for (i in seq_len(length(brk) - 1L)) {
  lo <- brk[i]; hi <- brk[i+1] - 1L
  res[[length(res)+1]] <- fit_band(d[birth_year >= lo & birth_year <= hi],
                                   sprintf("born %d-%d", lo, hi), "A: bands")
}

# --- Panel B: cumulative restriction to older cohorts -----------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PANEL B: cumulative restriction (born >= 1975 and <= C)\n")
cat("Dropping the youngest cohorts is a de facto age restriction.\n\n")
for (C in c(2009, 2000, 1995, 1990, 1985)) {
  res[[length(res)+1]] <- fit_band(d[birth_year >= 1975 & birth_year <= C],
                                   sprintf("born 1975-%d", C), "B: cumulative")
}

res <- rbindlist(Filter(Negate(is.null), res))

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SUMMARY\n\n")
print(res[, .(panel, band, n, mean_age, pct_non_HS, reduced_form,
              rf_se, rf_p, pp_per_sd)], row.names = FALSE, digits = 3)

cat("\n", strrep("-", 72), "\n", sep = "")
cat("READING THIS\n\n")
cat("Read pp_per_sd against pct_non_HS and mean_age, row by row. The censoring\n")
cat("story predicts they move together: bands with high non_HS rates and low\n")
cat("mean age carry the large estimates, and the estimate declines as the\n")
cat("sample shifts to older cohorts whose schooling is complete.\n\n")
cat("If Panel B declines toward something near the age-FE magnitude (roughly\n")
cat("+0.3 pp per SD), that is the diagnosis: the pooled 8.9 pp is cohort\n")
cat("censoring, the chapter can say so with a number attached, and the\n")
cat("cross-reference forward is earned.\n\n")
cat("If Panel B is FLAT, the magnitude is not a cohort artefact. The pooled\n")
cat("estimate would then be unexplained rather than provisionally explained,\n")
cat("and the text should not attribute it to age structure on the strength of\n")
cat("a check that failed to find it.\n\n")
cat("Note that narrow bands have less within-band exposure variation, so\n")
cat("Panel A standard errors widen mechanically toward the edges. Panel B is\n")
cat("the cleaner read because each row retains the full older sample.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/birth_cohort_final_spec.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/birth_cohort_final_spec.csv\n")
cat("11_birth_cohort_final_spec.R complete.\n")
