# =============================================================================
# 12_final_spec_no_age_fe.R  --  Final specification, with and without age
# Renumbered 2026-10-21 from 28_final_spec_no_age_fe.R (archive-and-renumber pass).
#                                fixed effects
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# WHERE THIS SITS
# The chapter's headline is the LEVELS specification: windowed instrument and
# windowed endogenous variable, pooled raw counts with no wave divisor, no
# base-population denominator, no age restriction. The sign-instability
# material (24, 25, 26) moves to its own chapter, and the birth-cohort check
# (18) becomes the pointer to it.
#
# Age fixed effects are therefore dropped from the headline and become part of
# what the later chapter takes up. This script estimates the specification both
# ways so the difference is documented rather than assumed, and so the
# cross-reference in the text rests on a number.
#
# WHAT IS RUN
#   {100k chunk, full sample} x {no age FE, age FE} x {y_HSplus, y_CollPlus}
#
# The chunk is included because that is how the estimation was requested, but
# linear models run on the full 8.5M in seconds, so both are reported. If they
# disagree materially the chunk is unrepresentative and that is worth knowing;
# if they agree, the full-sample numbers are the ones to report.
#
# A CAVEAT THAT BELONGS IN THE TEXT
# Without the age restriction, non_HS is roughly 30% of the sample against 7.8%
# under the old AGE >= 30 floor (27_multinomial_final_spec_test.R measured
# 30.4% on a 100k draw). Most of that excess is people young enough that
# schooling is not yet complete. Age fixed effects absorb the age main effect
# but not an interaction between age and exposure -- a 16-year-old coded non_HS
# cannot respond to exposure the way a 25-year-old can. Dropping age FE removes
# even the main-effect absorption, so the gap between the two columns here is a
# lower bound on how much age structure matters. Card's own answer (16-68 with
# at least one year of potential experience, to eliminate students) is the
# obvious remedy and is not implemented here.
#
# INFERENCE
# Restricted wild cluster bootstrap on every full-sample reduced form. CRVE
# p-values have run optimistic by four to six times throughout this project and
# are reported only as a benchmark. Effective clusters (Carter-Schnepel-
# Steigerwald G*) reported for each configuration.
#
# Writes outputs/tables/final_spec_no_age_fe.csv.
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

set.seed(20260928)
CHUNK_N <- 100000L
B_BOOT  <- 9999L

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
d[, `:=`(y_HSplus   = as.integer(Educ_Group != "non_HS"),
         y_CollPlus = as.integer(Educ_Group %in% c("College", "Grad")))]
d <- droplevels(d)

cat(sprintf("\nFull sample: %s obs | %d MSAs | %s identifying\n",
            format(nrow(d), big.mark = ","), uniqueN(d$MSA),
            format(d[exposure_missing == 0, .N], big.mark = ",")))
cat(sprintf("Outcome means: y_HSplus %.4f (non_HS %.1f%%) | y_CollPlus %.4f\n",
            d[, mean(y_HSplus)], 100 * d[, mean(Educ_Group == "non_HS")],
            d[, mean(y_CollPlus)]))
sd_z <- d[exposure_missing == 0, sd(z_lev)]
cat(sprintf("SD of z_lev among identifying observations: %.1f\n", sd_z))

chunk <- droplevels(d[sample.int(nrow(d), min(CHUNK_N, nrow(d)))])
cat(sprintf("Chunk: %s obs | %d MSAs | %s identifying\n",
            format(nrow(chunk), big.mark = ","), uniqueN(chunk$MSA),
            format(chunk[exposure_missing == 0, .N], big.mark = ",")))

# --- Machinery ---------------------------------------------------------------
BASE <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
         Mom_Educ + Dad_Educ + exposure_missing"

boot_p <- function(yt, zt, g, B, weights = "rademacher") {
  zz <- sum(zt^2); b <- sum(zt * yt) / zz
  s  <- rowsum(zt * (yt - b * zt), g, reorder = FALSE)[, 1]
  t_obs <- b / (sqrt(sum(s^2)) / zz)
  Sg <- rowsum(zt * yt, g, reorder = FALSE)[, 1]
  Q  <- rowsum(zt^2,    g, reorder = FALSE)[, 1]
  wv <- if (weights == "webb")
    c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2)) else c(-1, 1)
  V  <- matrix(sample(wv, length(Sg) * B, replace = TRUE), nrow = length(Sg))
  bs <- crossprod(V, Sg)[, 1] / zz
  Sc <- (V * Sg) - outer(Q, bs)
  ts <- bs / (sqrt(colSums(Sc^2)) / zz)
  list(t = t_obs, p = (1 + sum(abs(ts) >= abs(t_obs))) / (B + 1))
}

run_one <- function(dd, outcome, age_fe, label, do_boot) {
  fe   <- if (age_fe) "YEAR + MSA + AGE" else "YEAR + MSA"
  ctrl <- paste(BASE, "|", fe)

  fs <- feols(as.formula(sprintf("x_lev ~ z_lev + %s", ctrl)), data = dd, vcov = ~MSA)
  rf <- feols(as.formula(sprintf("%s ~ z_lev + %s", outcome, ctrl)), data = dd, vcov = ~MSA)
  iv <- feols(as.formula(sprintf("%s ~ %s | x_lev ~ z_lev", outcome, ctrl)),
              data = dd, vcov = ~MSA)
  ct <- iv$coeftable["fit_x_lev", ]

  zt <- residuals(feols(as.formula(paste("z_lev ~", ctrl)), data = dd))
  g  <- as.character(dd$MSA)
  gm <- tapply(zt^2, g, sum)
  cv2 <- mean((gm - mean(gm))^2) / mean(gm)^2
  Gstar <- length(gm) / (1 + cv2)

  bp <- bw <- NA_real_
  if (do_boot) {
    yt <- residuals(feols(as.formula(paste(outcome, "~", ctrl)), data = dd))
    bp <- boot_p(yt, zt, g, B_BOOT)$p
    bw <- boot_p(yt, zt, g, B_BOOT, "webb")$p
  }

  cat(sprintf("\n  %s\n", label))
  cat(sprintf("    first stage  %+.6e (t %6.2f)\n",
              coef(fs)["z_lev"], fs$coeftable["z_lev", "t value"]))
  cat(sprintf("    reduced form %+.6e (t %6.2f)  CRVE p %.4f%s\n",
              coef(rf)["z_lev"], rf$coeftable["z_lev", "t value"],
              rf$coeftable["z_lev", "Pr(>|t|)"],
              if (do_boot) sprintf("  boot p %.4f (Webb %.4f)", bp, bw) else ""))
  cat(sprintf("    2SLS         %+.6e SE %.4e  p %.4f\n",
              ct[["Estimate"]], ct[["Std. Error"]], ct[["Pr(>|t|)"]]))
  cat(sprintf("    %+.4f pp per identifying SD | G* %.1f of %d\n",
              100 * coef(rf)["z_lev"] * sd_z, Gstar, length(gm)))

  data.frame(spec = label, outcome = outcome, age_fe = age_fe,
             n = nobs(rf), first_stage = unname(coef(fs)["z_lev"]),
             reduced_form = unname(coef(rf)["z_lev"]),
             rf_se = rf$coeftable["z_lev", "Std. Error"],
             rf_p_crve = rf$coeftable["z_lev", "Pr(>|t|)"],
             rf_p_boot = bp, rf_p_boot_webb = bw,
             pp_per_sd = 100 * unname(coef(rf)["z_lev"]) * sd_z,
             tsls = ct[["Estimate"]], tsls_se = ct[["Std. Error"]],
             tsls_p = ct[["Pr(>|t|)"]], Gstar = Gstar,
             stringsAsFactors = FALSE)
}

res <- list()
for (oc in c("y_HSplus", "y_CollPlus")) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(sprintf("%s\n", oc))
  cat(strrep("-", 72), "\n", sep = "")
  cat("CHUNK (100k)\n")
  res[[length(res)+1]] <- run_one(chunk, oc, FALSE, "chunk | no age FE", FALSE)
  res[[length(res)+1]] <- run_one(chunk, oc, TRUE,  "chunk | age FE",    FALSE)
  cat("\n", strrep("-", 72), "\n", sep = "")
  cat("FULL SAMPLE\n")
  res[[length(res)+1]] <- run_one(d, oc, FALSE, "full | no age FE", TRUE)
  res[[length(res)+1]] <- run_one(d, oc, TRUE,  "full | age FE",    TRUE)
}
res <- rbindlist(res)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SUMMARY\n\n")
print(res[, .(spec, outcome, first_stage, reduced_form, rf_p_crve, rf_p_boot,
              pp_per_sd, Gstar)], row.names = FALSE, digits = 3)

cat("\n", strrep("-", 72), "\n", sep = "")
cat("READING THIS\n\n")
cat("Chunk against full sample first: if they disagree in sign or by more than\n")
cat("a factor of two, the chunk is unrepresentative and only the full-sample\n")
cat("rows should be used.\n\n")
cat("Then no-age-FE against age-FE. That gap is the quantity the chapter-three\n")
cat("cross-reference rests on, so it needs to be stated as a number rather\n")
cat("than as a direction. If the headline holds its sign and rough magnitude\n")
cat("across both, the age question is a refinement and the pointer forward is\n")
cat("honest. If the sign moves, the headline is conditional on excluding age\n")
cat("fixed effects and the text has to say so plainly -- deferring it to a\n")
cat("later chapter would not be enough.\n\n")
cat("Read rf_p_boot, not rf_p_crve.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/final_spec_no_age_fe.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/final_spec_no_age_fe.csv\n")
cat("12_final_spec_no_age_fe.R complete.\n")
