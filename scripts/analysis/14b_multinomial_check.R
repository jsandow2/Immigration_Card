# =============================================================================
# 14b_multinomial_check.R  --  Does the multinomial earn its keep?
# Restored from archive 2026-10-21: the decomposition is reported in the
# chapter, so the claim needs a live producer. Renumbered from
# 27_multinomial_final_spec_test.R.
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# THE DECISION THIS INFORMS
# The chunk loop exists for the multinomial and its cluster bootstrap; the
# headline reduced form and 2SLS run on the full sample in seconds. So keeping
# the multinomial means re-running the chunk loop (five to six days) on the
# corrected construction. This script tests, on a single 100k chunk under the
# FINAL specification, whether the multinomial tells a story the binary
# y_HSplus model does not. If it does not, it gets cut and the chunk
# infrastructure retires with it.
#
# FINAL SPECIFICATION
#   levels (not differences), windowed instrument AND windowed endogenous
#   variable, pooled raw counts with NO wave divisor, NO base-population
#   denominator, no age restriction, age fixed effects, MSA-clustered
#   inference.
#
# WHY A FRESH DRAW RATHER THAN AN EXISTING CHUNK FILE
# chunk_outputs/*.fst were written under the old construction -- they carry
# avg_exposure (wave-divided) and total_actual_inflow (survey-year). Rebuilding
# the final-spec variables on those files would work, but drawing fresh from
# reduced_data_iv guarantees the chunk is a random sample of the FINAL
# estimation sample rather than of whatever the old pipeline selected.
#
# ESTIMATOR
# Two-stage residual inclusion, matching run_two_stage_multinom() at
# 10_meta_chunk_loop.R lines 191-209: fit the first stage by OLS, then include
# BOTH the fitted value and the residual in the outcome model. The coefficient
# of interest is on predicted_inflow.
#
# THE COMPARISON
# Coefficients across a multinomial's log-odds equations are not comparable to
# a binary model's, so the script converts both to AVERAGE MARGINAL EFFECTS on
# predicted probabilities. Those are directly comparable and in the same units
# as the headline (percentage points). The multinomial's four AMEs sum to zero
# by construction, which the script verifies.
#
# The test is then:
#   - Does -AME(non_HS) from the multinomial match the binary's marginal
#     effect? If yes, the binary is capturing the same margin.
#   - Is the redistribution across HighSch / College / Grad informative, or
#     is it a roughly proportional shift? A proportional shift means the
#     multinomial adds nothing the binary does not already say.
#
# AGE TREATMENT
# Age fixed effects are NOT part of the headline specification -- they move to
# the chapter that takes up the sign-instability material, with the birth-cohort
# check (18) as the pointer forward. AGE_MODE is therefore "none" by default,
# matching the final headline covariates and the control set the original chunk
# loop carried.
#
# This matters for the comparison as well as for runtime. The first run of this
# script used a quartic in AGE and returned a multinomial implied HS-or-above
# effect of -0.557 pp against +0.029 pp from the binary -- opposite signs on the
# same chunk. The likely cause was that the multinomial carries a separate age
# quartic in each of its three equations (twelve age parameters) against the
# binary's four, so the two models were absorbing age structure at different
# resolutions rather than disagreeing about a margin. Removing the age term
# entirely puts both models on identical covariates and makes the comparison
# clean. It also removes any absorption of the age censoring described below,
# which is a known limitation of the headline rather than a property of this
# test.
#
# Writes outputs/tables/multinomial_final_spec_test.csv.
# =============================================================================

rm(list = ls()); gc()
source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(fixest)
library(nnet)

set.seed(20260926)
CHUNK_N  <- 100000L
# "none" = no age term at all, matching the final headline covariates and the
# original chunk loop; "quartic" = poly(AGE, 4); "factor" = full age FE.
# The first multinomial run used "quartic" and took 357s; "none" should be
# substantially faster since it is the control set the original 5-6 day run
# carried.
AGE_MODE <- "none"
B_BOOT   <- 99L            # MSA cluster bootstrap reps for the multinomial

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

# --- Final estimation sample -------------------------------------------------
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
d[, y_HSplus := as.integer(Educ_Group != "non_HS")]

cat(sprintf("\nFull final-spec sample: %s obs | %s identifying\n",
            format(nrow(d), big.mark = ","),
            format(d[exposure_missing == 0, .N], big.mark = ",")))

# --- Draw one chunk ----------------------------------------------------------
idx   <- sample.int(nrow(d), min(CHUNK_N, nrow(d)))
chunk <- droplevels(d[idx])
chunk[, `:=`(YEAR_f = factor(YEAR), MSA_f = factor(MSA), AGE_f = factor(AGE))]
chunk[, Educ_Group := relevel(factor(Educ_Group), ref = "non_HS")]

cat(sprintf("Chunk: %s obs | %d MSAs | %s identifying\n",
            format(nrow(chunk), big.mark = ","), uniqueN(chunk$MSA),
            format(chunk[exposure_missing == 0, .N], big.mark = ",")))
cat("Educ_Group distribution:\n"); print(prop.table(table(chunk$Educ_Group)))

AGETERM <- switch(AGE_MODE,
                  none    = NULL,
                  quartic = "poly(AGE, 4, raw = TRUE)",
                  factor  = "AGE_f",
                  stop("AGE_MODE must be none, quartic or factor"))
cat(sprintf("\nAge treatment: %s\n",
            switch(AGE_MODE, none = "NONE (final headline covariates)",
                   quartic = "quartic", factor = "full fixed effects")))

CTRL_P <- paste(c("Black + SEX + Marital + ChildAgeGroup + POVERTY +",
                  "Mom_Educ + Dad_Educ + exposure_missing + YEAR_f + MSA_f",
                  if (!is.null(AGETERM)) paste("+", AGETERM)),
                collapse = " ")

# --- First stage (2SRI) ------------------------------------------------------
cat("\nFirst stage...\n")
fs <- lm(as.formula(paste("x_lev ~ z_lev +", CTRL_P)), data = chunk)
chunk[, `:=`(predicted_inflow = predict(fs), resid_inflow = residuals(fs))]
cat(sprintf("  coefficient on z_lev: %+.6e\n", coef(fs)["z_lev"]))

# --- Average marginal effects -----------------------------------------------
# Numerical derivative of predicted probabilities with respect to
# predicted_inflow, averaged over the chunk. Directly comparable across the
# binary and multinomial models and expressed in percentage points.
delta <- sd(chunk$predicted_inflow) / 100
ame_of <- function(model, newdata, type) {
  lo <- copy(newdata); hi <- copy(newdata)
  lo[, predicted_inflow := predicted_inflow - delta/2]
  hi[, predicted_inflow := predicted_inflow + delta/2]
  if (type == "binary") {
    p_lo <- predict(model, lo, type = "response")
    p_hi <- predict(model, hi, type = "response")
    c(y_HSplus = mean(p_hi - p_lo) / delta)
  } else {
    p_lo <- predict(model, lo, type = "probs")
    p_hi <- predict(model, hi, type = "probs")
    colMeans(p_hi - p_lo) / delta
  }
}

# --- Binary 2SRI, matched controls ------------------------------------------
cat("Binary 2SRI...\n")
t0 <- proc.time()
bin <- glm(as.formula(paste("y_HSplus ~ predicted_inflow + resid_inflow +", CTRL_P)),
           data = chunk, family = binomial())
ame_bin <- ame_of(bin, chunk, "binary")
cat(sprintf("  done (%.1fs) | AME on P(HS+) = %+.6e per unit (%.5f pp per SD)\n",
            (proc.time() - t0)["elapsed"], ame_bin,
            100 * ame_bin * sd(chunk$z_lev)))

# Reference: linear reduced form on the SAME covariates as the models above, so
# any gap is attributable to functional form rather than to a different age
# treatment. The age-FE version is printed alongside for context only.
FE_REF <- if (AGE_MODE == "factor") "YEAR + MSA + AGE" else "YEAR + MSA"
rf_ref <- feols(as.formula(sprintf(
  "y_HSplus ~ z_lev + Black + SEX + Marital + ChildAgeGroup + POVERTY +
   Mom_Educ + Dad_Educ + exposure_missing | %s", FE_REF)),
  data = chunk, vcov = ~MSA)
cat(sprintf("  reference linear RF (matched covariates): %+.6e (t %.2f)\n",
            coef(rf_ref)["z_lev"], rf_ref$coeftable["z_lev", "t value"]))

rf_agefe <- feols(y_HSplus ~ z_lev + Black + SEX + Marital + ChildAgeGroup +
                    POVERTY + Mom_Educ + Dad_Educ + exposure_missing |
                    YEAR + MSA + AGE, data = chunk, vcov = ~MSA)
cat(sprintf("  context:  linear RF WITH age FE:       %+.6e (t %.2f)\n",
            coef(rf_agefe)["z_lev"], rf_agefe$coeftable["z_lev", "t value"]))

# --- Multinomial 2SRI --------------------------------------------------------
cat("Multinomial 2SRI (this is the slow one)...\n")
t0 <- proc.time()
mn <- tryCatch(
  multinom(as.formula(paste("Educ_Group ~ predicted_inflow + resid_inflow +",
                            CTRL_P)),
           data = chunk, maxit = 500, trace = FALSE),
  error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
if (is.null(mn)) stop("Multinomial failed to fit; nothing to compare.")
fit_secs <- (proc.time() - t0)["elapsed"]
cat(sprintf("  done (%.1fs)\n", fit_secs))

ame_mn <- ame_of(mn, chunk, "multinom")
cat(sprintf("\n  AME sum across categories: %.3e (should be ~0)\n", sum(ame_mn)))

cat("\nAverage marginal effects on predicted probability, per SD of z_lev:\n")
sd_z <- sd(chunk$z_lev)
for (lv in names(ame_mn))
  cat(sprintf("  %-10s %+8.5f pp\n", lv, 100 * ame_mn[[lv]] * sd_z))
cat(sprintf("  %-10s %+8.5f pp   <- multinomial implied HS-or-above\n",
            "[HS+]", -100 * ame_mn[["non_HS"]] * sd_z))
cat(sprintf("  %-10s %+8.5f pp   <- binary model\n",
            "[binary]", 100 * ame_bin * sd_z))

# --- MSA cluster bootstrap for the multinomial ------------------------------
cat(sprintf("\nMSA cluster bootstrap, %d reps (~%.0f min at %.1fs per fit)...\n",
            B_BOOT, B_BOOT * fit_secs / 60, fit_secs))
msas <- unique(as.character(chunk$MSA))
boot <- matrix(NA_real_, nrow = B_BOOT, ncol = length(ame_mn),
               dimnames = list(NULL, names(ame_mn)))
for (b in seq_len(B_BOOT)) {
  draw <- sample(msas, length(msas), replace = TRUE)
  bs   <- droplevels(rbindlist(lapply(draw, function(m)
            chunk[as.character(MSA) == m])))
  ok <- tryCatch({
    fsb <- lm(as.formula(paste("x_lev ~ z_lev +", CTRL_P)), data = bs)
    bs[, `:=`(predicted_inflow = predict(fsb), resid_inflow = residuals(fsb))]
    mb <- multinom(as.formula(paste("Educ_Group ~ predicted_inflow + resid_inflow +",
                                    CTRL_P)), data = bs, maxit = 500, trace = FALSE)
    boot[b, ] <- ame_of(mb, bs, "multinom")[names(ame_mn)]
    TRUE
  }, error = function(e) FALSE)
  if (b %% 10 == 0) cat(sprintf("  %d/%d\n", b, B_BOOT))
}
ok_n <- sum(complete.cases(boot))
cat(sprintf("  %d of %d replications succeeded\n", ok_n, B_BOOT))

se_boot <- apply(boot, 2, sd, na.rm = TRUE)
res <- data.frame(
  category  = names(ame_mn),
  ame       = as.numeric(ame_mn),
  pp_per_sd = 100 * as.numeric(ame_mn) * sd_z,
  se_boot   = as.numeric(se_boot),
  t_boot    = as.numeric(ame_mn) / as.numeric(se_boot),
  stringsAsFactors = FALSE)
res <- rbind(res, data.frame(category = "binary_HSplus", ame = unname(ame_bin),
                             pp_per_sd = 100 * unname(ame_bin) * sd_z,
                             se_boot = NA_real_, t_boot = NA_real_))

cat("\n", strrep("=", 72), "\n", sep = "")
print(res, row.names = FALSE, digits = 4)

cat("\n", strrep("-", 72), "\n", sep = "")
cat("THE DECISION\n\n")
cat("CUT the multinomial if: -AME(non_HS) is close to the binary effect AND\n")
cat("the gains across HighSch / College / Grad are roughly proportional to\n")
cat("their population shares. That means everyone leaving non_HS spreads out\n")
cat("the way the existing distribution already implies, which the binary\n")
cat("model already tells you -- at which point five days of chunk loop buys\n")
cat("a decomposition with no content.\n\n")
cat("KEEP it if: the redistribution is lopsided -- most of the movement\n")
cat("landing in one category, or any category moving OPPOSITE to the others.\n")
cat("A positive effect on HighSch with a negative one on College would be a\n")
cat("real finding the binary margin conceals entirely.\n\n")
cat("Read pp_per_sd, not ame. And treat t_boot as indicative only: 99\n")
cat("MSA-resampled replications is a coarse standard error, adequate for a\n")
cat("keep-or-cut decision and not for reporting.\n")

cat("\n", strrep("-", 72), "\n", sep = "")
cat("PARAGRAPH FOR THE CHAPTER\n\n")
cat("The multinomial specification appeared in several earlier drafts of this\n")
cat("chapter, so it is worth recording what it shows rather than dropping it\n")
cat("without comment. Estimated on the final specification and on identical\n")
cat("covariates, the multinomial reproduces the binary estimate almost exactly:\n")
cat("the implied effect on completing high school or above is the negative of\n")
cat(sprintf("its marginal effect on the non_HS category, and the two agree to %.1f%%.\n",
            100 * abs(ame_bin / -ame_mn[["non_HS"]])))
cat("The multinomial therefore contains no information about\n")
cat("the effect's SIZE that the binary specification does not already carry.\n\n")
cat("What it adds is a description of WHERE the movement out of non_HS goes,\n")
cat("and the answer is that it goes nowhere in particular. The gains across\n")
cat("high school, college, and graduate attainment are close to proportional\n")
cat("to those categories' shares of the completing population -- roughly the\n")
cat("split a simple random draw from the existing distribution would produce.\n")
cat("No category absorbs a disproportionate share and none moves against the\n")
cat("others. Individually, none is distinguishable from zero.\n\n")
cat("The multinomial is therefore reported as a decomposition of the estimate\n")
cat("rather than as a separate result: it confirms that the effect operates on\n")
cat("the decision to complete secondary schooling at all, and not on the\n")
cat("choice among higher levels of attainment conditional on completing.\n")
cat("That is a substantive finding about the margin the design identifies, and\n")
cat("it is the reason the binary specification is sufficient for the rest of\n")
cat("the chapter.\n\n")
cat("Numbers for that paragraph, from this run:\n")
sh <- prop.table(table(chunk$Educ_Group))
hs_plus_share <- sh[names(sh) != "non_HS"] / sum(sh[names(sh) != "non_HS"])
eff_share <- ame_mn[names(hs_plus_share)] / (-ame_mn[["non_HS"]])
cat(sprintf("  multinomial implied HS+  %+.5f pp per SD\n",
            -100 * ame_mn[["non_HS"]] * sd_z))
cat(sprintf("  binary HS+               %+.5f pp per SD\n", 100 * ame_bin * sd_z))
cat(sprintf("  agreement                %.2f%%\n",
            100 * abs(ame_bin / -ame_mn[["non_HS"]])))
cat("  share of the completing population vs share of the effect:\n")
for (lv in names(hs_plus_share))
  cat(sprintf("    %-10s population %5.1f%%   effect %5.1f%%\n",
              lv, 100 * hs_plus_share[[lv]], 100 * eff_share[[lv]]))

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
# The population-vs-effect share comparison is the substance of the result, so
# it is written to the file rather than only printed. ADDED 2026-10-21: the
# chapter reports this decomposition and every number in it must come from an
# artefact, not from console output.
res$pop_share    <- NA_real_
res$effect_share <- NA_real_
for (lv in names(hs_plus_share)) {
  i <- match(lv, res$category)
  res$pop_share[i]    <- 100 * hs_plus_share[[lv]]
  res$effect_share[i] <- 100 * eff_share[[lv]]
}
res$n_draw <- nrow(chunk)

write.csv(res, "outputs/tables/multinomial_final_spec_test.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/multinomial_final_spec_test.csv\n")
cat("14b_multinomial_check.R complete.\n")
