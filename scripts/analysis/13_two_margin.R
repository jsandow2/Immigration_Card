# =============================================================================
# 13_two_margin.R  --  Level and change of immigrant exposure, jointly
# Renumbered 2026-10-21 from 26_two_margin.R (archive-and-renumber pass).
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# WHY
# 25_level_difference_test.R found the two specifications give opposite signs,
# both at p ~ 0.001:
#
#   levels       2SLS +8.315e-07   (first stage +0.642)
#   differences  2SLS -5.559e-06   (first stage +0.563)
#
# with cor(level, difference) = -0.469. Part of that opposition is mechanical,
# but only part -- 78% of the difference's variance is independent of the
# level, so each carries information the other does not.
#
# Critically, this sign flip appears between two SIMPLE specifications: no wave
# divisor, no base-population denominator, no age restriction. The sign
# sensitivity is therefore a property of the underlying variation rather than
# an artefact of any normalisation choice, and no further simplification will
# resolve it.
#
# THE ARGUMENT FOR ESTIMATING BOTH
# They are not rival measures of one quantity. They are two margins of the
# environment a person grew up in:
#
#   LEVEL   how much immigration was present during ages 15-25 -- established
#           communities, networks, accumulated presence
#   CHANGE  how fast it was moving -- shock, disruption, adjustment
#
# A metro with high steady inflows and one with rapidly rising inflows are
# different places to be sixteen, and there is no reason their coefficients
# should share a sign. Entering only one forces the other's variation into the
# error term, where it is correlated with the included regressor at -0.469.
# That is textbook omitted-variable bias, and it is the most likely explanation
# for why the single-margin estimates disagree.
#
# NOTE ON THE DIFFERENCE MEASURE
# 25 verified numerically (max abs diff 2.3e-13) that averaging first
# differences across the window telescopes to [I(end) - I(start-1)] / n_window.
# So CHANGE means NET change across adolescence, not the path through it. Two
# people whose windows ended where they began score identically regardless of
# what happened in between. Describe it that way in the text.
#
# WHAT IS ESTIMATED
#   1. Reduced form, each margin alone (the 25 results, for reference)
#   2. Reduced form, both margins jointly
#   3. 2SLS, two endogenous variables instrumented by their two predicted
#      counterparts -- just-identified, no over-identification machinery
#
# Inference: cluster-robust at MSA plus a restricted wild cluster bootstrap on
# each reduced-form coefficient, since CRVE p-values have run optimistic by a
# factor of four to six throughout this project.
#
# Writes outputs/tables/two_margin.csv.
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

set.seed(20260924)
B_BOOT <- 9999L

cat("Loading analysis_data.RData...\n")
load("analysis_data.RData")
setDT(reduced_data_iv)

# --- Undivided pooled counts -------------------------------------------------
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

# --- First differences on a complete MSA x year grid ------------------------
yrs  <- seq(min(pred$YEAR), max(pred$YEAR))
grid <- CJ(MSA = sort(unique(pred$MSA)), YEAR = yrs)

pred <- merge(grid, pred[, .(MSA, YEAR, I_pred)], by = c("MSA","YEAR"), all.x = TRUE)
pred[is.na(I_pred), I_pred := 0]
setorder(pred, MSA, YEAR)
pred[, D_pred := I_pred - shift(I_pred, 1L, type = "lag"), by = MSA]

act <- merge(grid, act[, .(MSA, YEAR, I_act)], by = c("MSA","YEAR"), all.x = TRUE)
act[is.na(I_act), I_act := 0]
setorder(act, MSA, YEAR)
act[, D_act := I_act - shift(I_act, 1L, type = "lag"), by = MSA]

setnames(pred, "YEAR", "window_year")
setnames(act,  "YEAR", "window_year")

# --- Sample and window -------------------------------------------------------
d <- reduced_data_iv[immigr == 0 & MSA != 0 & MSA %in% selected_msas]
d[, y_HSplus := as.integer(Educ_Group != "non_HS")]

elig <- d[exposure_missing == 0, .(PERSON_ID, MSA, AGE, YEAR)]
elig[, birth_year := YEAR - AGE]
elig[, `:=`(start_year = birth_year + 15L,
            end_year   = pmin(birth_year + 25L, 2022L))]
elig <- elig[start_year <= end_year]
expanded <- elig[, .(window_year = seq(start_year, end_year)),
                 by = .(PERSON_ID, MSA)]

w <- merge(expanded, pred[, .(MSA, window_year, I_pred, D_pred)],
           by = c("MSA","window_year"), all.x = TRUE)
w <- merge(w, act[, .(MSA, window_year, I_act, D_act)],
           by = c("MSA","window_year"), all.x = TRUE)
agg <- w[, .(z_lev = mean(I_pred, na.rm = TRUE), z_dif = mean(D_pred, na.rm = TRUE),
             x_lev = mean(I_act,  na.rm = TRUE), x_dif = mean(D_act,  na.rm = TRUE)),
         by = PERSON_ID]

d <- merge(d, agg, by = "PERSON_ID", all.x = TRUE)
for (v in c("z_lev","z_dif","x_lev","x_dif")) d[is.na(get(v)), (v) := 0]
d <- droplevels(d)

cat(sprintf("\nSample: %s obs | %d MSAs | %s identifying\n",
            format(nrow(d), big.mark = ","), uniqueN(d$MSA),
            format(d[exposure_missing == 0, .N], big.mark = ",")))
cat(sprintf("cor(z_lev, z_dif) = %.4f  |  cor(x_lev, x_dif) = %.4f\n",
            d[exposure_missing == 0, cor(z_lev, z_dif)],
            d[exposure_missing == 0, cor(x_lev, x_dif)]))
sd_lev <- d[exposure_missing == 0, sd(z_lev)]
sd_dif <- d[exposure_missing == 0, sd(z_dif)]
cat(sprintf("SD among identifying obs: level %.1f | change %.1f\n", sd_lev, sd_dif))

CTRL <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
         Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA + AGE"

# --- Restricted wild cluster bootstrap for one coefficient ------------------
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

# Bootstrap a reduced-form coefficient, residualising on the controls AND on
# any other instrument in the model, so the test targets that margin alone.
rf_boot <- function(target, others) {
  rhs <- paste(c(others, sub("^\\s+", "", CTRL)), collapse = " + ")
  yt <- residuals(feols(as.formula(paste("y_HSplus ~", rhs)), data = d))
  zt <- residuals(feols(as.formula(paste(target,   "~", rhs)), data = d))
  g  <- as.character(d$MSA)
  list(rad = boot_p(yt, zt, g, B_BOOT),
       web = boot_p(yt, zt, g, B_BOOT, "webb"))
}

# --- 1-2. Reduced forms ------------------------------------------------------
cat("\n", strrep("-", 72), "\n", sep = "")
cat("REDUCED FORMS\n")

rf_rows <- list()
for (cfg in list(list(v = "z_lev", o = character(0), lab = "level alone"),
                 list(v = "z_dif", o = character(0), lab = "change alone"),
                 list(v = "z_lev", o = "z_dif",      lab = "level | change"),
                 list(v = "z_dif", o = "z_lev",      lab = "change | level"))) {
  rhs <- paste(c(cfg$v, cfg$o, sub("^\\s+", "", CTRL)), collapse = " + ")
  m  <- feols(as.formula(paste("y_HSplus ~", rhs)), data = d, vcov = ~MSA)
  ct <- m$coeftable[cfg$v, ]
  bb <- rf_boot(cfg$v, cfg$o)
  sdv <- if (cfg$v == "z_lev") sd_lev else sd_dif
  cat(sprintf("\n  %-16s %+.6e  SE %.4e  CRVE p %.4f  boot p %.4f (Webb %.4f)\n",
              cfg$lab, ct[["Estimate"]], ct[["Std. Error"]],
              ct[["Pr(>|t|)"]], bb$rad$p, bb$web$p))
  cat(sprintf("                   %+.4f pp per identifying SD\n",
              100 * ct[["Estimate"]] * sdv))
  rf_rows[[length(rf_rows)+1]] <- data.frame(
    model = "reduced form", spec = cfg$lab, margin = cfg$v,
    coef = ct[["Estimate"]], se = ct[["Std. Error"]],
    p_crve = ct[["Pr(>|t|)"]], p_boot = bb$rad$p, p_boot_webb = bb$web$p,
    pp_per_sd = 100 * ct[["Estimate"]] * sdv, stringsAsFactors = FALSE)
}

# --- 3. Two-endogenous 2SLS --------------------------------------------------
cat("\n", strrep("-", 72), "\n", sep = "")
cat("2SLS: two endogenous, two instruments (just-identified)\n")

m_iv <- feols(as.formula(sprintf(
  "y_HSplus ~ %s | x_lev + x_dif ~ z_lev + z_dif", CTRL)),
  data = d, vcov = ~MSA)
print(m_iv$coeftable[c("fit_x_lev", "fit_x_dif"), ])

cat("\nFirst stages:\n")
tryCatch(print(fitstat(m_iv, c("ivf", "ivwald"))),
         error = function(e) cat("  fitstat unavailable:", conditionMessage(e), "\n"))
cat("\nNote: with two endogenous regressors the per-equation F overstates\n")
cat("identification strength; Sanderson-Windmeijer conditional F is the\n")
cat("appropriate weak-instrument diagnostic and is not computed here.\n")

iv_rows <- do.call(rbind, lapply(c("fit_x_lev","fit_x_dif"), function(v) {
  ct <- m_iv$coeftable[v, ]
  data.frame(model = "2SLS (both)", spec = "joint", margin = v,
             coef = ct[["Estimate"]], se = ct[["Std. Error"]],
             p_crve = ct[["Pr(>|t|)"]], p_boot = NA_real_,
             p_boot_webb = NA_real_, pp_per_sd = NA_real_,
             stringsAsFactors = FALSE)
}))

res <- rbind(do.call(rbind, rf_rows), iv_rows)

cat("\n", strrep("=", 72), "\n", sep = "")
print(res, row.names = FALSE, digits = 4)

cat("\n", strrep("-", 72), "\n", sep = "")
cat("READING THIS\n\n")
cat("Compare 'level alone' to 'level | change', and 'change alone' to\n")
cat("'change | level'. If both coefficients survive jointly with the signs\n")
cat("they had alone, the two margins are separately identified and the sign\n")
cat("disagreement in 25 was omitted-variable bias -- each specification was\n")
cat("absorbing the other's variation through a -0.47 correlation. That is the\n")
cat("good case, and it is a more interesting result than either alone:\n")
cat("the level and the trajectory of immigration during adolescence have\n")
cat("opposite associations with completing high school.\n\n")
cat("If one collapses toward zero when both are included, that margin was\n")
cat("never carrying independent signal and the surviving one is the finding.\n\n")
cat("If both become unstable or the standard errors blow up, the -0.47\n")
cat("correlation is too strong to separate them at this effective cluster\n")
cat("count, and the honest report is one margin with the other disclosed as\n")
cat("an unresolved alternative.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/two_margin.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/two_margin.csv\n")
cat("13_two_margin.R complete.\n")
