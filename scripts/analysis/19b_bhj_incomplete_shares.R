# =============================================================================
# 19b_bhj_incomplete_shares.R  --  Does the estimate survive the Borusyak,
#                                  Hull and Jaravel sum-of-shares control?
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Added 2026-10-30.
#
# THE CONCERN
# \citet{borusyak2022quasi} distinguish "complete" from "incomplete" shares.
# When a unit's exposure shares sum to one, the shift-share instrument is a
# weighted AVERAGE of the shocks, so randomly assigned shocks produce a
# randomly assigned instrument. This design follows Card (2001) instead:
# share_{m,k} is MSA m's fraction of origin k's national population, so the
# shares sum to one across MSAs WITHIN an origin, not across origins within an
# MSA. The instrument is therefore a weighted SUM, and BHJ show that a unit
# with a larger sum of shares systematically receives a larger instrument value
# even under randomly assigned shocks -- which biases the estimate if that sum
# correlates with the error. Their prescription is to control for
#
#   S_m = sum_k share_{m,k}
#
# interacted with the fixed effects for the period of the shock (BHJ Practical
# Guide, Appendix A.3: with complete shares the control coincides with the
# period effects, but in the incomplete-shares case it must be interacted).
#
# WHY NOT JUST RENORMALISE THE SHARES
# Dividing by the within-MSA total gives the "local share" formulation, whose
# shares do sum to one. BHJ note (A.4) that controlling is the better solution,
# and here it is not merely better: S_m IS the scale information. Card's
# formulation makes the instrument a predicted HEADCOUNT -- a metro holding 20
# per cent of all Mexicans receives 0.20 times the national Mexican inflow.
# Renormalising turns it into a composition-weighted national average, under
# which two metros with the same origin mix and a hundredfold difference in
# immigrant population receive identical instrument values. For an outcome
# about how much immigration an individual encountered, the headcount is the
# object of interest, so the control is the principled fix.
#
# TWO VERSIONS, AND WHY BOTH ARE REPORTED
#   S_m x survey year.  Drops into the existing specification unchanged, since
#     YEAR effects are already present. But survey year is when a respondent
#     was enumerated, not the period of the shock, so this is not BHJ's control.
#   S_m x birth cohort.  The faithful analogue: the window W_i indexes which
#     arrival years enter an individual's exposure, and the cohort indexes the
#     window. This is the version that answers BHJ.
#
# Cohort slopes are entered as varying slopes (birth_year[[S_m]]), NOT as
# i(factor(birth_year), S_m). The two are algebraically identical, but the
# latter materialises 33 dense columns over 8.5 million rows and fails on a
# 32 GB machine ("cannot allocate vector of size 9.0 Gb"); varying slopes are
# absorbed in the demeaning instead. The double bracket omits the cohort fixed
# effect: single brackets would add it, and with YEAR effects already in the
# model a full set of cohort dummies spans the same space as a full set of AGE
# dummies, which is the over-control specification 12_final_spec_no_age_fe.R
# examined and rejected.
#
# HOW TO READ THE OUTPUT -- IN STANDARD ERRORS, NOT PER CENT
# The cohort control absorbs about half the identifying variation, so the
# controlled estimates are much less precise and their point estimates move
# accordingly. A percentage change is misleading here: the reduced form and the
# 2SLS move in OPPOSITE directions, which is arithmetic rather than evidence,
# since 2SLS = reduced form / first stage for a just-identified estimator and
# the first stage falls by half. The script therefore reports each movement in
# units of the controlled specification's own standard error.
#
# Writes outputs/tables/bhj_incomplete_shares.csv
#    and outputs/tables/table_bhj_incomplete_shares.tex
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")
source("scripts/analysis/00b_table_helper.R")

library(data.table)
library(fixest)

stopifnot(exists("write_dissertation_table"))
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

if (!exists("reduced_data_iv") || !exists("shift_share_data")) {
  cat("Loading analysis_data.RData...\n")
  load("analysis_data.RData")
}
setDT(shift_share_data); setDT(reduced_data_iv)

# --- S_m ---------------------------------------------------------------------
# shift_share_data is MSA x BPLD x YEAR and the share is constant within
# MSA x BPLD, so it is deduplicated before summing across origins. This uses
# the same share object the instrument itself was built from rather than
# rebuilding it from the microdata.
Sm <- unique(shift_share_data[, .(MSA, BPLD, share)])[, .(S_m = sum(share)), by = MSA]

# The share grid spans every MSA observed in the 1990 wave, including the
# non-metropolitan residual MSA == 0, which holds a slice of nearly every
# origin and so carries by far the largest S_m. The estimation sample excludes
# it, so the distribution that describes the design is the one over the
# selected metros.
Sm_est <- Sm[MSA %in% selected_msas & MSA != 0]

cat("\n", strrep("=", 72), "\n", sep = "")
cat("SUM OF SHARES\n\n")
cat(sprintf("  MSAs in the share grid        %d\n", nrow(Sm)))
cat(sprintf("  MSAs in the estimation sample %d\n", nrow(Sm_est)))
cat(sprintf("\n  S_m over the estimation MSAs:\n"))
cat(sprintf("    min %.3f   p25 %.3f   median %.3f   p75 %.3f   max %.3f   mean %.3f\n",
            min(Sm_est$S_m), quantile(Sm_est$S_m, .25), median(Sm_est$S_m),
            quantile(Sm_est$S_m, .75), max(Sm_est$S_m), mean(Sm_est$S_m)))
cat(sprintf("  largest S_m in the full grid (incl. MSA 0): %.3f\n", max(Sm$S_m)))

# Metropolitan names, to check the tail of S_m is the gateway ranking and not
# an artefact. MSA_Name is not carried on the microdata; it lives in
# top175_list_named, which is one table per year, so it is stacked first.
top_name <- NA_character_
if (exists("top175_list_named")) {
  msa_names <- unique(rbindlist(lapply(top175_list_named,
    function(x) as.data.table(x)[, .(MSA = as.character(MSA), MSA_Name)])))
  Sm_top <- merge(Sm_est[order(-S_m)][1:10][, .(MSA = as.character(MSA), S_m)],
                  msa_names, by = "MSA", all.x = TRUE)[order(-S_m)]
  cat("\n  Largest S_m among the estimation MSAs:\n")
  print(Sm_top, row.names = FALSE)
  if (!is.na(Sm_top$MSA_Name[1])) top_name <- Sm_top$MSA_Name[1]
} else {
  cat("\n  top175_list_named not available; S_m tail reported without names.\n")
}

# --- Appendix exhibit: the tail of S_m, named ---------------------------------
# The sum-of-shares control is easier to read against named metropolitan areas
# than against a bare maximum of 48. This table exists so the claim that S_m
# orders metropolitan areas by gateway intensity is supported by an exhibit
# rather than asserted, and so a reader can see that the largest value is the
# dominant immigrant gateway rather than a fault in the share construction.
if (!is.na(top_name)) {
  .f3 <- function(x) formatC(x, format = "f", digits = 3)
  fwrite(Sm_top, "outputs/tables/sm_gateways.csv")

  gw <- data.frame(
    `Metropolitan area` = c(Sm_top$MSA_Name, "@MIDRULE@",
                            "Median of the 175 estimation metropolitan areas"),
    `$S_m$` = c(.f3(Sm_top$S_m), "", .f3(median(Sm_est$S_m))),
    check.names = FALSE, stringsAsFactors = FALSE)

  write_dissertation_table(
    gw,
    caption = "Metropolitan Areas with the Largest Sum of Exposure Shares",
    label = "tab:sm_gateways",
    notes = paste0(
      "$S_m = \\sum_k s_{m,k}$ is the sum, across the ",
      nrow(unique(shift_share_data[, .(BPLD)])), " origins entering the inflow ",
      "series, of metropolitan area $m$'s 1990 share of each origin's national ",
      "population. Because the shares follow \\citet{card2001immigrant} and are ",
      "normalised by each origin's national total, they sum to one across ",
      "metropolitan areas within an origin and need not sum to anything in ",
      "particular across origins within a metropolitan area; $S_m$ therefore ",
      "exceeds one for any metropolitan area holding a substantial fraction of ",
      "many origins' national populations, and is best read as an index of how ",
      "much of total U.S. immigration a metropolitan area absorbs. The ordering ",
      "is the familiar ranking of immigrant gateways, which is the reason the ",
      "largest value is far above the median rather than evidence of an error ",
      "in the shares. Honolulu enters the ten largest despite its size because ",
      "it holds a large share of several Pacific-origin populations. $S_m$ is ",
      "the quantity controlled for in Table~\\ref{tab:bhj_incomplete_shares}, ",
      "following \\citet{borusyak2022quasi}; the ten metropolitan areas listed ",
      "here are those whose exposure that control most reduces."),
    align = "lc",
    file = "outputs/tables/table_sm_gateways.tex"
  )
}

# --- Estimation sample, identical to 09_main_results.R -----------------------
d <- reduced_data_iv[immigr == 0 & MSA != 0 & MSA %in% selected_msas]
d[, `:=`(
  birth_year = YEAR - AGE,
  y_HSplus   = as.integer(Educ_Group != "non_HS"),
  z = avg_exposure_nodiv,
  x = avg_actual_exposure_nodiv
)]
d[Sm, on = "MSA", S_m := i.S_m]
if (anyNA(d$S_m))
  stop("S_m failed to merge for ", sum(is.na(d$S_m)), " rows -- check that MSA ",
       "has the same type in shift_share_data and reduced_data_iv.")
d <- droplevels(d)

id   <- d[exposure_missing == 0]
sd_z <- id[, sd(z)]
cat(sprintf("\n  observations %s   identifying %s   SD of instrument %.1f\n",
            format(nrow(d), big.mark = ","), format(nrow(id), big.mark = ","), sd_z))

CTRL <- paste("Black + SEX + Marital + ChildAgeGroup + POVERTY +",
              "Mom_Educ + Dad_Educ + exposure_missing")

# --- How much instrument variation each specification leaves ------------------
# Residual SD of the instrument after the specification's own absorption,
# over the identifying observations. This is what sets how precisely the
# coefficient can be estimated, and it is the quantity against which the
# movement in the point estimates has to be judged.
retained <- function(fe) {
  m <- feols(as.formula(sprintf("z ~ 1 | %s", fe)), data = id)
  s <- sd(resid(m)); rm(m); gc()
  s
}

SPECS <- list(
  base   = list(fe = "YEAR + MSA",
                label = "Baseline"),
  syear  = list(fe = "MSA + YEAR[S_m]",
                label = "$S_m \\times$ survey year"),
  cohort = list(fe = "YEAR + MSA + birth_year[[S_m]]",
                label = "$S_m \\times$ birth cohort")
)

cat("\n", strrep("-", 72), "\n", sep = "")
cat("INSTRUMENT VARIATION RETAINED (identifying observations)\n\n")
for (nm in names(SPECS)) {
  SPECS[[nm]]$sd_retained <- retained(SPECS[[nm]]$fe)
  cat(sprintf("  %-28s %8.1f   (%.1f%% of the raw %.1f)\n",
              gsub("\\$|\\\\times", "", SPECS[[nm]]$label),
              SPECS[[nm]]$sd_retained,
              100 * SPECS[[nm]]$sd_retained / sd_z, sd_z))
}

# --- Estimation ---------------------------------------------------------------
est <- function(fe, label) {
  fs <- feols(as.formula(sprintf("x ~ z + %s | %s", CTRL, fe)), d, vcov = ~MSA)
  rf <- feols(as.formula(sprintf("y_HSplus ~ z + %s | %s", CTRL, fe)), d, vcov = ~MSA)
  iv <- feols(as.formula(sprintf("y_HSplus ~ %s | %s | x ~ z", CTRL, fe)), d, vcov = ~MSA)
  ct <- iv$coeftable["fit_x", ]
  out <- list(
    label       = label,
    first_stage = unname(coef(fs)[["z"]]),
    fs_se       = fs$coeftable["z", "Std. Error"],
    rf          = unname(coef(rf)[["z"]]),
    rf_se       = rf$coeftable["z", "Std. Error"],
    rf_t        = rf$coeftable["z", "t value"],
    rf_p        = rf$coeftable["z", "Pr(>|t|)"],
    tsls        = ct[["Estimate"]],
    tsls_se     = ct[["Std. Error"]],
    tsls_t      = ct[["t value"]],
    tsls_p      = ct[["Pr(>|t|)"]]
  )
  out$pp_per_sd <- 100 * out$rf * sd_z
  cat(sprintf("\n  %s\n", label))
  cat(sprintf("    first stage   %+.4f   SE %.4f\n", out$first_stage, out$fs_se))
  cat(sprintf("    reduced form  %+.4e  SE %.3e  t %6.2f  p %.4f\n",
              out$rf, out$rf_se, out$rf_t, out$rf_p))
  cat(sprintf("    2SLS          %+.4e  SE %.3e  t %6.2f  p %.4f\n",
              out$tsls, out$tsls_se, out$tsls_t, out$tsls_p))
  cat(sprintf("    %+.2f pp per SD of exposure\n", out$pp_per_sd))
  rm(fs, rf, iv); gc()
  out
}

cat("\n", strrep("-", 72), "\n", sep = "")
cat("ESTIMATES\n")
for (nm in names(SPECS)) {
  r <- est(SPECS[[nm]]$fe, gsub("\\$|\\\\times", "", SPECS[[nm]]$label))
  SPECS[[nm]] <- c(SPECS[[nm]], r[setdiff(names(r), "label")])
}

# --- Movement, in standard errors of the controlled specification ------------
cat("\n", strrep("-", 72), "\n", sep = "")
cat("MOVEMENT AGAINST THE BASELINE\n")
cat("Read in standard errors, not per cent: the controls cost identifying\n")
cat("variation, so the controlled estimates are less precise, and the reduced\n")
cat("form and 2SLS necessarily move in opposite directions when the first\n")
cat("stage falls.\n\n")
b <- SPECS$base
for (nm in c("syear", "cohort")) {
  s <- SPECS[[nm]]
  s$rf_move_se   <- (s$rf   - b$rf)   / s$rf_se
  s$tsls_move_se <- (s$tsls - b$tsls) / s$tsls_se
  SPECS[[nm]] <- s
  cat(sprintf("  %-28s reduced form %+5.2f SE   2SLS %+5.2f SE   first stage %.3f -> %.3f\n",
              gsub("\\$|\\\\times", "", s$label),
              s$rf_move_se, s$tsls_move_se, b$first_stage, s$first_stage))
}

# --- CSV ---------------------------------------------------------------------
res <- rbindlist(lapply(names(SPECS), function(nm) {
  s <- SPECS[[nm]]
  data.table(spec = nm,
             label = gsub("\\$|\\\\times", "", s$label),
             fe = s$fe,
             sd_retained = s$sd_retained,
             sd_retained_pct = 100 * s$sd_retained / sd_z,
             first_stage = s$first_stage, fs_se = s$fs_se,
             reduced_form = s$rf, rf_se = s$rf_se, rf_t = s$rf_t, rf_p = s$rf_p,
             tsls = s$tsls, tsls_se = s$tsls_se, tsls_t = s$tsls_t, tsls_p = s$tsls_p,
             pp_per_sd = s$pp_per_sd,
             rf_move_se   = if (is.null(s$rf_move_se))   NA_real_ else s$rf_move_se,
             tsls_move_se = if (is.null(s$tsls_move_se)) NA_real_ else s$tsls_move_se)
}), fill = TRUE)
fwrite(res, "outputs/tables/bhj_incomplete_shares.csv")
cat("\nWrote outputs/tables/bhj_incomplete_shares.csv\n")

# --- LaTeX -------------------------------------------------------------------
fe2 <- function(x) formatC(x, format = "e", digits = 2)
f3  <- function(x) formatC(x, format = "f", digits = 3)
f2  <- function(x) formatC(x, format = "f", digits = 2)
f1  <- function(x) formatC(x, format = "f", digits = 1)
stars <- function(p) if (p < .01) "$^{***}$" else if (p < .05) "$^{**}$" else
                     if (p < .10) "$^{*}$" else ""
col <- function(nm) {
  s <- SPECS[[nm]]
  c(paste0(fe2(s$rf), stars(s$rf_p)), paste0("(", fe2(s$rf_se), ")"),
    paste0(fe2(s$tsls), stars(s$tsls_p)), paste0("(", fe2(s$tsls_se), ")"),
    f3(s$first_stage),
    "",
    f2(s$pp_per_sd),
    f1(s$sd_retained), f1(100 * s$sd_retained / sd_z),
    "",
    if (is.null(s$rf_move_se))   "---" else f2(s$rf_move_se),
    if (is.null(s$tsls_move_se)) "---" else f2(s$tsls_move_se))
}

tbl <- data.frame(
  ` ` = c("Reduced form", "\\quad MSA-clustered SE",
          "2SLS", "\\quad MSA-clustered SE",
          "First-stage coefficient",
          "@MIDRULE@",
          "Effect per SD (pp)",
          "Instrument SD retained",
          "\\quad per cent of raw SD",
          "@MIDRULE@",
          "Movement in reduced form (SE)",
          "Movement in 2SLS (SE)"),
  Baseline = col("base"),
  check.names = FALSE, stringsAsFactors = FALSE
)
tbl[[SPECS$syear$label]]  <- col("syear")
tbl[[SPECS$cohort$label]] <- col("cohort")

write_dissertation_table(
  tbl,
  caption = paste("Sensitivity to the \\citet{borusyak2022quasi}",
                  "Sum-of-Shares Control"),
  label = "tab:bhj_incomplete_shares",
  notes = paste0(
    "$S_m = \\sum_k s_{m,k}$ is metropolitan area $m$'s 1990 exposure shares ",
    "summed across origins; its distribution across the panel is reported in ",
    "Table~\\ref{tab:sm_gateways}. Column (1) is the headline specification of ",
    "Table~\\ref{tab:education_main}. Columns (2) and (3) add $S_m$ interacted ",
    "with survey-year and with birth-cohort effects respectively, the latter as ",
    "cohort-specific slopes without cohort fixed effects. ``Instrument SD ",
    "retained'' is the standard deviation of the instrument surviving each ",
    "column's own absorption, over the identifying observations. Movement is ",
    "reported in standard errors of the column's own specification: the controls ",
    "absorb identifying variation, so the controlled estimates are less precise, ",
    "and because the estimator is just-identified a change in the first stage ",
    "moves the reduced form and the 2SLS in opposite directions. Standard errors ",
    "are clustered on metropolitan area. $^{***}$ $p<0.01$, $^{**}$ $p<0.05$, ",
    "$^{*}$ $p<0.10$."),
  align = "lccc",
  file = "outputs/tables/table_bhj_incomplete_shares.tex"
)

cat("\n", strrep("-", 72), "\n", sep = "")
