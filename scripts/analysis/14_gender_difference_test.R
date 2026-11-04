# =============================================================================
# 14_gender_difference_test.R  --  Does the effect differ between men and women?
# Renumbered 2026-10-21 from 31b_gender_difference_test.R (archive-and-renumber pass).
# Part of the modular Immigration_Card_Analysis.R pipeline.
# REBUILT 2026-09-30 on the final specification.
#
# SUPERSEDES 31_gender_heterogeneity.R, which should be archived. That script
# reported a pooled `exposure x female` interaction as the test of the
# difference. It is not one: the pooled model forces every other coefficient --
# Black, Marital, ChildAgeGroup, POVERTY, Mom_Educ, Dad_Educ, and both sets of
# fixed effects -- to be common across gender, while separate subsample fits let
# all of them differ. Returns to parental education and marital status plainly
# do differ by gender, and that misspecification loads onto the interaction
# term. On the college margin it returned +6.470e-03 against an actual
# subsample difference of -1.642e-03: opposite sign, four times the magnitude.
#
# WHAT THIS DOES INSTEAD
# Estimates men and women separately, each on its own controls and fixed
# effects, then tests the difference with a standard error that accounts for
# the two subsamples sharing MSA clusters. The naive alternative,
# sqrt(se_men^2 + se_women^2), assumes independence; men and women in the same
# metro are in the same cluster, so the two subsample estimators are correlated
# and that formula misstates the SE of the difference. The influence
# contribution of cluster g to the difference is
#
#     psi_g = s_g^W / ZZ_W  -  s_g^M / ZZ_M
#
# with both subsamples entering the SAME psi_g, which is how the dependence
# enters. Var(b_W - b_M) = sum_g psi_g^2. Both the correct and the naive SE are
# reported so the size of the difference is visible.
#
# INFERENCE. MSA-clustered, conventional p-values. No bootstrap -- that is the
# later chapter's material. The conclusions do not need it: on the previous
# specification the HS-margin difference carried a cluster-robust t of -0.513
# and the college difference -6.685, and the bootstrap only confirmed both.
#
# A THREAT SPECIFIC TO THIS CONTRAST, WHICH THE POOLED ESTIMATE DOES NOT FACE
# Women complete secondary schooling at higher rates and somewhat earlier than
# men. With no age restriction the sample includes cohorts whose attainment is
# still being realised, so the censoring is GENDER-DIFFERENTIAL -- it does not
# fall equally on the two subsamples. That can manufacture or mask a difference
# in the exposure coefficient. The pooled headline is immune, because there the
# censoring is absorbed by exposure_missing, YEAR and MSA rather than
# interacting with the contrast being estimated.
#
# So this script reports non_HS and college rates by sex within age bands
# BEFORE the estimates. If those rates diverge by sex among the young cohorts,
# any measured gender difference is partly composition and should be described
# that way. Read the diagnostic first.
#
# Writes outputs/tables/gender_difference_test.csv.
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
           "exposure_missing", "Educ_Group", "SEX", "AGE")
.miss <- setdiff(.need, names(all_natives_175))
if (length(.miss))
  stop("analysis_data.RData is stale or incomplete -- missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R first.")
rm(.need, .miss)

d <- all_natives_175[MSA != 0]
d[, `:=`(
  y_HSplus   = as.integer(Educ_Group != "non_HS"),
  y_CollPlus = as.integer(Educ_Group %in% c("College", "Grad")),
  z = avg_exposure_nodiv,
  x = avg_actual_exposure_nodiv,
  sex_i = as.integer(as.character(SEX))          # 1 = male, 2 = female
)]
stopifnot(all(d$sex_i %in% c(1L, 2L)))
d[, female := as.integer(sex_i == 2L)]
d <- droplevels(d)

sd_z <- d[exposure_missing == 0, sd(z)]

cat(sprintf("\n%s\n", strrep("=", 70)))
cat("Section 31b: gender difference in the exposure effect\n")
cat(sprintf("%s\n\n", strrep("=", 70)))
cat(sprintf("  observations   %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  men            %s\n", format(d[female == 0, .N], big.mark = ",")))
cat(sprintf("  women          %s\n", format(d[female == 1, .N], big.mark = ",")))
cat(sprintf("  identifying    %s\n", format(d[exposure_missing == 0, .N], big.mark = ",")))
cat(sprintf("  MSA clusters   %d\n", uniqueN(d$MSA)))
cat(sprintf("  SD of instrument among identifying: %.1f\n", sd_z))

# --- READ THIS FIRST: is the censoring gender-differential? ------------------
cat("\n", strrep("-", 70), "\n", sep = "")
cat("CENSORING BY SEX AND AGE BAND\n")
cat("The gender contrast is only clean if incomplete schooling falls equally\n")
cat("on men and women. Compare the non_HS columns within each band.\n\n")
d[, age_band := cut(AGE, c(0, 17, 21, 25, 30, 40, 68, Inf),
                    labels = c("<=17","18-21","22-25","26-30","31-40","41-68","69+"))]
cens <- dcast(
  d[, .(pct_non_HS = 100 * mean(Educ_Group == "non_HS"),
        pct_coll   = 100 * mean(y_CollPlus), n = .N),
    by = .(age_band, sex = fifelse(female == 1L, "women", "men"))],
  age_band ~ sex, value.var = c("n", "pct_non_HS", "pct_coll"))
cens[, `:=`(gap_non_HS = pct_non_HS_women - pct_non_HS_men,
            gap_coll   = pct_coll_women   - pct_coll_men)]
print(cens[order(age_band)], row.names = FALSE, digits = 3)
cat("\ngap_non_HS is women minus men. A gap that is large among the young\n")
cat("bands and small among the old ones is differential censoring, not a\n")
cat("difference in completed attainment.\n")
d[, age_band := NULL]

# --- Machinery ---------------------------------------------------------------
# SEX is dropped from the subsample controls (constant within each), and
# everything else is free to differ across the two fits -- which is the point.
CTRL_SUB <- "Black + Marital + ChildAgeGroup + POVERTY +
             Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

# Age-fixed-effects variant, run alongside as a DIAGNOSTIC (added 2026-09-30).
#
# WHY. The censoring table above shows the non_HS gap between women and men is
# largest among the young bands (-4.20 at 18-21) and shrinks with age (-1.27 at
# 41-68), so incomplete schooling does not fall equally on the two subsamples.
# That alone can generate a gender difference in the exposure coefficient, and
# in the direction observed: the undivided inflow series declines as n_waves
# falls, so younger cohorts carry LOWER z while also having HIGHER non_HS. Women
# hold a larger non_HS advantage at low z than at high z, which flattens their
# slope relative to men's and yields a smaller positive coefficient. Exactly
# what the baseline rows show.
#
# Age fixed effects absorb that cohort gradient. This is a CONTROL, not a sample
# restriction -- no observations are dropped, so it does not reintroduce the age
# cut that has been ruled out.
#
# HOW TO READ THE PAIR. If the difference collapses under age fixed effects, it
# was differential censoring and the gender material comes out of Chapter 2
# cleanly, costing nothing. If it survives, the difference is real and becomes
# evidence for the later chapter's argument that age structure has to be
# modelled rather than absorbed -- a gender contrast that only appears once
# cohort effects are handled is a stronger case for that than the pooled
# magnitude is.
CTRL_SUB_AGE <- "Black + Marital + ChildAgeGroup + POVERTY +
                 Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA + AGE"

msas     <- sort(unique(d$MSA))
msa_lvls <- as.character(msas)
csum <- function(v, g) {
  out <- tapply(v, factor(g, levels = msa_lvls), sum)
  out[is.na(out)] <- 0
  as.numeric(out)
}

# Residualise a subsample and return its aligned cluster-level pieces.
prep <- function(dd, outcome, ctrl = CTRL_SUB) {
  dd    <- droplevels(dd)
  other <- sub("^\\s+", "", ctrl)
  yt <- residuals(feols(as.formula(paste(outcome, "~", other)), data = dd))
  zt <- residuals(feols(as.formula(paste("z ~",     other)), data = dd))
  g  <- as.character(dd$MSA)

  S  <- csum(zt * yt, g); Q <- csum(zt^2, g)
  ZZ <- sum(Q); b <- sum(S) / ZZ
  s  <- csum(zt * (yt - b * zt), g)

  list(b = b, se = sqrt(sum(s^2)) / ZZ, s = s, ZZ = ZZ, n = nrow(dd),
       base = mean(dd[[outcome]]),
       sd_idn = dd[exposure_missing == 0, sd(z)])
}

run <- function(outcome, ctrl = CTRL_SUB, spec = "baseline (no age FE)") {
  cat("\n", strrep("-", 70), "\n", sep = "")
  cat(sprintf("%s  --  %s\n", outcome, spec))

  M <- prep(d[female == 0], outcome, ctrl)
  W <- prep(d[female == 1], outcome, ctrl)

  diff   <- W$b - M$b
  psi    <- W$s / W$ZZ - M$s / M$ZZ
  se_ok  <- sqrt(sum(psi^2))
  se_bad <- sqrt(M$se^2 + W$se^2)
  t_obs  <- diff / se_ok
  p_obs  <- 2 * pnorm(-abs(t_obs))

  cat(sprintf("  men        %+.6e  (SE %.4e)  %+7.3f pp/SD  base %.4f  n %s\n",
              M$b, M$se, 100 * M$b * sd_z, M$base, format(M$n, big.mark = ",")))
  cat(sprintf("  women      %+.6e  (SE %.4e)  %+7.3f pp/SD  base %.4f  n %s\n",
              W$b, W$se, 100 * W$b * sd_z, W$base, format(W$n, big.mark = ",")))
  cat(sprintf("  difference %+.6e\n", diff))
  cat(sprintf("    SE accounting for shared MSAs  %.6e   t %+7.3f   p %.4f\n",
              se_ok, t_obs, p_obs))
  cat(sprintf("    SE assuming independence       %.6e   (ratio %.3f)\n",
              se_bad, se_bad / se_ok))
  cat(sprintf("    95%% CI on the difference       [%+.4e, %+.4e]\n",
              diff - 1.96 * se_ok, diff + 1.96 * se_ok))

  data.frame(outcome = outcome, spec = spec,
             coef_men = M$b, se_men = M$se, base_men = M$base,
             coef_women = W$b, se_women = W$se, base_women = W$base,
             difference = diff, se_shared = se_ok, se_naive = se_bad,
             t = t_obs, p = p_obs,
             pp_sd_men = 100 * M$b * sd_z, pp_sd_women = 100 * W$b * sd_z,
             stringsAsFactors = FALSE)
}

SPEC_BASE <- "baseline (no age FE)"
SPEC_AGE  <- "with age FE (diagnostic)"
res <- rbindlist(c(
  lapply(c("y_HSplus", "y_CollPlus"), run, spec = SPEC_BASE),
  lapply(c("y_HSplus", "y_CollPlus"), run,
         ctrl = CTRL_SUB_AGE, spec = SPEC_AGE)
))
setorder(res, outcome, spec)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("SUMMARY\n\n")
print(res[, .(outcome, spec, coef_men, coef_women, difference, se_shared, t, p)],
      row.names = FALSE, digits = 3)

cat("\nThe age-FE rows are a diagnostic, not a competing estimate. Compare the\n")
cat("difference and its t within each outcome:\n\n")
for (oc in unique(res$outcome)) {
  # Exact match, not %like%: "baseline (no age FE)" contains "age FE", so a
  # substring match on "age FE" returns BOTH rows and sprintf then vectorises
  # over them, printing each comparison twice with mismatched values.
  b <- res[outcome == oc & spec == SPEC_BASE, ]
  a <- res[outcome == oc & spec == SPEC_AGE,  ]
  stopifnot(nrow(b) == 1L, nrow(a) == 1L)
  cat(sprintf("  %-12s baseline %+.4e (t %+6.2f)  ->  with age FE %+.4e (t %+6.2f)",
              oc, b$difference, b$t, a$difference, a$t))
  cat(sprintf("   retains %.0f%%\n", 100 * a$difference / b$difference))
}

cat("\n", strrep("-", 70), "\n", sep = "")
cat("READING THIS\n\n")
cat("`difference` equals coef_women - coef_men by construction, which the\n")
cat("superseded pooled interaction in 31 did not. Compare se_naive to\n")
cat("se_shared: the ratio is how wrong an independence assumption would be,\n")
cat("and it reflects how strongly the two subsample estimates move together\n")
cat("across metros. A large ratio means the common metro-level variation\n")
cat("cancels in the difference, which makes the contrast more precisely\n")
cat("estimated than either level -- not a mistake, the same logic that makes\n")
cat("a difference-in-differences sharper than either arm.\n\n")
cat("Then go back to the censoring table. If non_HS gaps between women and men\n")
cat("are large among the young bands, a measured difference in the exposure\n")
cat("coefficient is partly differential censoring rather than a difference in\n")
cat("response, and the text should say so. If the gaps are stable across\n")
cat("bands, the contrast is clean and can be reported as a behavioural\n")
cat("difference.\n\n")
cat("THE AGE-FE ROWS ARE THE ADJUDICATION. Age fixed effects absorb the cohort\n")
cat("gradient that makes the censoring gender-differential, without dropping\n")
cat("any observations.\n\n")
cat("  Difference collapses under age FE  ->  it was censoring. Gender comes\n")
cat("  out of Chapter 2 and nothing is lost, because there was no behavioural\n")
cat("  finding to report.\n\n")
cat("  Difference survives under age FE   ->  it is real, and it becomes\n")
cat("  evidence for the later chapter: a gender contrast that is only visible\n")
cat("  once cohort effects are handled is a sharper argument that age\n")
cat("  structure must be modelled than the pooled magnitude is. Report it\n")
cat("  there rather than in Chapter 2, since it would be conditional on a\n")
cat("  specification Chapter 2 does not use.\n\n")
cat("A CAUTION FOR THE WRITE-UP. Under age fixed effects the difference is\n")
cat("precisely estimated while neither gender's LEVEL is distinguishable from\n")
cat("zero -- the common metro-level variation cancels in the contrast but not\n")
cat("in either arm. So the supportable claim is that the RESPONSE DIFFERS BY\n")
cat("GENDER. It is not that immigration raises attainment for men and lowers\n")
cat("it for women; neither of those is established. Check the level t-values\n")
cat("before writing any directional sentence about either group.\n\n")
cat("Note what a null on the HS margin would rule out. An opportunity-cost\n")
cat("channel -- immigration raising the return to young native labour and\n")
cat("pulling marginal students out of school -- predicts the effect should\n")
cat("concentrate in men, whose teen employment is far more exposed to\n")
cat("construction, manufacturing and agriculture. Equal effects point instead\n")
cat("to school-side mechanisms: class size, resource dilution, peer\n")
cat("composition, which operate on everyone in the affected schools.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/gender_difference_test.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/gender_difference_test.csv\n")
cat("14_gender_difference_test.R complete.\n")
