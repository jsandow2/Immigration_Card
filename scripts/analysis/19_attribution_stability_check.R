# =============================================================================
# 19_attribution_stability_check.R  --  Robustness: how much does MSA
# Renumbered 2026-10-21 from 34_attribution_stability_check.R (archive-and-renumber pass).
#                                       attribution move with time since
#                                       arrival?
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# THE LIMITATION
# IPUMS reports one geography per respondent: metropolitan area of residence
# AT SURVEY. There is no "MSA at time of arrival" variable. The actual-inflow
# series therefore credits an arrival to the metro where the person was
# eventually enumerated. For someone who arrived in 1995 and was surveyed in
# 2019, those may differ. The instrument is unaffected -- it uses NATIONAL
# arrivals by origin, which do not depend on where anyone later moves -- so
# only the endogenous variable carries the error.
#
# WHY IT MIGHT NOT BE BENIGN
# Measurement error in an endogenous regressor is ordinarily absorbed by the
# first stage. The reason to check is that this error may not be independent
# of the instrument: secondary migration follows established co-national
# communities, which is what the 1990 shares measure. If arrivals who land
# elsewhere later move toward high-share metros, the error correlates with the
# instrument.
#
# WHY NOT SIMPLY RESTRICT THE SERIES AND RE-ESTIMATE
# That was tried in 33_arrival_recency_check.R and it does not work. The
# survey calendar is not annual -- waves are 1990, 2000, 2006-2019 and
# 2021-2023 -- so a rule like "observed within two years of arrival" removes
# whole arrival years rather than misattributed individuals. Arrival years
# 1991-1997 and 2001-2003 vanish entirely, exposure windows spanning the
# 1990s collapse to near-zero, and the resulting first-stage decline measures
# the calendar rather than attribution. See that script's header.
#
# WHAT THIS SCRIPT DOES INSTEAD
# Asks the question only where the data can answer it, and descriptively
# rather than by re-estimating. For arrival years 2006-2019, where annual ACS
# coverage makes a recency rule feasible, build the (MSA, arrival year) inflow
# vector twice:
#
#   NEAR  counts only respondents observed within RECENT_MAX years of arrival
#   ALL   counts every respondent reporting that arrival year
#
# If attribution error is small, the two agree. If secondary migration moves
# people systematically toward metros with large 1990 co-national communities,
# then ALL will overstate inflows to high-share metros relative to NEAR, and
# the discrepancy will correlate with the shares. That correlation is the
# quantity of interest and it is measured directly.
#
# Writes outputs/tables/attribution_stability.csv
# =============================================================================

# Clears the session but KEEPS bpld_lookup. It is built from the IPUMS DDI in
# 08_load_analysis_subsets.R and cannot be rebuilt without re-running 08.
rm(list = setdiff(ls(), "bpld_lookup")); gc()
source("scripts/analysis/00_setup_analysis.R")

library(data.table)

RECENT_MAX <- 2L                 # years since arrival for the NEAR series
YR_LO      <- 2006L              # first arrival year with annual coverage
YR_HI      <- 2019L              # last before the 2020 gap

cat("Loading analysis_data.RData...\n")
load("analysis_data.RData")
setDT(reduced_data_iv)

arr <- reduced_data_iv[immigr == 1 & YRIMMIG > 0 & MSA != 0 &
                         MSA %in% selected_msas,
                       .(wave = as.integer(YEAR), YRIMMIG = as.integer(YRIMMIG),
                         MSA = as.integer(MSA), BPLD = as.integer(BPLD))]
arr[, yrs_since := wave - YRIMMIG]

# --- Confirm the restriction is feasible on this year range -----------------
# The premise of the whole script. If some arrival year in the range lacks a
# wave within RECENT_MAX years, it drops out of NEAR entirely and the
# comparison is contaminated in exactly the way 33 was.
waves <- sort(unique(arr$wave))
cat(sprintf("\nSurvey waves present: %s\n", paste(waves, collapse = ", ")))
feas <- data.table(YRIMMIG = YR_LO:YR_HI)
feas[, n_near_waves := vapply(YRIMMIG, function(y)
  sum(waves >= y & waves <= y + RECENT_MAX), integer(1))]
cat(sprintf("\nArrival years %d-%d, waves available within %d years of arrival:\n",
            YR_LO, YR_HI, RECENT_MAX))
print(feas, row.names = FALSE)
if (feas[, any(n_near_waves == 0L)])
  stop("Some arrival year in the range has no wave within ", RECENT_MAX,
       " years. The NEAR series would drop it entirely and this comparison ",
       "would repeat the flaw in 33_arrival_recency_check.R.")
cat("\n  Every arrival year in the range is observable under the rule, so the\n")
cat("  comparison below reflects WHO is counted, not WHICH YEARS survive.\n")

a <- arr[YRIMMIG >= YR_LO & YRIMMIG <= YR_HI]
cat(sprintf("\nArrival records in range: %s | within %d years: %s (%.1f%%)\n",
            format(nrow(a), big.mark = ","), RECENT_MAX,
            format(a[yrs_since >= 0 & yrs_since <= RECENT_MAX, .N], big.mark = ","),
            100 * a[, mean(yrs_since >= 0 & yrs_since <= RECENT_MAX)]))

# --- The two inflow vectors -------------------------------------------------
cell_all  <- a[, .(infl_all = .N), by = .(MSA, YRIMMIG)]
cell_near <- a[yrs_since >= 0 & yrs_since <= RECENT_MAX,
               .(infl_near = .N), by = .(MSA, YRIMMIG)]
cells <- merge(cell_all, cell_near, by = c("MSA", "YRIMMIG"), all.x = TRUE)
cells[is.na(infl_near), infl_near := 0L]

# Shares within arrival year, which is the quantity the instrument predicts.
# Raw counts differ by construction because NEAR uses fewer waves; the
# distribution ACROSS metros is what attribution error would distort.
cells[, `:=`(sh_all  = infl_all  / sum(infl_all),
             sh_near = infl_near / sum(infl_near)), by = YRIMMIG]

cat("\n", strrep("=", 72), "\n", sep = "")
cat("1. DO THE TWO SERIES AGREE ON WHERE ARRIVALS GO?\n\n")
cat(sprintf("  (MSA, arrival year) cells: %s\n", format(nrow(cells), big.mark = ",")))
cat(sprintf("  correlation of within-year shares : %.4f\n",
            cells[, cor(sh_all, sh_near)]))
cat(sprintf("  rank correlation                  : %.4f\n",
            cells[, cor(sh_all, sh_near, method = "spearman")]))
cat("\n  By arrival year:\n")
print(cells[, .(cells = .N, cor = round(cor(sh_all, sh_near), 4),
                mean_abs_diff_pp = round(100 * mean(abs(sh_all - sh_near)), 4)),
            by = YRIMMIG][order(YRIMMIG)], row.names = FALSE)

# --- Step 2: does the discrepancy correlate with the 1990 shares? -----------
# The question that matters. A metro's share of a given origin's 1990
# population is what drives its predicted inflow, so if ALL overstates
# arrivals to high-share metros relative to NEAR, the attribution error is
# correlated with the instrument and the first stage is partly measuring
# secondary migration rather than initial settlement.
if (!exists("shift_share_data")) {
  cat("\n  shift_share_data not in analysis_data.RData -- skipping step 2.\n")
} else {
  ss <- as.data.table(shift_share_data)
  # Total 1990 share held by each metro, summed over origins and weighted by
  # each origin's national size, which is what makes a metro "high-share" for
  # the inflows actually arriving.
  w1990 <- ss[, .(share = share[1]), by = .(MSA, BPLD)]
  natsz <- a[, .(nat = .N), by = BPLD]
  w1990 <- merge(w1990, natsz, by = "BPLD")
  msa_share <- w1990[, .(share_index = sum(share * nat) / sum(nat)), by = MSA]
  msa_share[, MSA := as.integer(MSA)]

  # THE DISCREPANCY MUST BE SCALE-FREE BEFORE IT IS CORRELATED WITH SIZE.
  # An absolute percentage-point gap is mechanically larger for a metro
  # holding 9 per cent of arrivals than for one holding 0.5 per cent, even
  # under proportionally identical error, so correlating it against a share
  # index measures size twice. Three quantities are reported instead:
  #
  #   gap_pp_yr   mean absolute gap per arrival year, in points -- the
  #               magnitude, comparable across years but still scale-bound
  #   gap_rel     that gap relative to the metro's OWN share -- scale-free,
  #               and the quantity to correlate against the share index
  #   gap_signed  signed, so direction is visible: positive means the pooled
  #               series credits the metro MORE than the near-arrival series
  #
  # The signed relative gap is the adverse quantity. If it is positive and
  # rising in the share index, the pooled series is systematically moving
  # arrivals toward metros with large 1990 co-national communities.
  nyr <- cells[, uniqueN(YRIMMIG)]
  dd <- cells[, .(sh_all_mean  = mean(sh_all),
                  gap_pp_yr    = 100 * mean(abs(sh_all - sh_near)),
                  gap_signed   = 100 * mean(sh_all - sh_near)), by = MSA]
  dd[, `:=`(gap_rel        = gap_pp_yr   / (100 * sh_all_mean),
            gap_rel_signed = gap_signed  / (100 * sh_all_mean))]
  dd <- merge(dd, msa_share, by = "MSA")

  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("2. DOES THE DISCREPANCY TRACK THE 1990 SHARES?\n\n")
  cat(sprintf("  metros compared: %d | arrival years averaged over: %d\n",
              nrow(dd), nyr))
  cat("\n  Correlation with the 1990 share index:\n")
  cat(sprintf("    absolute gap per year   (scale-bound) : %+.4f  [rank %+.4f]\n",
              dd[, cor(gap_pp_yr, share_index)],
              dd[, cor(gap_pp_yr, share_index, method = "spearman")]))
  cat(sprintf("    RELATIVE gap            (scale-free)  : %+.4f  [rank %+.4f]\n",
              dd[, cor(gap_rel, share_index)],
              dd[, cor(gap_rel, share_index, method = "spearman")]))
  cat(sprintf("    RELATIVE gap, signed    (scale-free)  : %+.4f  [rank %+.4f]\n",
              dd[, cor(gap_rel_signed, share_index)],
              dd[, cor(gap_rel_signed, share_index, method = "spearman")]))
  cat("\n  Read the scale-free rows. The first is reported only so the\n")
  cat("  difference between it and them is visible: a large correlation on the\n")
  cat("  absolute gap alongside a small one on the relative gap means the\n")
  cat("  apparent relationship was metro size, not attribution.\n")

  cat("\n  Share-index quartile, to see whether the gateways stand apart:\n")
  dd[, q := cut(share_index, quantile(share_index, 0:4/4), include.lowest = TRUE,
                labels = c("Q1 smallest", "Q2", "Q3", "Q4 largest"))]
  print(dd[, .(metros = .N,
               med_share_index = round(median(share_index), 5),
               med_gap_pp_yr   = round(median(gap_pp_yr), 4),
               med_gap_rel     = round(median(gap_rel), 4),
               med_gap_signed  = round(median(gap_rel_signed), 4)),
           by = q][order(q)], row.names = FALSE)
  cat("\n  A POSITIVE correlation is the adverse case: metros with large 1990\n")
  cat("  co-national communities would be gaining arrivals in the pooled\n")
  cat("  series that the near-arrival series does not credit to them, which\n")
  cat("  is secondary migration toward established communities and would mean\n")
  cat("  the first stage partly measures that movement.\n")
  cat("  A correlation near zero means the attribution error is orthogonal to\n")
  cat("  the instrument, which is the condition under which it is absorbed by\n")
  cat("  the first stage and does not bias two-stage least squares.\n")

  cat("\n  Ten metros with the largest RELATIVE signed gap (pooled over-credits):\n")
  print(head(dd[order(-gap_rel_signed),
                .(MSA, share_index = round(share_index, 5),
                  sh_all_pct = round(100 * sh_all_mean, 3),
                  gap_pp_yr = round(gap_pp_yr, 4),
                  gap_rel_signed = round(gap_rel_signed, 4))], 10),
        row.names = FALSE)
  cat("\n  Ten metros with the largest ABSOLUTE gap, for comparison:\n")
  print(head(dd[order(-gap_pp_yr),
                .(MSA, share_index = round(share_index, 5),
                  sh_all_pct = round(100 * sh_all_mean, 3),
                  gap_pp_yr = round(gap_pp_yr, 4),
                  gap_rel_signed = round(gap_rel_signed, 4))], 10),
        row.names = FALSE)
  cat("\n  If the two lists share their top entries, the concentration is real.\n")
  cat("  If the absolute list is the large metros and the relative list is not,\n")
  cat("  the concentration was size.\n")

  write.csv(dd, "outputs/tables/attribution_by_msa.csv", row.names = FALSE)
  cat("\n  Saved: outputs/tables/attribution_by_msa.csv\n")

  # =========================================================================
  # 3. THE TEST THAT ACTUALLY BINDS: WITHIN METRO, ACROSS ARRIVAL YEARS
  # =========================================================================
  # Step 2 asks whether high-share metros are more affected on average. That
  # is the CROSS-METRO dimension, and the specification carries MSA fixed
  # effects, which absorb it. If New York's attribution error is a constant
  # 5 per cent of its share in every year, the fixed effect removes it
  # entirely and it is harmless.
  #
  # What is not absorbed is variation WITHIN a metro across arrival years.
  # If a metro's error is larger in the years when its predicted inflow is
  # larger, that error moves with the identifying variation and survives the
  # fixed effects. That is the condition under which measurement error in the
  # endogenous regressor stops being benign.
  #
  # Both series are demeaned by metro, so what is correlated is deviation
  # from each metro's own average -- the dimension the design uses.
  pred <- as.data.table(msa_year_predicted_inflows)
  setnames(pred, "YEAR", "YRIMMIG")
  pred[, MSA := as.integer(MSA)]

  wi <- merge(cells[, .(MSA, YRIMMIG, gap = sh_all - sh_near)],
              pred[, .(MSA, YRIMMIG, z = total_predicted_inflow_nodiv)],
              by = c("MSA", "YRIMMIG"))
  # Demeaning removes each metro's LEVEL but not its VARIANCE, so a pooled
  # correlation on the demeaned series is still dominated by the metros with
  # the widest swings -- the same scale confound as step 2, one level down.
  # Standardising within metro removes it: every metro then contributes a
  # unit-variance series and the pooled figure is an equal-weighted average
  # of within-metro relationships rather than a size-weighted one.
  wi[, `:=`(gap_w = gap - mean(gap), z_w = z - mean(z)), by = MSA]
  wi[, `:=`(gap_s = if (sd(gap_w) > 0) gap_w / sd(gap_w) else NA_real_,
            z_s   = if (sd(z_w)   > 0) z_w   / sd(z_w)   else NA_real_), by = MSA]

  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("3. WITHIN METRO, ACROSS ARRIVAL YEARS -- THE BINDING TEST\n\n")
  cat(sprintf("  (MSA, year) cells matched to the instrument: %s of %s\n",
              format(nrow(wi), big.mark = ","), format(nrow(cells), big.mark = ",")))
  if (nrow(wi) < 0.9 * nrow(cells))
    cat("  *** Fewer than 90%% matched. Check the MSA/YEAR keys before reading on.\n")

  cat("\n  Pooled across cells:\n")
  cat(sprintf("    demeaned only  (size-weighted, NOT the one to read) : %+.4f\n",
              wi[, cor(gap_w, z_w)]))
  cat(sprintf("    STANDARDISED within metro (scale-free)             : %+.4f\n",
              wi[!is.na(gap_s) & !is.na(z_s), cor(gap_s, z_s)]))
  cat(sprintf("    rank correlation                                   : %+.4f\n",
              wi[, cor(gap_w, z_w, method = "spearman")]))
  cat("\n    The first is reported only so the gap between it and the second is\n")
  cat("    visible. Demeaning removes levels, not variances, so metros with the\n")
  cat("    widest swings still dominate it.\n")

  # Per metro, so a few large metros cannot carry the pooled figure.
  per <- wi[, .(n = .N, r = if (.N > 3 && sd(z_w) > 0 && sd(gap_w) > 0)
                              cor(gap_w, z_w) else NA_real_), by = MSA]
  per <- per[!is.na(r)]
  cat(sprintf("\n  Estimated separately for each of %d metros:\n", nrow(per)))
  cat(sprintf("    mean within-metro correlation   : %+.4f\n", per[, mean(r)]))
  cat(sprintf("    median                          : %+.4f\n", per[, median(r)]))
  cat(sprintf("    metros with r > +0.3            : %d (%.1f%%)\n",
              per[r > 0.3, .N], 100 * per[, mean(r > 0.3)]))
  cat(sprintf("    metros with r < -0.3            : %d (%.1f%%)\n",
              per[r < -0.3, .N], 100 * per[, mean(r < -0.3)]))

  # A tail count means nothing without the rate chance alone would produce.
  # Under a true correlation of zero, Fisher's z has standard error
  # 1/sqrt(n-3) with n the arrival years available to that metro, so the
  # expected share beyond +/-0.3 follows directly. Computed per metro and
  # averaged, because metros differ slightly in how many years they have.
  per <- merge(per, wi[, .(n_yrs = .N), by = MSA], by = "MSA")
  per[, se_z := 1 / sqrt(pmax(n_yrs - 3L, 1L))]
  exp_pos <- per[, mean(pnorm(atanh(0.3), 0, se_z, lower.tail = FALSE))]
  exp_neg <- per[, mean(pnorm(atanh(-0.3), 0, se_z))]
  cat(sprintf("\n    expected beyond +0.3 under a true zero : %.1f%%\n", 100 * exp_pos))
  cat(sprintf("    expected beyond -0.3 under a true zero : %.1f%%\n", 100 * exp_neg))
  cat(sprintf("    median arrival years per metro         : %.0f\n",
              per[, median(n_yrs)]))
  cat("\n    Compare the observed tails against these, not against zero. With\n")
  cat("    roughly a dozen years per metro a correlation of 0.3 is well inside\n")
  cat("    sampling noise, so a tail count near the expected rate is evidence\n")
  cat("    of no relationship rather than of a weak one.\n")
  cat("\n  If these cluster near zero with the positive and negative tails\n")
  cat("  roughly balanced, the attribution error is noise with respect to the\n")
  cat("  identifying variation and the first stage absorbs it. A mean well\n")
  cat("  above zero would mean the error rises with predicted inflow within\n")
  cat("  metros, which the MSA fixed effects cannot remove.\n")

  write.csv(wi, "outputs/tables/attribution_within_msa.csv", row.names = FALSE)
  cat("\n  Saved: outputs/tables/attribution_within_msa.csv\n")
}

# --- Summary for the chapter table ----------------------------------------
# The three steps produce one number each that the chapter actually reports.
# Assembled here rather than recomputed in the table script, so the figures in
# the text cannot drift from the ones printed above.
summ <- data.frame(
  quantity = c(
    "Correlation of within-year metro shares, near vs all",
    "Mean absolute discrepancy per cell (percentage points)",
    "Correlation of signed relative gap with 1990 share index",
    "Correlation of gap with instrument, within metro (standardised)",
    "Mean of per-metro correlations",
    "Metros with correlation beyond +0.3 (per cent)",
    "Expected beyond +0.3 under a true zero (per cent)"),
  value = c(
    cells[, cor(sh_all, sh_near)],
    100 * cells[, mean(abs(sh_all - sh_near))],
    if (exists("dd")) dd[, cor(gap_rel_signed, share_index)] else NA_real_,
    if (exists("wi")) wi[!is.na(gap_s) & !is.na(z_s), cor(gap_s, z_s)] else NA_real_,
    if (exists("per")) per[, mean(r)] else NA_real_,
    if (exists("per")) 100 * per[, mean(r > 0.3)] else NA_real_,
    if (exists("exp_pos")) 100 * exp_pos else NA_real_),
  stringsAsFactors = FALSE)
summ$arrival_years <- paste0(YR_LO, "-", YR_HI)
summ$recent_max    <- RECENT_MAX
write.csv(summ, "outputs/tables/attribution_summary.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/attribution_summary.csv\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(cells, "outputs/tables/attribution_stability.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/attribution_stability.csv\n")

cat("\n", strrep("-", 72), "\n", sep = "")
cat("READING THIS\n\n")
cat("This is a descriptive check and re-estimates nothing, deliberately. The\n")
cat("recency rule cannot be applied to the full series without removing whole\n")
cat("arrival years, so the estimation-based version of this test is not\n")
cat("available on this survey calendar.\n\n")
cat("Step 1 establishes the magnitude: whether measuring arrivals close to\n")
cat("arrival changes where they appear to have gone at all.\n\n")
cat("Step 2 asks whether the error is larger in high-share metros. It is the\n")
cat("weaker of the two correlation tests, because the specification carries\n")
cat("MSA fixed effects, which absorb anything that differs across metros but\n")
cat("is constant within them.\n\n")
cat("STEP 3 IS THE ONE THAT BINDS. Identification comes from variation within\n")
cat("a metro across cohorts, so what matters is whether the attribution error\n")
cat("moves with the instrument along that dimension. Near-zero there is what\n")
cat("licenses the claim that the error is absorbed by the first stage and\n")
cat("leaves two-stage least squares consistent.\n\n")
cat("The 2006-2019 window is more favourable to the design than the full\n")
cat("sample: these are the years with the densest coverage and the shortest\n")
cat("gap between arrival and enumeration. A clean result here does not\n")
cat("license the same claim for 1990s arrivals, and the text should say so.\n")
cat("19_attribution_stability_check.R complete.\n")
