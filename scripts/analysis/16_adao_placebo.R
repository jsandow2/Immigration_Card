# =============================================================================
# 16_adao_placebo.R  --  Adao, Kolesar & Morales (2019) placebo test
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 19_adao_placebo.R.
# Renumbered 2026-10-21 from 20_adao_placebo.R (archive-and-renumber pass).
# REBUILT 2026-09-29 on the final specification.
#
# WHAT THIS TESTS, AND WHY IT IS CHAPTER 2's INFERENCE
# Adao, Kolesar & Morales show that shift-share designs induce correlation in
# residuals ACROSS regions with similar share vectors -- regions exposed to the
# same national shocks move together even when they are far apart. Conventional
# region-clustered standard errors do not capture that, so they can be
# understated. This is the canonical robustness check for the design and the
# one a referee in this literature expects. Alternative inference procedures
# are taken up in the later chapter.
#
# Two things are computed.
#
#   1. PLACEBO IV. The same specification, run on outcomes determined at birth
#      and therefore incapable of responding to immigration during adolescence:
#      Black, Asian, Aboriginal. A shift-share instrument that "predicts" these
#      is picking up something other than the exposure it claims to measure.
#      VETSTAT is deliberately excluded -- military service is post-exposure and
#      not a clean placebo.
#
#      Read the MAGNITUDE, not just the p-value. The script reports each placebo
#      effect in the same percentage-points-per-SD units as the headline and as
#      a percentage of it. A placebo effect that is a tenth the size of the real
#      one is reassuring whatever its t-statistic; one that is half the size is
#      damning whatever its t-statistic.
#
#   2. MSA-VERSUS-STATE CLUSTERING. If cross-regional correlation is inflating
#      precision, clustering at a coarser level should widen the standard error.
#      A state/MSA SE ratio near or below one means MSA clustering is not
#      understating uncertainty on this dimension; a ratio well above one means
#      it is.
#
# WHAT CHANGED IN THE REBUILD
#   - Endogenous variable is avg_actual_exposure_nodiv and the instrument is
#     avg_exposure_nodiv. Previously total_actual_inflow (survey year) and
#     avg_exposure (wave-divided) -- both superseded; see the Step 4b note in
#     07_exposure_iv.R.
#   - Converted from AER::ivreg + sandwich::vcovCL to fixest::feols, which
#     absorbs the fixed effects and clusters natively.
#   - Now SELF-SUFFICIENT. It previously inherited `samp` from Section 19's
#     environment and read pooled objects from after_sec11.RData. Section 11 is
#     the chunk-loop pooling, which is retired, and the after_sec19 checkpoint
#     carries pre-respecification objects. This script now refreshes from
#     analysis_data.RData and rebuilds the sample under the same seed, so it
#     reproduces Section 19's draw without depending on session state.
#   - The dead Section 11 cross-check is replaced by a direct MSA-versus-state
#     comparison estimated on the full sample under the final specification.
#
# Writes outputs/tables/adao_placebo_results.csv.
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(fixest)

cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
all_natives_175 <- as.data.table(.tmp_env$all_natives_175)
selected_msas   <- .tmp_env$selected_msas
rm(.tmp_env); gc()

.need <- c("avg_exposure_nodiv", "avg_actual_exposure_nodiv", "exposure_missing",
           "Educ_Group", "MSA", "STATEFIP")
.miss <- setdiff(.need, names(all_natives_175))
if (length(.miss))
  stop("analysis_data.RData is stale or incomplete -- missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R first.")
rm(.need, .miss)

cat(sprintf("\n%s\n", strrep("=", 60)))
cat("Section 20: Adao et al. (2019) Placebo Test\n")
cat(sprintf("%s\n\n", strrep("=", 60)))

# --- Outcome and placebo variables, built BEFORE the sample is drawn ---------
# Ordering matters: samp is a copy, so anything added to all_natives_175
# afterwards does not reach it. This is the failure that broke 19 twice.
all_natives_175[, y_HSplus := as.integer(Educ_Group != "non_HS")]

PLACEBOS <- c("Black", "Asian", "Aboriginal")
avail    <- PLACEBOS[PLACEBOS %in% names(all_natives_175)]
if (!length(avail))
  stop("None of the placebo outcomes are present: ", paste(PLACEBOS, collapse = ", "))
if (length(avail) < length(PLACEBOS))
  cat(sprintf("NOTE: placebo outcomes not found and skipped: %s\n",
              paste(setdiff(PLACEBOS, avail), collapse = ", ")))

# --- Sample: same seed and construction as Section 19 ------------------------
set.seed(42)
n_sample   <- 500000
sample_idx <- sample(nrow(all_natives_175), n_sample)
samp       <- droplevels(all_natives_175[sample_idx])

cat(sprintf("Sample: %s observations | %d MSAs | %s identifying\n",
            format(nrow(samp), big.mark = ","), uniqueN(samp$MSA),
            format(samp[exposure_missing == 0, .N], big.mark = ",")))
sd_z <- samp[exposure_missing == 0, sd(avg_exposure_nodiv)]
cat(sprintf("SD of instrument among identifying observations: %.1f\n\n", sd_z))

# --- Estimation ---------------------------------------------------------------
# Control set excludes the race indicators throughout, so that it is identical
# across the placebos and across the reference outcome -- otherwise the Asian
# placebo would carry Black as a control while the Black placebo could not, and
# the magnitudes would not be comparable.
CTRL_PLA <- "SEX + Marital + ChildAgeGroup + POVERTY +
             Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

run_iv <- function(outcome, dat, cluster_var = "MSA") {
  if (uniqueN(dat[[outcome]], na.rm = TRUE) < 2) return(NULL)
  f <- as.formula(sprintf(
    "%s ~ %s | avg_actual_exposure_nodiv ~ avg_exposure_nodiv",
    outcome, CTRL_PLA))
  m <- tryCatch(feols(f, data = dat,
                      cluster = as.formula(paste0("~", cluster_var))),
                error = function(e) { cat(sprintf("FAILED: %s\n",
                                                  conditionMessage(e))); NULL })
  if (is.null(m)) return(NULL)
  m$coeftable["fit_avg_actual_exposure_nodiv", ]
}

# Reference: the real outcome on the same control set, so placebo magnitudes
# can be expressed relative to something.
cat("--- Reference outcome (same control set) ---\n")
ref <- run_iv("y_HSplus", samp)
ref_pp <- ref[["Estimate"]] * sd_z * 100
cat(sprintf("  %-12s est %+.5e | t %+7.3f | p %.4f | %+.4f pp/SD\n\n",
            "y_HSplus", ref[["Estimate"]], ref[["t value"]],
            ref[["Pr(>|t|)"]], ref_pp))

cat("--- Placebo IV regressions ---\n\n")
rows <- list()
for (y in avail) {
  ct_msa   <- run_iv(y, samp, "MSA")
  ct_state <- run_iv(y, samp, "STATEFIP")
  if (is.null(ct_msa) || is.null(ct_state)) {
    cat(sprintf("  %-12s SKIPPED\n", y)); next
  }
  pp    <- ct_msa[["Estimate"]] * sd_z * 100
  ratio <- ct_state[["Std. Error"]] / ct_msa[["Std. Error"]]
  cat(sprintf("  %-12s est %+.5e | t(MSA) %+7.3f (p %.4f) | t(State) %+7.3f | SE ratio %.3f\n",
              y, ct_msa[["Estimate"]], ct_msa[["t value"]],
              ct_msa[["Pr(>|t|)"]], ct_state[["t value"]], ratio))
  cat(sprintf("  %-12s %+.4f pp/SD  =  %.1f%% of the reference effect\n",
              "", pp, 100 * abs(pp / ref_pp)))
  rows[[length(rows) + 1]] <- data.frame(
    outcome = y, estimate = ct_msa[["Estimate"]],
    msa_se = ct_msa[["Std. Error"]], msa_t = ct_msa[["t value"]],
    msa_p = ct_msa[["Pr(>|t|)"]],
    state_se = ct_state[["Std. Error"]], state_t = ct_state[["t value"]],
    state_p = ct_state[["Pr(>|t|)"]], se_ratio = ratio,
    pp_per_sd = pp, pct_of_reference = 100 * abs(pp / ref_pp),
    stringsAsFactors = FALSE)
}
placebo_df <- rbindlist(rows)

cat(sprintf("\n%s\n", strrep("=", 60)))
cat("PLACEBO TEST RESULTS\n")
cat(sprintf("%s\n\n", strrep("=", 60)))
cat(sprintf("Placebo outcomes tested:        %d\n", nrow(placebo_df)))
cat(sprintf("Rejections at 5%% (MSA cluster):  %d of %d\n",
            placebo_df[, sum(msa_p < 0.05)], nrow(placebo_df)))
cat(sprintf("Rejections at 5%% (State cluster):%d of %d\n",
            placebo_df[, sum(state_p < 0.05)], nrow(placebo_df)))
cat(sprintf("Mean |t| (MSA):                 %.3f\n", placebo_df[, mean(abs(msa_t))]))
cat(sprintf("Mean |t| (State):               %.3f\n", placebo_df[, mean(abs(state_t))]))
cat(sprintf("Mean SE ratio (State/MSA):      %.3f\n", placebo_df[, mean(se_ratio)]))
cat(sprintf("Largest placebo as %% of reference effect: %.1f%%\n\n",
            placebo_df[, max(pct_of_reference)]))
print(placebo_df[, .(outcome, estimate, msa_t, msa_p, state_t, se_ratio,
                     pp_per_sd, pct_of_reference)],
      row.names = FALSE, digits = 3)

# --- MSA versus state clustering on the FULL sample, final specification -----
# Replaces a block that read pooled objects from after_sec11.RData. Section 11
# is the chunk-loop pooling, now retired, and its figures were on the
# superseded construction. Estimated directly instead -- linear models on 8.5M
# rows take seconds.
cat("\n", strrep("-", 60), "\n", sep = "")
cat("MSA vs STATE clustering, full sample, final specification\n\n")
d_full <- all_natives_175[MSA != 0 & MSA %in% selected_msas]
d_full <- droplevels(d_full)

clustering_rows <- list()
for (oc in c("y_HSplus")) {
  m_msa <- feols(as.formula(sprintf(
    "%s ~ %s | avg_actual_exposure_nodiv ~ avg_exposure_nodiv", oc, CTRL_PLA)),
    data = d_full, cluster = ~MSA)
  m_st  <- feols(as.formula(sprintf(
    "%s ~ %s | avg_actual_exposure_nodiv ~ avg_exposure_nodiv", oc, CTRL_PLA)),
    data = d_full, cluster = ~STATEFIP)
  e  <- m_msa$coeftable["fit_avg_actual_exposure_nodiv", ]
  s  <- m_st$coeftable["fit_avg_actual_exposure_nodiv", ]
  r  <- s[["Std. Error"]] / e[["Std. Error"]]
  # Saved as well as printed, added 2026-10-05. These figures were console-only,
  # so the chapter could not report the state-clustered standard error beside
  # the MSA-clustered one and the ratio had to be asserted rather than shown.
  clustering_rows[[length(clustering_rows) + 1L]] <- data.frame(
    outcome        = oc,
    n              = nrow(d_full),
    estimate       = e[["Estimate"]],
    msa_se         = e[["Std. Error"]],
    msa_t          = e[["t value"]],
    msa_p          = e[["Pr(>|t|)"]],
    msa_clusters   = uniqueN(d_full$MSA),
    state_se       = s[["Std. Error"]],
    state_t        = s[["t value"]],
    state_p        = s[["Pr(>|t|)"]],
    state_clusters = uniqueN(d_full$STATEFIP),
    se_ratio       = r,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-10s est %+.5e\n", oc, e[["Estimate"]]))
  cat(sprintf("             SE(MSA)   %.5e  t %7.3f  (%d clusters)\n",
              e[["Std. Error"]], e[["t value"]], uniqueN(d_full$MSA)))
  cat(sprintf("             SE(State) %.5e  t %7.3f  (%d clusters)\n",
              s[["Std. Error"]], s[["t value"]], uniqueN(d_full$STATEFIP)))
  cat(sprintf("             SE ratio (State/MSA) = %.3f\n", r))
  cat("             Interpret narrowly. A ratio at or below one shows no sign\n")
  cat("             of state-level cross-regional correlation inflating\n")
  cat("             precision, but it is NOT a clean test of it: there are far\n")
  cat("             fewer state than MSA clusters, and cluster-robust variance\n")
  cat("             is more downward-biased the fewer clusters there are. So a\n")
  cat("             smaller state SE is partly what the cluster count alone\n")
  cat("             would produce, and the comparison conflates that with the\n")
  cat("             correlation structure Adao et al. are concerned with.\n\n")
}

cat(strrep("-", 60), "\n", sep = "")
cat("READING THIS\n\n")
cat("The placebo outcomes are fixed at birth. Any apparent effect on them is\n")
cat("the instrument correlating with something other than formative-years\n")
cat("exposure. Judge it on pct_of_reference first: a placebo a tenth the size\n")
cat("of the real effect is reassuring whatever its t-statistic, and one half\n")
cat("the size is damning whatever its t-statistic. With 175 nominal clusters\n")
cat("the t-statistics themselves are not to be taken at face value.\n\n")
cat("On the clustering comparison, the defensible claim is narrow: a state/MSA\n")
cat("ratio at or below one shows no EVIDENCE of state-level cross-regional\n")
cat("correlation inflating precision. It does not establish the absence of it.\n")
cat("Coarsening the cluster definition also cuts the cluster count -- 45 states\n")
cat("against 175 metros here -- and cluster-robust variance is more\n")
cat("downward-biased the fewer clusters there are, so part of any reduction is\n")
cat("attributable to the cluster count rather than to the correlation\n")
cat("structure. A ratio well ABOVE one would be the informative case, since\n")
cat("the cluster-count effect pushes the other way.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(placebo_df, "outputs/tables/adao_placebo_results.csv", row.names = FALSE)
cat(sprintf("\nSaved: outputs/tables/adao_placebo_results.csv\n"))

# Full-sample clustering comparison on the final specification. Separate file
# from the placebo results because it is a different sample (all 8.5M rows,
# not the draw the placebos run on) and a different quantity (the headline
# estimate, not a placebo), and conflating the two in one table would invite
# exactly that misreading.
clustering_df <- rbindlist(clustering_rows)
write.csv(clustering_df, "outputs/tables/clustering_comparison.csv",
          row.names = FALSE)
cat("Saved: outputs/tables/clustering_comparison.csv\n")

cat(sprintf("\n%s\n", strrep("=", 60)))
cat("Section 20 complete.\n")
cat(sprintf("%s\n", strrep("=", 60)))
# NO CHECKPOINT WRITTEN, as of 2026-10-04.
#
# This wrote after_sec20.RData, which nothing reads. The checkpoint chain was
# dismantled between 2026-09-30 and 2026-10-04 -- every surviving script
# loads analysis_data.RData directly -- and after that the save was writing
# large objects to disk on every run for no consumer. Worse, a stale
# checkpoint left on disk is precisely how 15_rotemberg.R came to decompose
# a superseded estimate in September. The only checkpoint still in the
# Chapter 2 pipeline is after_sec22, written by 21 and read by 24.

cat("\n16_adao_placebo.R complete (no checkpoint).\n")
