# =============================================================================
# 09b_exposure_distribution.R  --  What does the treatment variable actually
#                                  look like, and where does its variation sit?
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Added 2026-10-30.
#
# WHY THIS EXISTS
# Every headline magnitude in the chapter is reported per standard deviation of
# exposure, because the exposure measure is a count of arrivals and carries no
# natural unit of its own. That SD is therefore the scaling constant behind the
# 4.74 percentage points in Table 1.1, behind the per-SD comparison to
# Hunt (2017), and behind every per-SD figure in the robustness section -- and
# it appears nowhere in the paper. A reader cannot currently find out what a
# one-standard-deviation increase in exposure is.
#
# WHAT IT REPORTS
#   Panel A  the distribution of the instrument and of the realized inflow,
#            at the individual level among identifying observations
#   Panel B  where the variation lives: how much survives removing metropolitan
#            -area means, and how much survives removing birth-cohort means.
#            Treatment is defined at MSA x birth cohort, so both margins carry
#            variation and the split is not obvious a priori.
#   Panel C  the counts that fix the design's dimensions
#
# WHY THE IDENTIFYING SUBSAMPLE
# exposure_missing == 1 means z is exactly zero by construction (the ages-15-25
# window opens before the 1990 base year). Those observations are retained in
# estimation -- they inform the fixed effects and control coefficients -- but
# including them in a distribution of the treatment would describe a spike at
# zero rather than the variation the estimate uses. sd_z in 09_main_results.R
# is taken on the identifying subsample for the same reason, so the SD reported
# here is the same number that scales the published coefficients. The pooled
# SD is reported alongside it so the difference is visible rather than implied.
#
# Writes outputs/tables/exposure_distribution.csv
#    and outputs/tables/table_exposure_distribution.tex
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")
source("scripts/analysis/00b_table_helper.R")

library(data.table)
library(fixest)

stopifnot(exists("write_dissertation_table"))
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
all_natives_175 <- as.data.table(.tmp_env$all_natives_175)
rm(.tmp_env); gc()

.need <- c("avg_exposure_nodiv", "avg_actual_exposure_nodiv",
           "exposure_missing", "AGE", "YEAR", "MSA")
.miss <- setdiff(.need, names(all_natives_175))
if (length(.miss))
  stop("analysis_data.RData is stale or incomplete -- missing: ",
       paste(.miss, collapse = ", "),
       "\nRe-run 06_shift_share_iv.R and 07_exposure_iv.R first.")
rm(.need, .miss)

# Same sample construction as 09_main_results.R and 10_exposure_missing_check.R.
d <- all_natives_175[MSA != 0]
d[, `:=`(
  birth_year = YEAR - AGE,
  z = avg_exposure_nodiv,
  x = avg_actual_exposure_nodiv
)]

id <- d[exposure_missing == 0]

# --- Panel A: distribution ---------------------------------------------------
qs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
dist_stats <- function(v) {
  q <- quantile(v, probs = qs, na.rm = TRUE, names = FALSE)
  c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
    p10 = q[1], p25 = q[2], p50 = q[3], p75 = q[4], p90 = q[5])
}
a_z <- dist_stats(id$z)
a_x <- dist_stats(id$x)

# --- Panel B: where the variation lives --------------------------------------
# Residual SD after sweeping out one margin's means. Sweeping MSA means leaves
# the across-cohort variation within metros; sweeping cohort means leaves the
# across-metro variation within cohorts. Variance share is 1 - resid/total, so
# it answers "how much of the spread is this margin accounting for".
within_sd <- function(v, g) {
  t <- data.table(v = v, g = g)
  t[, r := v - mean(v), by = g]
  sd(t$r, na.rm = TRUE)
}
decomp <- function(v, msa, coh) {
  s  <- sd(v, na.rm = TRUE)
  wm <- within_sd(v, msa)   # within MSA  -> across cohorts
  wc <- within_sd(v, coh)   # within cohort -> across MSAs
  c(sd_within_msa = wm,
    sd_within_cohort = wc,
    vshare_msa = 100 * (1 - (wm^2) / (s^2)),
    vshare_cohort = 100 * (1 - (wc^2) / (s^2)))
}
b_z <- decomp(id$z, id$MSA, id$birth_year)
b_x <- decomp(id$x, id$MSA, id$birth_year)

# --- Panel B2: what survives the fixed effects -------------------------------
# The decomposition above sweeps one margin at a time. The estimating equation
# sweeps metropolitan-area and survey-year effects JOINTLY, so the variation
# actually available to identify the coefficient is what remains after both.
# That residual SD is the quantity a referee asks for when the headline is
# reported per standard deviation: 72 per cent of the instrument's raw variance
# is between metros, and the MSA effects absorb all of it.
#
# Fixed effects only, no controls. The controls are individual-level
# demographics (sex, marital status, poverty, parental education), and the
# treatment is constant within an MSA x birth-cohort cell, so they can absorb
# treatment variation only through their correlation with cell means -- a
# second-order channel. The fixed effects are what remove variation directly.
#
# Means are taken over the identifying observations, consistent with the rest
# of the table. The pooled-sample figure is printed below for comparison; it is
# not the relevant one, because the non-identifying observations carry z = 0 by
# construction and so contribute no variation of their own.
dm_id <- as.data.table(demean(X = id[, .(z, x)], f = id[, .(YEAR, MSA)]))
sd_fe <- c(z = sd(dm_id$z, na.rm = TRUE), x = sd(dm_id$x, na.rm = TRUE))

dm_pool <- as.data.table(demean(X = d[, .(z, x)], f = d[, .(YEAR, MSA)]))
sd_fe_pooled <- c(z = sd(dm_pool$z, na.rm = TRUE), x = sd(dm_pool$x, na.rm = TRUE))
rm(dm_pool); gc()

# --- Panel C: dimensions -----------------------------------------------------
n_cells  <- id[, uniqueN(paste(MSA, birth_year))]
n_msa    <- id[, uniqueN(MSA)]
n_cohort <- id[, uniqueN(birth_year)]
sd_pooled_z <- sd(d$z, na.rm = TRUE)

# --- Console report ----------------------------------------------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("EXPOSURE DISTRIBUTION\n\n")
cat(sprintf("  observations (pooled)        %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  observations (identifying)   %s\n", format(nrow(id), big.mark = ",")))
cat(sprintf("  MSA x birth-cohort cells     %s\n", format(n_cells, big.mark = ",")))
cat(sprintf("  metropolitan areas           %d\n", n_msa))
cat(sprintf("  birth cohorts                %d  (%d-%d)\n",
            n_cohort, id[, min(birth_year)], id[, max(birth_year)]))
cat(sprintf("\n  SD of instrument, identifying  %12.1f   <- scales every per-SD figure\n", a_z[["sd"]]))
cat(sprintf("  SD of instrument, pooled       %12.1f\n", sd_pooled_z))
cat(sprintf("  mean of instrument             %12.1f\n", a_z[["mean"]]))
cat(sprintf("  mean of realized inflow        %12.1f\n", a_x[["mean"]]))
cat(sprintf("  SD of realized inflow          %12.1f\n", a_x[["sd"]]))
cat(sprintf("\n  instrument SD within MSA       %12.1f  (%.1f%% of variance from MSA means)\n",
            b_z[["sd_within_msa"]], b_z[["vshare_msa"]]))
cat(sprintf("  instrument SD within cohort    %12.1f  (%.1f%% of variance from cohort means)\n",
            b_z[["sd_within_cohort"]], b_z[["vshare_cohort"]]))

cat("\n  NET OF MSA AND YEAR EFFECTS JOINTLY (the identifying variation)\n")
cat(sprintf("    instrument, identifying      %12.1f   (%.1f%% of its raw SD)\n",
            sd_fe[["z"]], 100 * sd_fe[["z"]] / a_z[["sd"]]))
cat(sprintf("    instrument, pooled           %12.1f\n", sd_fe_pooled[["z"]]))
cat(sprintf("    realized inflow, identifying %12.1f   (%.1f%% of its raw SD)\n",
            sd_fe[["x"]], 100 * sd_fe[["x"]] / a_x[["sd"]]))

# Rescaling the published coefficients onto the identifying SD. Read these
# against the per-SD figures in main_results.csv, which use the raw SD.
.rf   <- 7.24926434428727e-06   # headline reduced form, Table 1.1 column 1
.tsls <- 1.11957825376898e-05   # headline 2SLS,         Table 1.1 column 1
cat("\n  IMPLIED RESCALING OF THE HEADLINE (high-school margin)\n")
cat(sprintf("    reduced form: %5.2f pp per raw SD  ->  %5.2f pp per identifying SD\n",
            100 * .rf * a_z[["sd"]], 100 * .rf * sd_fe[["z"]]))
cat(sprintf("    2SLS:         %5.2f pp per raw SD  ->  %5.2f pp per identifying SD\n",
            100 * .tsls * a_x[["sd"]], 100 * .tsls * sd_fe[["x"]]))
cat("    The raw-SD figures are the ones published; these are the comparison.\n")

# --- CSV ---------------------------------------------------------------------
res <- data.table(
  quantity = c("mean", "sd", "p10", "p25", "p50", "p75", "p90",
               "sd_within_msa", "sd_within_cohort",
               "vshare_msa_pct", "vshare_cohort_pct",
               "sd_net_msa_year", "sd_net_msa_year_pooled",
               "n_identifying", "n_pooled", "sd_pooled",
               "n_cells", "n_msa", "n_cohort"),
  instrument = c(a_z, b_z, sd_fe[["z"]], sd_fe_pooled[["z"]],
                 nrow(id), nrow(d), sd_pooled_z,
                 n_cells, n_msa, n_cohort),
  realized   = c(a_x, b_x, sd_fe[["x"]], sd_fe_pooled[["x"]],
                 nrow(id), nrow(d), NA_real_,
                 n_cells, n_msa, n_cohort)
)
fwrite(res, "outputs/tables/exposure_distribution.csv")
cat("\nWrote outputs/tables/exposure_distribution.csv\n")

# --- LaTeX -------------------------------------------------------------------
f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
fn <- function(x) formatC(x, format = "d", big.mark = ",")

tbl <- data.frame(
  Statistic = c(
    "Mean", "Standard deviation",
    "10th percentile", "25th percentile", "Median",
    "75th percentile", "90th percentile",
    "@MIDRULE@",
    "SD within metropolitan area (across cohorts)",
    "SD within birth cohort (across metropolitan areas)",
    "Variance share, metropolitan-area means (per cent)",
    "Variance share, birth-cohort means (per cent)",
    "SD net of metropolitan-area and survey-year effects",
    "@MIDRULE@",
    "Identifying observations",
    "MSA $\\times$ birth-cohort cells",
    "Metropolitan areas",
    "Birth cohorts"),
  Instrument = c(
    f0(a_z[["mean"]]), f0(a_z[["sd"]]),
    f0(a_z[["p10"]]), f0(a_z[["p25"]]), f0(a_z[["p50"]]),
    f0(a_z[["p75"]]), f0(a_z[["p90"]]),
    "",
    f0(b_z[["sd_within_msa"]]), f0(b_z[["sd_within_cohort"]]),
    f1(b_z[["vshare_msa"]]), f1(b_z[["vshare_cohort"]]),
    f0(sd_fe[["z"]]),
    "",
    fn(nrow(id)), fn(n_cells), fn(n_msa), fn(n_cohort)),
  `Realized inflow` = c(
    f0(a_x[["mean"]]), f0(a_x[["sd"]]),
    f0(a_x[["p10"]]), f0(a_x[["p25"]]), f0(a_x[["p50"]]),
    f0(a_x[["p75"]]), f0(a_x[["p90"]]),
    "",
    f0(b_x[["sd_within_msa"]]), f0(b_x[["sd_within_cohort"]]),
    f1(b_x[["vshare_msa"]]), f1(b_x[["vshare_cohort"]]),
    f0(sd_fe[["x"]]),
    "",
    fn(nrow(id)), fn(n_cells), fn(n_msa), fn(n_cohort)),
  check.names = FALSE, stringsAsFactors = FALSE
)

write_dissertation_table(
  tbl,
  caption = "Distribution of Immigrant Exposure and of the Realized Inflow",
  label = "tab:exposure_distribution",
  notes = paste(
    "The instrument is the 1990 origin-group settlement share interacted with",
    "national arrival counts by origin, summed across origins and averaged over",
    "each individual's ages 15--25 window. The endogenous variable is the",
    "realized inflow into the metropolitan area itself, carrying no settlement",
    "share, averaged over the same window. Both are counts of arrivals, which is",
    "why effects are reported per standard deviation rather than per unit.",
    "Statistics are over the identifying observations, those whose window opens",
    "in or after the 1990 base year and so can carry non-zero exposure. The",
    "middle block reports how much spread survives removing each margin's means;",
    "the variance shares do not sum to 100 because the two margins are not",
    "orthogonal in an unbalanced panel. Its final row sweeps metropolitan-area",
    "and survey-year effects jointly, giving the variation available to identify",
    "the coefficients."),
  align = "lcc",
  file = "outputs/tables/table_exposure_distribution.tex"
)

cat("\n", strrep("-", 72), "\n", sep = "")
