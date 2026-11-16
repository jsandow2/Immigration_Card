# =============================================================================
# 16c_collapse_diagnostic.R  --  which cell definition reproduces the headline
#
# 16b_akm_standard_errors.R has to run the design at the level where the
# shift-share structure lives, because that is the only level AKM's variance
# estimator is defined on. That means collapsing 3.5 million records to cells,
# and the collapse is only trustworthy if it reproduces the individual-level
# estimate. So far it does not: MSA x cohort cells return the wrong sign.
#
# This script isolates the cause without building the share matrix, which is
# the expensive part of 16b. It runs the same 2SLS at the individual level and
# at two cell definitions, with and without controls.
#
# WHAT THE COMPARISONS SHOW
#   A vs D   Both z and x are constant within MSA x cohort, so a cell
#            definition at least that fine, carrying the same fixed effects,
#            must reproduce the individual estimate EXACTLY rather than
#            approximately. A mismatch means the cell definition has dropped a
#            dimension the estimate depends on.
#   C vs D   Whether survey year has to be a fixed effect in its own right.
#            MSA x cohort cells cannot carry one: averaging over survey years
#            is what defines the cell.
#   D vs E   Whether aggregating the controls to cell shares costs anything
#            once the fixed-effect structure is right.
#
# Writes nothing. Run it, read the table, then fix 16b.
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

d <- all_natives_175[exposure_missing == 0 & MSA %in% selected_msas]
d[, birth_year := YEAR - AGE]
d[, z := avg_exposure_nodiv]
d[, x := avg_actual_exposure_nodiv]
d[, y := as.integer(Educ_Group != "non_HS")]

# Both sides of the first stage are MSA x cohort objects. If either varies
# within a cell, a cell-mean regression is not algebraically equivalent to the
# individual one and no comparison below is meaningful.
chk <- d[, .(nz = uniqueN(round(z, 6)), nx = uniqueN(round(x, 6))),
         by = .(MSA, birth_year)]
cat(sprintf("\n  max distinct z within MSA x cohort: %d\n", chk[, max(nz)]))
cat(sprintf("  max distinct x within MSA x cohort: %d\n", chk[, max(nx)]))
if (chk[, max(nz)] > 1L || chk[, max(nx)] > 1L)
  stop("z or x varies within MSA x cohort; the collapse cannot be lossless.")
rm(chk)

CTRL <- "Black + SEX + Marital + ChildAgeGroup + POVERTY + Mom_Educ + Dad_Educ"

fmt <- function(tag, m, note = "") {
  ct <- m$coeftable["fit_x", ]
  cat(sprintf("  %-4s %-34s %+.6e  SE %.3e  t %6.2f  %s\n",
              tag, note, ct[["Estimate"]], ct[["Std. Error"]],
              ct[["t value"]], ""))
  invisible(ct[["Estimate"]])
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("INDIVIDUAL LEVEL\n")
cat(strrep("=", 78), "\n")
A <- fmt("A", feols(as.formula(sprintf("y ~ %s | YEAR + MSA | x ~ z", CTRL)),
                    data = d, vcov = ~MSA), "controls, YEAR + MSA FE")
B <- fmt("B", feols(y ~ 1 | YEAR + MSA | x ~ z, data = d, vcov = ~MSA),
         "no controls, YEAR + MSA FE")

# --- cell definitions --------------------------------------------------------
cellA <- d[, .(n = .N, y = mean(y), x = mean(x), z = mean(z)),
           by = .(MSA, birth_year)]
cellB <- d[, .(n = .N, y = mean(y), x = mean(x), z = mean(z)),
           by = .(MSA, birth_year, YEAR)]

cat("\n  MSA x cohort cells:              ", format(nrow(cellA), big.mark = ","),
    "\n  MSA x cohort x survey-year cells:", format(nrow(cellB), big.mark = ","),
    "\n", sep = "")

cat("\n", strrep("=", 78), "\n", sep = "")
cat("CELL LEVEL, NO CONTROLS\n")
cat(strrep("=", 78), "\n")
C <- fmt("C", feols(y ~ 1 | MSA | x ~ z, data = cellA, weights = ~n,
                    vcov = ~MSA), "MSA x cohort, MSA FE only")
D <- fmt("D", feols(y ~ 1 | MSA + YEAR | x ~ z, data = cellB, weights = ~n,
                    vcov = ~MSA), "MSA x cohort x year, MSA+YEAR FE")

cat("\n", strrep("=", 78), "\n", sep = "")
cat("VERDICT\n")
cat(strrep("=", 78), "\n")
cat(sprintf("  B (individual, no controls)        %+.6e\n", B))
cat(sprintf("  D (3-way cells, same FE)           %+.6e   ratio %.6f\n",
            D, D / B))
cat(sprintf("  C (2-way cells, MSA FE only)       %+.6e   ratio %.6f\n",
            C, C / B))
cat("\n")
if (abs(D / B - 1) < 1e-6) {
  cat("  D reproduces B. The collapse is lossless once survey year is carried\n")
  cat("  as a fixed effect, so 16b must use MSA x cohort x survey-year cells.\n")
} else {
  cat("  D does NOT reproduce B. The cell definition is not the whole story --\n")
  cat("  do not proceed to the share matrix until this is understood.\n")
}
cat(sprintf("\n  A (individual, with controls)      %+.6e\n", A))
cat("  Compare A with D to see what the controls are worth once the fixed\n")
cat("  effects are right.\n")
