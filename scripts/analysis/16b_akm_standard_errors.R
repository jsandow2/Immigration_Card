# =============================================================================
# 16b_akm_standard_errors.R  --  Adao, Kolesar & Morales (2019) exposure-robust
#                                standard errors
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# WHAT THIS COMPUTES, AND HOW IT DIFFERS FROM 16
# 16_adao_placebo.R runs a placebo IV and compares MSA- against state-clustered
# standard errors. Neither is the AKM correction. The placebo asks whether the
# instrument predicts outcomes it should not, which is a validity check. The
# clustering comparison coarsens the GEOGRAPHIC cluster, and the correlation
# AKM describe runs across units with similar SHARE VECTORS whether or not they
# are geographic neighbours. This script computes the variance estimator
# itself, using the authors' own R implementation (ShiftShareSE, Kolesar).
#
# THE UNIT OF ANALYSIS
# AKM's estimator is defined on a regression whose unit carries a vector of
# exposure shares. The headline regression is individual-level with 3.5 million
# records; the shift-share structure lives at MSA x birth cohort, where
# treatment is assigned. The windowed instrument collapses to that level
# exactly, because the 1990 shares do not vary with the cohort:
#
#   z_{m,c} = (1/|W_c|) sum_{t in W_c} sum_k s_{m,k} * raw_count_{k,t}
#           = sum_k s_{m,k} * Delta_{k,c}
#
# so the design is a shift-share at MSA x cohort with shocks indexed by origin
# x cohort.
#
# HOW THE CONTROLS ARE HANDLED, AND WHY IT IS NOT BY AGGREGATION
# The individual controls cannot be carried to the cell level as cell means.
# Cell means of parental education, poverty and child-age composition are MSA x
# cohort composition variables that co-move with immigrant exposure, so they
# absorb the between-cell variation the design is identified on: 16c measured
# the damage, with the cell estimate collapsing from +5.73e-05 to -1.39e-06
# once aggregated controls were added.
#
# Frisch-Waugh gives the right route. Residualising y and x on the controls and
# fixed effects AT THE INDIVIDUAL LEVEL and then collapsing reproduces the
# headline coefficient exactly, because z is constant within MSA x cohort:
#
#   z'y~ = sum_i z_i y~_i = sum_{(m,c)} z_{m,c} n_{m,c} mean(y~)_{m,c}
#
# and the same for x. Leaving z itself unresidualised is what preserves the
# W * Delta structure AKM needs, and costs nothing: the annihilator is
# idempotent and symmetric, so z'y~ = z~'y~. Step 4 asserts the reproduction
# and stops if it fails.
#
# WHAT TO READ OFF THE OUTPUT
# The transportable quantity is the RATIO of the AKM standard error to the
# conventionally clustered one on the SAME regression. A ratio near one means
# exposure correlation is not inflating precision here; a ratio well above one
# means the chapter's reported inference is optimistic by about that factor.
#
# Shocks are clustered on origin, so the 33 cohort-specific shocks drawn from
# one origin are treated as a single cluster; national inflows are serially
# correlated within origin. The same fit is run with shocks unclustered, to
# show how much of the answer is that choice.
#
# Writes outputs/tables/akm_se_results.csv.
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(fixest)

if (!requireNamespace("ShiftShareSE", quietly = TRUE))
  stop("Install ShiftShareSE first:  install.packages(\"ShiftShareSE\")")
library(ShiftShareSE)

# --- Inputs ------------------------------------------------------------------
cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
all_natives_175 <- as.data.table(.tmp_env$all_natives_175)
selected_msas   <- .tmp_env$selected_msas
rm(.tmp_env); gc()

cat("Loading the shift-share panel from the 06 checkpoint...\n")
.ss_path <- file.path(local_checkpoint_dir, "after_sec06.RData")
if (!file.exists(.ss_path))
  stop("Checkpoint not found at ", .ss_path,
       "\nThis script needs the share x shock panel that 06_shift_share_iv.R ",
       "writes. Re-run 06, or set local_checkpoint_dir to where it saved.")
.ss_env <- new.env()
load(.ss_path, envir = .ss_env)
if (!exists("shift_share_data", envir = .ss_env))
  stop("after_sec06.RData does not contain shift_share_data.")
shift_share_data <- as.data.table(.ss_env$shift_share_data)
rm(.ss_env); gc()

# --- Step 1: the estimation sample, matching 09_main_results.R ---------------
d <- all_natives_175[exposure_missing == 0 & MSA %in% selected_msas]
d[, birth_year := YEAR - AGE]
d[, z := avg_exposure_nodiv]
d[, x := avg_actual_exposure_nodiv]
d[, y_HSplus   := as.integer(Educ_Group != "non_HS")]
d[, y_CollPlus := as.integer(Educ_Group %in% c("College", "Grad"))]

for (.yv in c("y_HSplus", "y_CollPlus")) {
  if (uniqueN(d[[.yv]]) < 2L)
    stop(.yv, " is constant. Educ_Group levels present: ",
         paste(sort(unique(as.character(d$Educ_Group))), collapse = ", "))
}

cat(sprintf("\n  individual records      %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  MSAs                    %d\n", uniqueN(d$MSA)))
cat(sprintf("  birth cohorts           %d (%d-%d)\n", uniqueN(d$birth_year),
            min(d$birth_year), max(d$birth_year)))

# The collapse is exact only if both sides of the first stage are constant
# within the cell.
.chk <- d[, .(nz = uniqueN(round(z, 6)), nx = uniqueN(round(x, 6))),
          by = .(MSA, birth_year)]
if (.chk[, max(nz)] > 1L || .chk[, max(nx)] > 1L)
  stop("z or x varies within MSA x cohort; the collapse would not be lossless.")
rm(.chk)

# --- Step 2: individual-level benchmark --------------------------------------
CTRL <- "Black + SEX + Marital + ChildAgeGroup + POVERTY + Mom_Educ + Dad_Educ |
         YEAR + MSA"

cat("\nIndividual-level benchmark (headline specification):\n")
bench <- list()
for (yv in c("y_HSplus", "y_CollPlus")) {
  iv <- feols(as.formula(sprintf("%s ~ %s | x ~ z", yv, CTRL)),
              data = d, vcov = ~MSA)
  ct <- iv$coeftable["fit_x", ]
  bench[[yv]] <- c(est = ct[["Estimate"]], se = ct[["Std. Error"]])
  cat(sprintf("  %-11s 2SLS %+.6e   SE(MSA) %.4e   t %6.2f\n",
              yv, ct[["Estimate"]], ct[["Std. Error"]], ct[["t value"]]))
}

# --- Step 3: residualise at the individual level, then collapse --------------
cat("\nResidualising on controls and fixed effects...\n")
d[, rx := resid(feols(as.formula(sprintf("x ~ %s", CTRL)), data = d))]
for (yv in c("y_HSplus", "y_CollPlus"))
  set(d, j = paste0("r_", yv),
      value = resid(feols(as.formula(sprintf("%s ~ %s", yv, CTRL)), data = d)))

cell <- d[, .(n = .N,
              r_y_HSplus   = mean(r_y_HSplus),
              r_y_CollPlus = mean(r_y_CollPlus),
              rx = mean(rx), z = mean(z)),
          by = .(MSA, birth_year)]
setorder(cell, MSA, birth_year)
cat(sprintf("  cells (MSA x cohort)    %s\n", format(nrow(cell), big.mark = ",")))
cat(sprintf("  median cell size        %.0f\n", cell[, median(n)]))

# --- Step 4: the collapse must reproduce the headline ------------------------
cat("\nCell-level reproduction of the headline estimate:\n")
cell_fit <- list()
for (yv in c("y_HSplus", "y_CollPlus")) {
  m <- feols(as.formula(sprintf("r_%s ~ 1 | rx ~ z", yv)),
             data = cell, weights = ~n, vcov = ~MSA)
  ct <- m$coeftable["fit_rx", ]
  cell_fit[[yv]] <- ct
  r <- ct[["Estimate"]] / bench[[yv]]["est"]
  cat(sprintf("  %-11s cell %+.6e   individual %+.6e   ratio %.6f\n",
              yv, ct[["Estimate"]], bench[[yv]]["est"], r))
  if (abs(r - 1) > 1e-4)
    stop("The cell-level collapse does not reproduce the individual estimate ",
         "for ", yv, " (ratio ", signif(r, 6), "). Frisch-Waugh makes this an ",
         "identity, so a mismatch means the residualisation or the weighting ",
         "is wrong, not that aggregation costs something.")
}

# --- Step 5: the share matrix and the shock vector ---------------------------
cohorts <- sort(unique(cell$birth_year))
msas    <- sort(unique(cell$MSA))

shares <- unique(shift_share_data[MSA %in% msas, .(MSA, BPLD, share)])
origins <- sort(unique(shares$BPLD))
if (shares[, .N, by = .(MSA, BPLD)][, max(N)] > 1L)
  stop("shift_share_data holds more than one share for some MSA-origin pair.")
cat(sprintf("\n  origins                 %d\n", length(origins)))

shocks_raw <- unique(shift_share_data[, .(BPLD, YEAR, raw_count)])
shock_list <- rbindlist(lapply(cohorts, function(cc) {
  yrs <- seq(cc + 15L, min(cc + 25L, 2022L))
  shocks_raw[YEAR %in% yrs, .(birth_year = cc, shock = mean(raw_count)),
             by = BPLD]
}))
shock_list <- shock_list[CJ(BPLD = origins, birth_year = cohorts),
                         on = .(BPLD, birth_year)]
shock_list[is.na(shock), shock := 0]

share_mat <- dcast(shares, MSA ~ BPLD, value.var = "share", fill = 0)
.pos <- match(msas, share_mat$MSA)
if (anyNA(.pos))
  stop(sum(is.na(.pos)), " estimation MSAs have no row in the 1990 share matrix.")
share_mat <- share_mat[.pos]
S <- as.matrix(share_mat[, -1L])
rownames(S) <- as.character(msas)
S <- S[, as.character(origins), drop = FALSE]

# There are more origins than MSAs, so the share matrix cannot have full column
# rank: at most 175 of the 230 origins carry exposure variation that is not a
# combination of the others. This is reported because it bears on how much
# independent shock-level variation the variance can draw on, and because
# ShiftShareSE will warn about it.
#
# EVERY ORIGIN IS KEPT ALL THE SAME. Rank deficiency means a dropped origin's
# share vector lies in the span of the retained ones, not that its shock is
# redundant: z = sum_k s_k * Delta_k needs each origin's own term, so removing
# a column deletes s_j * Delta_j from the instrument. Reconstructing z from a
# reduced set would require folding the dropped shocks into the retained
# columns, which would make the shocks linear combinations rather than national
# inflows -- and the shocks being the actual quasi-random objects is the whole
# premise of the method. Dropping collinear columns is legitimate inside the
# variance calculation, which is where the package does it.
.rank_S <- qr(S)$rank
cat(sprintf("  share matrix rank       %d of %d origins (%d MSAs)\n",
            .rank_S, ncol(S), nrow(S)))
cat(sprintf("  origins with no independent exposure: %d (retained; the\n",
            ncol(S) - .rank_S))
cat("    collinearity warning from ShiftShareSE is expected)\n")
setorder(shock_list, birth_year, BPLD)

n_cell <- nrow(cell); n_shock <- length(origins) * length(cohorts)
cat(sprintf("  share matrix            %s x %s (%.2f GB dense)\n",
            format(n_cell, big.mark = ","), format(n_shock, big.mark = ","),
            n_cell * n_shock * 8 / 1024^3))

W <- matrix(0, nrow = n_cell, ncol = n_shock)
for (i in seq_along(cohorts)) {
  rows <- which(cell$birth_year == cohorts[i])
  cols <- ((i - 1L) * length(origins) + 1L):(i * length(origins))
  W[rows, cols] <- S[as.character(cell$MSA[rows]), , drop = FALSE]
}
shock_vec <- shock_list$shock

z_recon <- as.vector(W %*% shock_vec)
.relv <- abs(z_recon - cell$z) / pmax(abs(cell$z), 1e-8)
cat(sprintf("\n  max relative error, W %%*%% shocks vs z:  %.3e\n", max(.relv)))
cat(sprintf("  correlation(recon, z):                  %.8f\n",
            cor(z_recon, cell$z)))
if (max(.relv) > 1e-6)
  stop("The reconstructed instrument does not match avg_exposure_nodiv. The ",
       "AKM variance would be computed on the wrong exposure structure.")

sector_cluster <- rep(seq_along(origins), times = length(cohorts))

# --- Step 6: AKM inference ---------------------------------------------------
KEY <- c(ehw = "EHW", cluster = "Reg. cluster", akm = "AKM", akm0 = "AKM0")
pick <- function(f, slot, m) {
  v <- f[[slot]]; k <- KEY[[m]]
  if (is.null(v) || !(k %in% names(v))) return(NA_real_)
  unname(v[[k]])
}

results <- list()
for (yv in c("y_HSplus", "y_CollPlus")) {
  cat(sprintf("\n%s\nAKM inference: %s\n%s\n", strrep("=", 72), yv,
              strrep("=", 72)))
  ct_cl <- cell_fit[[yv]]

  .args <- list(as.formula(sprintf("r_%s ~ 1 | rx", yv)),
                X = cell$z, data = cell, W = W,
                region_cvar = cell$MSA, sector_cvar = sector_cluster,
                method = c("ehw", "region_cluster", "akm", "akm0"))
  if ("weights" %in% names(formals(ShiftShareSE::ivreg_ss))) {
    .args$weights <- cell$n
  } else if (yv == "y_HSplus") {
    cat("  NOTE: this ivreg_ss takes no weights argument; cells enter",
        "unweighted.\n")
  }
  fit <- do.call(ShiftShareSE::ivreg_ss, .args)

  # The same fit treating each origin-cohort shock as independent. This is not
  # an alternative estimate to weigh against the one above: the instrument
  # averages national inflows over an eleven-year window, so adjacent cohorts
  # share ten of their eleven years and independence across cohorts within an
  # origin is ruled out by construction. It is computed to show what that
  # assumption does, and is not reported in the chapter.
  .args_nc <- .args; .args_nc$sector_cvar <- NULL
  fit_nc <- do.call(ShiftShareSE::ivreg_ss, .args_nc)

  base_se <- pick(fit, "se", "cluster")
  if (is.na(base_se)) base_se <- ct_cl[["Std. Error"]]

  cat(sprintf("  estimate %+.6e  (headline %+.6e)\n",
              ct_cl[["Estimate"]], bench[[yv]]["est"]))
  cat(sprintf("\n  %-24s %-12s %-10s %s\n", "method", "SE", "ratio", "95% CI"))
  cat(sprintf("  %-24s %.4e   %-8s   [%+.3e, %+.3e]\n", "MSA-clustered (feols)",
              ct_cl[["Std. Error"]], "--",
              ct_cl[["Estimate"]] - 1.96 * ct_cl[["Std. Error"]],
              ct_cl[["Estimate"]] + 1.96 * ct_cl[["Std. Error"]]))
  for (m in c("cluster", "ehw", "akm", "akm0")) {
    s <- pick(fit, "se", m)
    cat(sprintf("  %-24s %.4e   %-8.3f   [%+.3e, %+.3e]\n", KEY[[m]], s,
                s / base_se, pick(fit, "ci.l", m), pick(fit, "ci.r", m)))
  }
  cat("\n  Shocks unclustered (each origin-cohort its own shock):\n")
  for (m in c("akm", "akm0")) {
    s <- pick(fit_nc, "se", m)
    cat(sprintf("  %-24s %.4e   %-8.3f   [%+.3e, %+.3e]\n",
                paste0(KEY[[m]], ", no clust"), s, s / base_se,
                pick(fit_nc, "ci.l", m), pick(fit_nc, "ci.r", m)))
  }

  results[[yv]] <- data.table(
    outcome = yv,
    est_individual = bench[[yv]]["est"], se_individual = bench[[yv]]["se"],
    est_cell = ct_cl[["Estimate"]], se_cluster = ct_cl[["Std. Error"]],
    se_cluster_akm = pick(fit, "se", "cluster"),
    se_ehw = pick(fit, "se", "ehw"),
    se_akm = pick(fit, "se", "akm"), se_akm0 = pick(fit, "se", "akm0"),
    se_akm_noclust = pick(fit_nc, "se", "akm"),
    se_akm0_noclust = pick(fit_nc, "se", "akm0"),
    akm_ci_l = pick(fit, "ci.l", "akm"), akm_ci_r = pick(fit, "ci.r", "akm"),
    # AKM0 inverts a test, so its interval is not the estimate plus and minus a
    # multiple of its standard error and has to be carried separately.
    akm0_ci_l = pick(fit, "ci.l", "akm0"), akm0_ci_r = pick(fit, "ci.r", "akm0"),
    cluster_ci_l = pick(fit, "ci.l", "cluster"),
    cluster_ci_r = pick(fit, "ci.r", "cluster"))
}

out <- rbindlist(results)
out[, `:=`(ratio_akm = se_akm / se_cluster, ratio_akm0 = se_akm0 / se_cluster,
           n_cells = n_cell, n_origins_used = length(origins),
           n_msas = nrow(S), rank_shares = .rank_S)]
fwrite(out, "outputs/tables/akm_se_results.csv")
cat("\nWrote outputs/tables/akm_se_results.csv\n")
print(out)
