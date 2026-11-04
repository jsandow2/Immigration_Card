# =============================================================================
# 18_rate_spec_loo.R  --  Does the rate specification survive dropping the
# Renumbered 2026-10-21 from 20j_rate_spec_loo.R (archive-and-renumber pass).
#                          metros that carry it?
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# THE QUESTION
# 20i found that expressing inflows as rates rather than counts more than
# doubles the effective cluster count (G* 2.8 -> 6.3), breaks New York's 73.3%
# grip on the variance, and moves the t-statistic from 1.78 to 3.08. That looks
# like progress until you read which metros replaced New York:
#
#   35620 New York        16.4%
#   31080 Los Angeles     14.9%
#   32580 McAllen         13.2%   <- Rio Grande Valley border metro
#   41860 San Francisco   11.9%
#   33100 Miami           10.8%
#
# and in the PERWT-denominator variant, 15180 Brownsville-Harlingen at 11.8%.
#
# WHY THE BORDER METROS ARE A PROBLEM AND NOT JUST A CURIOSITY
#
#   1. The denominator does this mechanically. Dividing inflow by base
#      population inflates leverage wherever inflow is large relative to size.
#      Nowhere in the country fits that better than the lower Rio Grande
#      Valley, so the rate specification recruits those metros by construction
#      rather than because they are informative.
#
#   2. It reintroduces the Mexico concentration from Section 19 through a
#      different door. McAllen's and Brownsville's 1990 origin shares are
#      overwhelmingly Mexican, so their predicted inflows are close to the
#      national Mexican arrival series scaled by a large share. If those two
#      metros jointly carry ~25% of the identifying variation, the instrument
#      is substantially ONE shock -- exactly the failure mode the Rotemberg
#      decomposition exists to detect.
#
#   3. Border metros are not the population the claim is about. Native
#      educational attainment there is entangled with cross-border commuting,
#      colonias, and a Mexican-origin native population whose schooling
#      trajectory differs systematically from the national one.
#
# WHAT THIS SCRIPT DECIDES
# Leverage is not the same as dependence. A cluster can carry a large share of
# the variance and still not move the point estimate when removed. This script
# separates the two:
#
#   - Leave-one-out over the top contributors to the rate specification
#   - Drop both border metros together (the substantive test)
#   - Drop New York and Los Angeles together (the contrast: if the estimate is
#     equally fragile to removing the large gateways, fragility is a property
#     of the design, not of the border metros specifically)
#   - The count specification without New York, as an anchor, so the two
#     specifications are judged on the same page
#
# HOW TO READ IT -- POINT ESTIMATES, NOT P-VALUES
# At G* around 6 no standard error here is trustworthy, and that includes the
# ones this script prints. The t-statistics are reported for comparability with
# 20i and should not be interpreted as tests. The decision rests on whether the
# POINT ESTIMATE is stable across removals:
#
#   Stable estimate  -> the concentration is leverage. The rate specification
#                       is defensible, with the concentration disclosed.
#   Estimate moves   -> the apparent gain in 20i was those metros. The rate
#                       specification should not carry the chapter.
#
# A t-statistic that collapses while the estimate holds is expected and is NOT
# evidence against the spec: dropping a high-leverage cluster removes variation
# and widens the SE by construction.
#
# Nothing written to the pipeline.
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

cat("Loading analysis_data.RData...\n")
load("analysis_data.RData")
setDT(reduced_data_iv)

# Names for the metros that appear in the diagnostics. Verified against
# data/crosswalks/cbsa_2013_official.xls on 2026-09-24.
MSA_NAMES <- c(
  "35620" = "New York",        "31080" = "Los Angeles",
  "32580" = "McAllen",         "41860" = "San Francisco",
  "33100" = "Miami",           "15180" = "Brownsville",
  "21340" = "El Paso",         "29700" = "Laredo",
  "16980" = "Chicago",         "47900" = "Washington DC",
  "26420" = "Houston",         "19100" = "Dallas",
  "42660" = "Seattle",         "38060" = "Phoenix",
  "19740" = "Denver",          "40140" = "Riverside",
  "41940" = "San Jose",        "12540" = "Bakersfield",
  "41700" = "San Antonio",     "45300" = "Tampa"
)
nm <- function(code) {
  code <- as.character(code)
  unname(ifelse(code %in% names(MSA_NAMES), MSA_NAMES[code], code))
}

# All four Texas-Mexico border metros in the 2013 delineation, not just the two
# that surfaced in 20i's top five. They share near-identical 1990 origin
# composition, so their predicted inflows are close to collinear by
# construction -- dropping only the visible ones would leave the concentration
# partly intact and understate the problem.
BORDER_METROS  <- c(32580, 15180, 21340, 29700)  # McAllen, Brownsville, El Paso, Laredo
GATEWAY_METROS <- c(35620, 31080)                # New York, Los Angeles

# --- Base-year size: raw respondent count, natives 4x ------------------------
# Raw counts, matching the units of the raw-count numerator built in
# 06_shift_share_iv.R. See 20i for the full argument.
pop <- reduced_data_iv[MSA != 0, .(
  n_raw = sum(fifelse(immigr == 0L, 4, 1))
), by = .(MSA, YEAR)]
setorder(pop, MSA, YEAR)

pop_1990 <- pop[YEAR == 1990, .(MSA, b = n_raw)]
earliest <- pop[, .SD[1], by = MSA][, .(MSA, f = n_raw)]
pop_base <- merge(earliest, pop_1990, by = "MSA", all.x = TRUE)
pop_base[, base_raw := fifelse(is.na(b), f, b)]
pop_base <- pop_base[, .(MSA, base_raw)]

# --- Estimation sample -------------------------------------------------------
AGE_FLOOR <- 30L
d <- reduced_data_iv[immigr == 0 & MSA != 0 & AGE >= AGE_FLOOR]
d <- merge(d, msa_year_predicted_inflows, by = c("MSA", "YEAR"), all.x = TRUE)
d[is.na(total_predicted_inflow), total_predicted_inflow := 0]
d <- merge(d, pop_base, by = "MSA", all.x = TRUE)
d[, y_HSplus := as.integer(Educ_Group != "non_HS")]
d <- d[MSA %in% selected_msas]

d[, `:=`(
  actual_rate   = 1000 * total_actual_inflow / base_raw,
  exposure_rate = 1000 * avg_exposure        / base_raw
)]

cat(sprintf("\nBaseline sample: %d obs | %d MSAs | %d identifying\n",
            nrow(d), uniqueN(d$MSA), d[exposure_missing == 0, .N]))

# --- Estimation and diagnostics ---------------------------------------------
CTRL <- "Black + SEX + Marital + ChildAgeGroup + POVERTY +
         Mom_Educ + Dad_Educ + exposure_missing | YEAR + MSA"

fit_one <- function(dd, endog, instr, label, verbose = TRUE) {
  dd <- droplevels(dd)
  f  <- as.formula(sprintf("y_HSplus ~ %s | %s ~ %s", CTRL, endog, instr))
  m  <- feols(f, data = dd, vcov = ~MSA)
  ct <- m$coeftable[paste0("fit_", endog), ]

  z_t <- residuals(feols(as.formula(paste(instr, "~", CTRL)), data = dd))
  msa <- as.character(dd$MSA)

  gamma <- tapply(z_t^2, msa, sum)
  G     <- length(gamma)
  cv2   <- mean((gamma - mean(gamma))^2) / mean(gamma)^2
  Gstar <- G / (1 + cv2)

  score <- tapply(z_t * residuals(m), msa, sum)
  sh    <- sort(score^2 / sum(score^2), decreasing = TRUE)

  if (verbose) {
    cat(sprintf("\n=== %s ===\n", label))
    cat(sprintf("  n %d | MSAs %d | identifying %d\n",
                nrow(dd), G, dd[exposure_missing == 0, .N]))
    cat(sprintf("  Estimate %.5e  SE %.5e  t %.3f  p %.4f\n",
                ct[["Estimate"]], ct[["Std. Error"]],
                ct[["t value"]], ct[["Pr(>|t|)"]]))
    cat(sprintf("  CV^2 %.2f | G* %.1f | top cluster %.1f%% (%s) | top 5 %.1f%%\n",
                cv2, Gstar, 100 * sh[1], nm(names(sh)[1]), 100 * sum(sh[1:5])))
    cat(sprintf("  Largest: %s\n",
                paste(sprintf("%s %.1f%%", nm(names(sh)[1:5]), 100 * sh[1:5]),
                      collapse = " | ")))
  }

  list(row = data.frame(
         spec = label, n = nrow(dd), msas = G,
         estimate = ct[["Estimate"]], se = ct[["Std. Error"]],
         t = ct[["t value"]], p = ct[["Pr(>|t|)"]],
         cv2 = cv2, Gstar = Gstar,
         top1 = unname(sh[1]), top1_msa = nm(names(sh)[1]),
         stringsAsFactors = FALSE),
       shares = sh)
}

# --- Baseline ----------------------------------------------------------------
cat("\n", strrep("-", 72), "\n", sep = "")
cat("RATE SPECIFICATION -- BASELINE\n")
base <- fit_one(d, "actual_rate", "exposure_rate", "Rate: all 175 MSAs")
base_est <- base$row$estimate
base_ci  <- base_est + c(-1.96, 1.96) * base$row$se

# How much is actually at stake in the border group, and which of the four
# survive the top-175 filter at all.
present <- BORDER_METROS[BORDER_METROS %in% unique(d$MSA)]
absent  <- setdiff(BORDER_METROS, present)
cat("\nBorder-metro group\n")
cat(sprintf("  in sample     : %s\n",
            if (length(present)) paste(nm(present), collapse = ", ") else "none"))
cat(sprintf("  not in top-175: %s\n",
            if (length(absent)) paste(nm(absent), collapse = ", ") else "none"))
if (length(present)) {
  bs <- base$shares[names(base$shares) %in% as.character(present)]
  cat(sprintf("  shares        : %s\n",
              paste(sprintf("%s %.1f%%", nm(names(bs)), 100 * bs), collapse = " | ")))
  cat(sprintf("  combined      : %.1f%% of baseline variance\n", 100 * sum(bs)))
  cat(sprintf("  observations  : %d (%.2f%% of sample)\n",
              d[MSA %in% present, .N], 100 * d[MSA %in% present, .N] / nrow(d)))
}

# --- Leave-one-out over the top contributors --------------------------------
top_codes <- names(base$shares)[1:6]

cat("\n", strrep("-", 72), "\n", sep = "")
cat("RATE SPECIFICATION -- LEAVE-ONE-OUT OVER TOP CONTRIBUTORS\n")
loo <- lapply(top_codes, function(cd) {
  fit_one(d[as.character(MSA) != cd], "actual_rate", "exposure_rate",
          sprintf("Rate: drop %s", nm(cd)))
})

# --- Joint drops -------------------------------------------------------------
cat("\n", strrep("-", 72), "\n", sep = "")
cat("RATE SPECIFICATION -- JOINT DROPS\n")
drop_border  <- fit_one(d[!MSA %in% BORDER_METROS],
                        "actual_rate", "exposure_rate",
                        "Rate: drop McAllen + Brownsville")
drop_gateway <- fit_one(d[!MSA %in% GATEWAY_METROS],
                        "actual_rate", "exposure_rate",
                        "Rate: drop New York + Los Angeles")

# --- Count specification anchor ---------------------------------------------
cat("\n", strrep("-", 72), "\n", sep = "")
cat("COUNT SPECIFICATION -- ANCHOR (different units; compare stability only)\n")
cnt_base <- fit_one(d, "total_actual_inflow", "avg_exposure",
                    "Count: all 175 MSAs")
cnt_dropny <- fit_one(d[MSA != 35620], "total_actual_inflow", "avg_exposure",
                      "Count: drop New York")

# --- Summary -----------------------------------------------------------------
res <- do.call(rbind, c(list(base$row),
                        lapply(loo, `[[`, "row"),
                        list(drop_border$row, drop_gateway$row,
                             cnt_base$row, cnt_dropny$row)))

cat("\n", strrep("-", 72), "\n", sep = "")
cat("SUMMARY\n\n")
print(res[, c("spec", "msas", "estimate", "se", "t", "Gstar", "top1", "top1_msa")],
      row.names = FALSE, digits = 3)

# Stability of the rate estimate, expressed relative to the baseline CI.
rate_rows <- res[grepl("^Rate", res$spec), ]
rate_rows$pct_change  <- 100 * (rate_rows$estimate - base_est) / base_est
rate_rows$inside_ci   <- rate_rows$estimate >= base_ci[1] &
                         rate_rows$estimate <= base_ci[2]

cat("\n\nRate-specification point-estimate stability\n")
cat(sprintf("Baseline %.5e, 95%% CI [%.5e, %.5e]\n",
            base_est, base_ci[1], base_ci[2]))
print(rate_rows[, c("spec", "estimate", "pct_change", "inside_ci")],
      row.names = FALSE, digits = 3)

cat("\n", strrep("-", 72), "\n", sep = "")
cat("READING THIS\n\n")
cat("Judge on the point estimate, not on t or p. Standard errors at G* ~ 6 are\n")
cat("not trustworthy and are printed only for comparability with 20i.\n\n")
cat("If every rate estimate stays inside the baseline CI and the percentage\n")
cat("changes are modest, the top metros carry LEVERAGE but not DEPENDENCE: the\n")
cat("rate specification is defensible provided the concentration is disclosed.\n\n")
cat("If dropping McAllen and Brownsville moves the estimate materially while\n")
cat("dropping New York and Los Angeles does not, the rate specification's\n")
cat("apparent improvement in 20i was purchased from two border metros whose\n")
cat("identifying variation is close to a single national Mexican shock, and it\n")
cat("should not carry the chapter.\n\n")
cat("If BOTH joint drops move the estimate materially, fragility is a property\n")
cat("of the design rather than of any particular metro, and neither\n")
cat("specification supports a precise claim -- in which case the honest report\n")
cat("is the estimate with an explicit statement that cluster-robust inference\n")
cat("is unavailable at this effective cluster count.\n")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs/tables/rate_spec_loo.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/rate_spec_loo.csv\n")
cat("18_rate_spec_loo.R complete.\n")
