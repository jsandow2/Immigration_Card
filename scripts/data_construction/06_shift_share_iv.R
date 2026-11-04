# =============================================================================
# 06_shift_share_iv.R  --  Section 6: Shift-share IV construction
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# =============================================================================

source("scripts/data_construction/00_setup.R")
load(file.path(local_checkpoint_dir, "after_sec05.RData"))
cat("Loaded checkpoint: after_sec05.RData\n")

# Section 6: Shift-share IV construction (data.table rewrite)
#
# Shift-share predicted inflow for MSA m in year t:
#   predicted_inflow_{m,t} = sum_k share_{m,k,1990} * national_inflow_{k,t}
#
# where share_{m,k,1990} is the 1990 share of immigrants from origin k in
# MSA m, and national_inflow_{k,t} is the national inflow of immigrants
# from origin k in year t.
#
# Construction steps:
#   1. Build 1990 shares by MSA and BPLD
#   2. Build national inflows by year and BPLD
#   3. Create full MSA x BPLD x YEAR grid
#   4. Multiply shares x shifts to get predicted inflows
#   5. Aggregate to MSA x YEAR level
#
# METHODOLOGY NOTE (2026-09-11): national_inflow and actual_inflow (built
# below) use raw respondent counts pooled across every survey wave, not
# PERWT-weighted sums restricted to a single closest wave. Three design
# choices worth being explicit about:
#
#   - Raw counts, not PERWT. PERWT is calibrated to reproduce each survey
#     year's own population totals -- it is not designed to measure
#     historical arrival flows. Summing PERWT across multiple waves for
#     the same (BPLD, YRIMMIG) cell independently rescales each wave to
#     the full national population and then adds those already-inflated
#     estimates together, which can produce implausible totals (e.g. a
#     single arrival cohort appearing larger than that origin's entire
#     US-resident stock). Raw respondent counts, pooled across waves,
#     avoid this: each wave samples different real individuals, so
#     pooling adds genuine additional observations of the same cohort
#     rather than re-inflating the same estimate multiple times over.
#     reduced_data has already been thinned in 04_sample_construction.R so
#     1990/2000 match the other years' sampling density, which is what
#     makes raw counts comparable across waves in the first place.
#
#   - Pooled across all waves, not restricted to the closest one. This
#     lets the inflow series use the full range of arrival years present
#     in the data (back to the early 1900s) rather than dropping cohorts
#     that fall far from any single survey wave. The accepted tradeoff:
#     accuracy is best for recent arrival cohorts and degrades gradually
#     further back in time, as mortality and return migration shrink the
#     pool of survivors available to be sampled, and as self-reported
#     YRIMMIG becomes more prone to recall/rounding error. This is a
#     gradual effect, not a sharp cutoff, and becomes material mainly for
#     cohorts roughly 20-30+ years before a given survey wave -- rounding
#     in particular is uncommon for recent arrivals and mostly shows up
#     that far back. Native-born individuals whose exposure window falls
#     in that degraded range are intended to be caught by
#     exposure_missing in 07_exposure_iv.R.
#
#   - Divided by the number of eligible waves. Pooling raw counts across
#     waves does not double-count people (each wave samples different real
#     individuals), but without a divisor it mechanically inflates recent
#     arrival cohorts relative to older ones for a reason that has nothing
#     to do with actual immigration: a cohort that arrived in 2006 can be
#     observed by every ACS wave from 2006 through the present (many
#     chances to be sampled), while a cohort that arrived in 2023 can only
#     be observed by the 2023 wave (one chance). Left uncorrected this
#     produces a spurious ~15x mechanical decline in measured inflow across
#     the most recent part of the series. n_waves = the count of survey
#     waves with YEAR >= YRIMMIG (i.e., waves that could have observed that
#     arrival cohort at all); national_inflow/actual_inflow = raw pooled
#     count / n_waves. Arrival years with fewer than 2 eligible waves --
#     under the current data, only YRIMMIG == 2023, observed by just the
#     2023 wave itself -- are dropped entirely rather than divided, since a
#     single-wave count is also only ~50% complete (people who will move
#     into a covered MSA, age into the sample, etc. after that single
#     survey can never appear) and averaging that one thin, partial wave
#     with nothing else doesn't fix that.
# -----------------------------------------------------------------------------

# Ensure reduced_data is a data.table (it already is at this point in the
# pipeline, but setDT is idempotent so this is safe)
setDT(reduced_data)

# Step 1: 1990 immigrant shares by MSA and origin (BPLD)
#
# WHY THE BASE YEAR IS 1990. Two conditions have to hold at once. The shares
# must be measured after the Immigration and Nationality Act of 1965 reshaped
# origin composition, so that they describe the settlement of the population
# whose later arrivals are being instrumented, and long enough after it for
# those patterns to have consolidated geographically. They must also stay
# close to the period the exposure window averages over, since an eleven-year
# window needs shares near it to retain a first stage. 1990 satisfies both.
# The argument as it reaches the reader is in the chapter's identification
# section.
#
# Share definition follows Card (2001): share_{m,k,1990} is the fraction of
# all U.S. immigrants from origin k in 1990 who were located in MSA m.
# Denominator is the national total from origin k, NOT the total immigrants
# within MSA m. This is the national share formulation and is less vulnerable
# to the endogeneity concern that local immigrant composition reflects local
# demand conditions. Computed from raw respondent counts within the single
# 1990 wave (not PERWT), for consistency with national_inflow/actual_inflow
# below.
#
# The Top-175 MSA ranking in 05_stable_msa_selection.R is PERWT-weighted, a
# separate and much lower-stakes choice: 173 of the 175 metros overlap with
# the raw-count ranking either way.

immigrants_1990 <- reduced_data[immigr == 1 & YEAR == 1990,
                                .(count = .N),
                                by = .(MSA, BPLD)]

immigrants_1990[, national_total := sum(count), by = BPLD]
immigrants_1990[, share := count / national_total]

immigrant_shares <- immigrants_1990[, .(MSA, BPLD, share)]

# Step 2: National inflows by year and origin (BPLD)
#
# Raw respondent count, pooled across every survey wave and divided by the
# number of eligible waves, by (arrival year, origin). See the METHODOLOGY
# NOTE above for why this is a raw, divided count rather than a
# PERWT-weighted, single-wave figure.
#
# arrivals is built once here and reused below for actual_inflows_raw, so
# the wave/eligibility logic (and the immigr==1 & YRIMMIG>0 filter) is only
# written once.
waves <- sort(unique(reduced_data$YEAR))
arrivals <- reduced_data[immigr == 1 & YRIMMIG > 0,
                         .(wave = as.integer(YEAR), YRIMMIG = as.integer(YRIMMIG),
                           BPLD = as.integer(BPLD), MSA = MSA)]

# n_waves = number of survey waves that could possibly have observed a given
# arrival cohort (waves surveyed in or after the arrival year). Arrival
# years with fewer than 2 eligible waves -- currently just YRIMMIG == 2023,
# observed only by the 2023 wave itself -- are dropped rather than divided;
# see METHODOLOGY NOTE above.
eligible_waves <- data.table(YRIMMIG = sort(unique(arrivals$YRIMMIG)))
eligible_waves[, n_waves := vapply(YRIMMIG, function(y) sum(waves >= y), integer(1))]
eligible_waves <- eligible_waves[n_waves >= 2L]

national_inflows <- arrivals[, .(raw_count = .N), by = .(YEAR = YRIMMIG, BPLD)]
national_inflows <- merge(national_inflows, eligible_waves,
                          by.x = "YEAR", by.y = "YRIMMIG")
national_inflows[, national_inflow := raw_count / n_waves]

# Step 3: Full MSA x BPLD grid from 1990 MSAs and observed origins
all_msas_1990 <- unique(reduced_data[YEAR == 1990, .(MSA)])
all_bplds     <- unique(immigrant_shares[, .(BPLD)])

# Cross join: every MSA x every BPLD
msa_bpld_grid <- CJ(MSA = all_msas_1990$MSA, BPLD = all_bplds$BPLD)

# Step 4: Merge shares; unobserved MSA-origin pairs get share = 0
valid_pairs <- merge(msa_bpld_grid, immigrant_shares,
                     by = c("MSA", "BPLD"), all.x = TRUE)
valid_pairs[is.na(share), share := 0]

# Step 5: Expand to all years and compute predicted inflows
valid_years   <- unique(national_inflows$YEAR)
years_dt      <- data.table(YEAR = valid_years)

# Cross join valid_pairs x valid_years
# Temporary key column to enable cartesian merge
valid_pairs[, .cross := 1L]
years_dt[,   .cross := 1L]
shift_share_data <- merge(valid_pairs, years_dt, by = ".cross", allow.cartesian = TRUE)
shift_share_data[, .cross := NULL]
valid_pairs[,    .cross := NULL]

# Merge national inflows
shift_share_data <- merge(shift_share_data, national_inflows,
                          by = c("YEAR", "BPLD"), all.x = TRUE)
shift_share_data[is.na(national_inflow), national_inflow := 0]
shift_share_data[is.na(raw_count),       raw_count       := 0]
shift_share_data[, predicted_inflow := share * national_inflow]

# UNDIVIDED variant (added 2026-09-29). raw_count is the pooled respondent
# count before the eligible-waves division, so weighting it by the 1990 shares
# gives the shift-share instrument with no wave normalisation.
#
# This MUST be constructed here, at MSA x arrival year, rather than recovered
# downstream: 07_exposure_iv.R averages over each individual's ages-15-25
# window, and n_waves differs across the years in that window, so no single
# factor applied to the averaged variable can reverse the division. Analysis
# scripts 24-28 each reversed it at this level in-script; emitting it from the
# pipeline retires those duplicate builds.
shift_share_data[, predicted_inflow_nodiv := share * raw_count]

# Aggregate to MSA x YEAR level
msa_year_predicted_inflows <- shift_share_data[,
                                               .(total_predicted_inflow       = sum(predicted_inflow,       na.rm = TRUE),
                                                 total_predicted_inflow_nodiv = sum(predicted_inflow_nodiv, na.rm = TRUE)),
                                               by = .(MSA, YEAR)]

# Actual inflows by year, MSA, and origin (BPLD) -- same raw, pooled-across-
# waves, divided-by-eligible-waves logic as national_inflows above, broken
# out by each respondent's MSA of residence at the time they were surveyed.
# This is the only geography IPUMS/ACS provides; there is no "MSA at time
# of arrival" variable. Reuses arrivals/eligible_waves built above.
actual_inflows_raw <- arrivals[, .(raw_count = .N), by = .(YEAR = YRIMMIG, MSA, BPLD)]
actual_inflows_raw <- merge(actual_inflows_raw, eligible_waves,
                            by.x = "YEAR", by.y = "YRIMMIG")
actual_inflows_raw[, actual_inflow := raw_count / n_waves]

actual_inflows_msa_year <- actual_inflows_raw[,
                                              .(total_actual_inflow       = sum(actual_inflow, na.rm = TRUE),
                                                total_actual_inflow_nodiv = sum(raw_count,     na.rm = TRUE)),
                                              by = .(MSA, YEAR)]

# Free memory from intermediate objects no longer needed. immigrant_shares
# and national_inflows are also dropped here (and left out of the checkpoint
# below) -- both are fully folded into shift_share_data already and are not
# referenced again downstream.
rm(immigrants_1990, all_msas_1990, all_bplds, msa_bpld_grid,
   valid_pairs, actual_inflows_raw, years_dt,
   immigrant_shares, national_inflows,
   waves, arrivals, eligible_waves)
gc()

# -----------------------------------------------------------------------------

save(reduced_data, selected_msas, top_year_pairs, top175_list_named,
     shift_share_data, msa_year_predicted_inflows, actual_inflows_msa_year,
     file = file.path(local_checkpoint_dir, "after_sec06.RData"))
cat("Checkpoint saved:", file.path(local_checkpoint_dir, "after_sec06.RData"), "\n")
cat("06_shift_share_iv.R complete.\n")
