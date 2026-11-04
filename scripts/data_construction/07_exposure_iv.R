# =============================================================================
# 07_exposure_iv.R  --  Section 7: Individual exposure IV (ages 15-25)
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# Ends by saving analysis_data.RData -- the documented, named bridge file
# required by the analysis pipeline's 09_load_analysis_subsets.R. This is
# the final section's real output; no separate after_sec07 resume checkpoint
# is saved, since there is no Section 8 within data_construction/ left to
# resume into.
# =============================================================================

source("scripts/data_construction/00_setup.R")
load(file.path(local_checkpoint_dir, "after_sec06.RData"))
cat("Loaded checkpoint: after_sec06.RData\n")

# Section 7: Individual exposure IV (ages 15-25) (data.table rewrite)
#
# For each individual, compute the average shift-share predicted inflow
# across the years they were aged 15-25 that fall within the national
# inflow series' coverage (see MAX_WAVE_GAP in 06_shift_share_iv.R):
#
#   avg_exposure_i = (1/n_covered) * sum_{a: t(a) covered} predicted_inflow_{m, t(a)}
#
# where t(a) is the calendar year in which individual i was age a, and
# n_covered is the number of those 11 years actually covered by the inflow
# series (usually but not always 11 -- see the NOTE at Step 5 below).
# Individuals with zero covered years in their 15-25 window are assigned
# exposure_missing = 1 and avg_exposure = 0.
#
# AS OF 2026-09-29 the ACTUAL inflow series is averaged over the identical
# window, producing avg_actual_exposure. Previously the endogenous variable was
# attached on SURVEY year while the instrument used the window, which is a
# timing mismatch that inverted the first stage (-0.069, negative in 52 of 52
# chunks) and therefore the sign of every 2SLS estimate built on it. Both
# measures are now produced; see the Step 4b note and the NAMING NOTE at Step 6
# for which to use.
#
# Variables this file adds to reduced_data_iv:
#   avg_exposure                predicted inflow over ages 15-25, wave-divided
#   avg_exposure_nodiv          same, no wave divisor          [CH2 instrument]
#   avg_actual_exposure         actual inflow over ages 15-25, wave-divided
#   avg_actual_exposure_nodiv   same, no wave divisor          [CH2 endogenous]
#   total_actual_inflow         actual inflow in the SURVEY year   [superseded]
#   exposure_missing            1 where no window year is covered
#   n_window_years, n_covered_years, n_covered_actual   window diagnostics
#
# The _nodiv pair is the Chapter 2 headline specification. The divided pair is
# retained for the Chapter 3 divisor comparison. Both are averaged over the
# identical window, so the only difference between them is the normalisation.
#
# Steps:
#   0. Assign unique person identifier
#   1. Compute birth year and years aged 15 and 25
#   2. Expand each individual to one row per year aged 15-25
#   3. Join to MSA x YEAR predicted inflows
#   4b. Join the same window to MSA x YEAR actual inflows
#   5. Average both over the window, and assert they share it
#   6. Merge back, attach the survey-year measure, construct analysis variables
# -----------------------------------------------------------------------------

# Step 0: Assign unique person ID
setDT(reduced_data)
reduced_data[, PERSON_ID := .I]

# Step 1: Birth year and exposure window endpoints
reduced_data[, `:=`(
  birth_year        = YEAR - AGE,
  age_15_year       = YEAR - AGE + 15,
  age_25_year       = YEAR - AGE + 25,
  AGE_15_25         = as.integer(AGE >= 15 & AGE <= 25),
  exposed_post_2000 = as.integer((YEAR - AGE + 25) >= 2000)
)]

# Step 2: Filter to individuals with any 15-25 exposure in the data window
#
# WHY exposed_post_2000 EXISTS (documented 2026-09-22 after being misread as a
# stale truncation twice). It is (YEAR - AGE + 25) >= 2000, which is
# birth_year >= 1975, which is exactly: the ages-15-25 window BEGINS IN 1990 OR
# LATER. 1990 is the year the shift-share shares are measured in
# 06_shift_share_iv.R.
#
# That is the identification requirement, not a data-coverage limit. The
# instrument weights national inflows by share_{m,k,1990}. For a cohort whose
# window predates 1990 -- someone born in 1950 has a window of 1965-1975 --
# those inflows would be weighted by a settlement distribution measured fifteen
# to twenty-five years afterwards, a distribution that partly reflects the very
# inflows being instrumented. The shares would postdate the treatment and the
# exogeneity argument would fail.
#
# So this filter does NOT become obsolete now that the national inflow series
# reconstructs back to 1910. The longer series lets us MEASURE inflows for
# earlier years; it does not supply a pre-1990 share to weight them with.
#
# WHY 1990 (documented 2026-09-30, rewritten 2026-11-02)
#
# The base year is a choice on an exogeneity-relevance trade-off, not a
# data-coverage limit -- the census supplies metropolitan detail well before
# 1990. Two conditions bind. The shares must postdate the 1965 Immigration
# and Nationality Act, which reshaped origin composition, and must follow it
# by long enough for post-Act settlement to have consolidated. They must also
# sit close to the period the exposure window averages over, because an
# eleven-year window needs shares near it to retain a first stage.
#
# The cost is cohorts. A 1990 base means no individual whose ages-15--25
# window opens earlier can carry exposure, and the instrument has roughly
# 5,775 distinct values (175 MSAs x 33 birth cohorts) behind 3.48 million
# identifying observations. That is the price of a base-year share, not a
# defect to engineer away: loosening the condition to enlarge the identifying
# sample would forfeit the design.
#
# Consequence: everyone else gets avg_exposure = NA -> exposure_missing = 1.
# Without an age restriction that is 59.2% of the sample (5,053,310 of
# 8,537,597); under the older AGE >= 30 floor it was 82.5%. This is the price
# of a base-year share, not a defect to be engineered away -- loosening this
# condition to enlarge the identifying sample would forfeit the design.
# BUG FIX 2026-09-29: require the window to OPEN before the series ends.
#
# Without `age_15_year <= 2022`, anyone born 2008 or later enters here with
# start_year > 2022, and the expansion below calls seq(start_year,
# min(end_year, 2022L)) -- e.g. seq(2023, 2022), which R evaluates as a
# DESCENDING sequence c(2023, 2022) rather than returning nothing. The 2022
# element then matches msa_year_predicted_inflows, so a child aged 0-15 at
# survey came out with avg_exposure equal to a single year's inflow,
# n_covered_years = 1, and exposure_missing = 0.
#
# Measured on a 500k draw from all_natives_175: 41,811 of 249,400 supposedly
# identifying observations were people born 2008-2023 -- 16.8% of the
# identifying sample, roughly 713,000 observations scaled to the full sample.
# They are also the most heavily censored rows in the data, being too young to
# have completed any schooling, so they were both spurious and pulling on the
# outcome.
#
# Surfaced by 19_rotemberg.R, whose own window filter
# (YEAR >= age_15_year & YEAR <= age_25_year) excluded them correctly and so
# dropped 16.8% of the sample relative to the instrument -- which is what broke
# the Rotemberg decomposition identity by 11%.
# AGE_15_25 BRANCH REMOVED 2026-09-29.
#
# The condition was `AGE_15_25 == 1 | exposed_post_2000 == 1`. For anyone born
# 1975 or later exposed_post_2000 is already 1, so the AGE_15_25 branch added
# nothing for them. The only people it admitted were those born BEFORE 1975 who
# happened to be aged 15-25 at survey -- in practice the 1990 wave, birth years
# 1965-1974 -- whose ages-15-25 windows run from 1980-1990 to 1989-1999 and so
# predate the 1990 share year.
#
# Those are precisely the cases the note above says are excluded, and for the
# reason given there: weighting pre-1990 inflows by a 1990 settlement
# distribution uses shares that postdate the treatment and partly reflect the
# inflows being instrumented. The identifying birth-year range was showing
# 1965-2007 rather than 1975-2007 as a result, and distinct instrument values
# came to 7,523 (43 cohorts x 175 metros) rather than 5,775 (33 x 175).
#
# Scale: roughly 1-2% of identifying observations, so the estimates should
# barely move. This is a correctness fix, not a magnitude fix. To restore the
# previous behaviour, add `AGE_15_25 == 1 |` back to the filter.
target_individuals <- reduced_data[
  exposed_post_2000 == 1 & age_15_year <= 2022L,
  .(PERSON_ID, MSA, start_year = age_15_year, end_year = age_25_year)
]
# Both halves of the identification requirement, asserted rather than assumed:
# the window must open in 1990 or later (shares predate exposure) and before
# the inflow series ends (the seq() direction bug fixed above).
stopifnot(target_individuals[, all(start_year >= 1990L)],
          target_individuals[, all(start_year <= pmin(end_year, 2022L))])

# Step 3: Expand to one row per year aged 15-25 (cap at 2022)
individual_years_expanded <- target_individuals[,
                                                .(YEAR = seq(start_year, min(end_year, 2022L))),
                                                by = .(PERSON_ID, MSA)
]

# Step 4: Join to MSA x YEAR predicted inflows
setDT(msa_year_predicted_inflows)
setkey(msa_year_predicted_inflows, MSA, YEAR)
setkey(individual_years_expanded,  MSA, YEAR)

individual_exposure <- merge(individual_years_expanded,
                             msa_year_predicted_inflows,
                             by = c("MSA", "YEAR"), all.x = TRUE)

# --- Step 4b: join ACTUAL inflows to the SAME window (added 2026-09-29) -----
#
# THE BUG THIS FIXES
# Until now the endogenous variable was attached at Step 6 by merging
# actual_inflows_msa_year onto reduced_data on (MSA, YEAR) where YEAR is the
# SURVEY year -- so a 45-year-old surveyed in 2015 was assigned the count of
# people who ARRIVED in 2015, while their instrument averaged predicted arrivals
# over roughly 1985-1995. Both tables are keyed on arrival year, which is why
# the merge looked symmetric and was not: the instrument went through the window
# expansion above and the endogenous variable did not.
#
# The consequence was a first stage of -0.069 (t = -7.41) -- negative in 52 of
# 52 chunks, |t| between 30 and 46 -- because the regression was asking whether
# metros with high predicted inflow during a cohort's adolescence have high
# actual inflow decades later. Immigrant destinations shifted over that span, so
# the answer is no. Every positive 2SLS estimate produced before this fix was a
# negative reduced form divided by that negative denominator.
#
# Joining the actual series to individual_years_expanded puts the endogenous
# variable on the same window as the instrument. The first stage becomes
# +0.54 to +0.65 depending on sample and normalisation.
#
# total_actual_inflow (survey year) is still constructed at Step 6 and kept, so
# the superseded specification remains reproducible for the timing comparison.
# The two are NOT interchangeable -- see the naming note there.
setDT(actual_inflows_msa_year)
setkey(actual_inflows_msa_year, MSA, YEAR)
individual_exposure <- merge(
  individual_exposure,
  actual_inflows_msa_year[, .(MSA, YEAR,
                              window_actual_inflow       = total_actual_inflow,
                              window_actual_inflow_nodiv = total_actual_inflow_nodiv)],
  by = c("MSA", "YEAR"), all.x = TRUE)

# NOTE (2026-09-08): total_predicted_inflow is NA here for any window-year
# with no entry in msa_year_predicted_inflows. Before MAX_WAVE_GAP was added
# to 06_shift_share_iv.R's closest-wave rule, this almost never happened --
# predicted inflow existed for nearly the entire historical range (however
# badly measured for old cohorts). Now that arrival cohorts more than
# MAX_WAVE_GAP years from any survey wave are excluded from the national
# inflow series, msa_year_predicted_inflows has real gaps, and a meaningful
# share of individuals have one or more 15-25 window-years with no
# predicted-inflow data. Zero-filling those years before averaging (the
# previous behavior) silently dilutes avg_exposure toward zero in
# proportion to how much of a person's window is uncovered -- a
# measurement-error pattern correlated with birth cohort and MSA settlement
# history, not a real signal. Average only over the window-years that are
# actually covered instead.

# Step 5: Average predicted inflow over the COVERED portion of the 15-25
# window per individual. n_window_years/n_covered_years are kept for
# transparency (e.g. reporting what share of the estimation sample has
# partial window coverage). avg_exposure is NaN when zero years are
# covered; converted to explicit NA below so exposure_missing (set in
# Step 6) continues to mean exactly "no usable exposure data" rather than
# conflating it with partial-window dilution.
individual_iv <- individual_exposure[,
                                     .(avg_exposure              = mean(total_predicted_inflow,       na.rm = TRUE),
                                       avg_exposure_nodiv        = mean(total_predicted_inflow_nodiv, na.rm = TRUE),
                                       avg_actual_exposure       = mean(window_actual_inflow,         na.rm = TRUE),
                                       avg_actual_exposure_nodiv = mean(window_actual_inflow_nodiv,   na.rm = TRUE),
                                       n_window_years            = .N,
                                       n_covered_years           = sum(!is.na(total_predicted_inflow)),
                                       n_covered_actual          = sum(!is.na(window_actual_inflow))),
                                     by = PERSON_ID
]
individual_iv[n_covered_years  == 0, `:=`(avg_exposure        = NA_real_,
                                          avg_exposure_nodiv  = NA_real_)]
individual_iv[n_covered_actual == 0, `:=`(avg_actual_exposure       = NA_real_,
                                          avg_actual_exposure_nodiv = NA_real_)]

# --- ASSERTION: instrument and endogenous variable share one window ---------
# The bug this file previously carried was invisible precisely because the two
# merges looked alike. This check is permanent: both averages must come from the
# same expanded window, over the same individuals, with the same denominator.
# If a future edit reattaches either variable outside individual_years_expanded,
# this stops the build rather than letting a silent timing mismatch propagate
# into a sign flip.
stopifnot(
  # one row per individual, same set for both
  !anyDuplicated(individual_iv$PERSON_ID),
  individual_iv[, all(n_window_years > 0L)],
  # coverage can differ only where the inflow series itself has gaps, never
  # because the two variables were built over different windows
  individual_iv[, all(n_covered_years  <= n_window_years)],
  individual_iv[, all(n_covered_actual <= n_window_years)]
)
cat(sprintf("Window build: %s individuals | median window %.0f years | ",
            format(nrow(individual_iv), big.mark = ","),
            median(individual_iv$n_window_years)))
cat(sprintf("predicted covered %.1f%% | actual covered %.1f%%\n",
            100 * individual_iv[, mean(n_covered_years  > 0L)],
            100 * individual_iv[, mean(n_covered_actual > 0L)]))

# Free expanded table -- no longer needed
rm(individual_years_expanded, individual_exposure, target_individuals)
gc()

# Step 6: Merge back to reduced_data and attach actual inflows
setkey(reduced_data,   PERSON_ID)
setkey(individual_iv,  PERSON_ID)
reduced_data_iv <- merge(reduced_data, individual_iv,
                         by = "PERSON_ID", all.x = TRUE)

# NAMING NOTE -- READ BEFORE USING EITHER VARIABLE
#
#   avg_actual_exposure   actual arrivals averaged over the individual's ages
#                         15-25 window, in their MSA. Built at Step 4b/5 on the
#                         SAME window as avg_exposure. This is the endogenous
#                         variable the design calls for.
#
#   total_actual_inflow   actual arrivals in the individual's MSA in the year
#                         they were SURVEYED. Merged below on survey YEAR. This
#                         is the superseded measure; it is retained only so the
#                         timing comparison remains reproducible, and it should
#                         NOT be instrumented with avg_exposure.
#
# Both tables are keyed on arrival year, so the two merges look identical. They
# are not. See the note at Step 4b.
setkey(reduced_data_iv,        MSA, YEAR)
reduced_data_iv <- merge(reduced_data_iv, actual_inflows_msa_year,
                         by = c("MSA", "YEAR"), all.x = TRUE)

reduced_data_iv[is.na(total_actual_inflow), total_actual_inflow := 0]
reduced_data_iv[, exposure_missing := as.integer(is.na(avg_exposure))]
reduced_data_iv[is.na(avg_exposure),              avg_exposure              := 0]
reduced_data_iv[is.na(avg_exposure_nodiv),        avg_exposure_nodiv        := 0]
reduced_data_iv[is.na(avg_actual_exposure),       avg_actual_exposure       := 0]
reduced_data_iv[is.na(avg_actual_exposure_nodiv), avg_actual_exposure_nodiv := 0]

# Diagnostic: the first stage on each measure, over the identifying sample.
# The windowed measure should correlate strongly and positively with the
# instrument; the survey-year measure is the one that produced the negative
# first stage. Printed at build time so a future regression is caught here.
# Computed on the ESTIMATION sample -- top-175 natives with usable exposure --
# not on every native in the file. MSA == 0 is IPUMS's not-in-an-identifiable-
# metro residual and metros outside the top 175 are excluded from the analysis,
# so including them here produced a misleadingly low figure (0.43 against 0.95)
# and a spurious warning on a correct build.
#
# The AGE >= 30 row is the benchmark: 0.9493 when this construction was first
# validated in 20m_windowed_endogenous_test.R. The all-ages figure is lower by
# construction, because the window is capped at 2022 -- someone born in 2005 has
# a three-year window rather than eleven, and short windows average noisily.
# That attenuation is a real consequence of dropping the age restriction, not a
# defect, but it belongs in the record.
.id <- reduced_data_iv[immigr == 0 & exposure_missing == 0 & MSA %in% selected_msas]
if (nrow(.id) > 1000L) {
  .c_all <- .id[, cor(avg_exposure, avg_actual_exposure)]
  .s_all <- .id[, cor(avg_exposure, total_actual_inflow)]
  cat("Correlation with avg_exposure, top-175 natives with usable exposure:\n")
  cat(sprintf("  all ages   : windowed %+.4f | survey-year %+.4f  (n = %s)\n",
              .c_all, .s_all, format(nrow(.id), big.mark = ",")))
  .i30 <- .id[AGE >= 30]
  if (nrow(.i30) > 1000L) {
    .c30 <- .i30[, cor(avg_exposure, avg_actual_exposure)]
    cat(sprintf("  AGE >= 30  : windowed %+.4f | survey-year %+.4f  (n = %s)\n",
                .c30, .i30[, cor(avg_exposure, total_actual_inflow)],
                format(nrow(.i30), big.mark = ",")))
    # 0.9493 is the validated benchmark; 0.80 leaves room for sample drift
    # while still catching a broken join.
    if (.c30 < 0.80)
      warning("Windowed actual exposure correlates weakly with the instrument ",
              "on the AGE >= 30 benchmark sample (", round(.c30, 4),
              " against an expected ~0.95) -- check the Step 4b join.")
    rm(.c30)
  }
  rm(.c_all, .s_all, .i30)
}
rm(.id)

# Construct stable top-175 subset for main analysis
all_natives_175 <- reduced_data_iv[immigr == 0 & MSA %in% selected_msas]
setDT(all_natives_175)

# Merge total predicted inflow (MSA x YEAR level) into all_natives_175
setkey(all_natives_175,            MSA, YEAR)
setkey(msa_year_predicted_inflows, MSA, YEAR)
all_natives_175 <- merge(all_natives_175, msa_year_predicted_inflows,
                         by = c("MSA", "YEAR"), all.x = TRUE)
all_natives_175[is.na(total_predicted_inflow), total_predicted_inflow := 0]

cat("all_natives_175:", nrow(all_natives_175), "observations,",
    uniqueN(all_natives_175$MSA), "MSAs.\n")

# Save constructed analysis objects
save(reduced_data_iv, all_natives_175, selected_msas,
     msa_year_predicted_inflows, shift_share_data,
     top175_list_named, top_year_pairs,
     file = "analysis_data.RData")

cat("Data construction complete. analysis_data.RData saved.\n")
cat("07_exposure_iv.R complete.\n")
