# =============================================================================
# 04_sample_construction.R  --  Section 4: Sample construction
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# =============================================================================

source("scripts/data_construction/00_setup.R")
load(file.path(local_checkpoint_dir, "after_sec03.RData"))
cat("Loaded checkpoint: after_sec03.RData\n")

# Section 4: Sample construction
#
# Analysis sample: all foreign-born + 25% random sample of U.S.-born by year.
# Sampling is reproducible via set.seed(123).
# Output: reduced_data.RData
# -----------------------------------------------------------------------------

set.seed(123)

# --- Thin 1990 and 2000 to match the ACS years' sampling density, and
#     rescale PERWT to compensate ---
#
# 1990 and 2000 are 5% Census samples; every other year in this extract
# (2006-2023) is a 1% ACS sample. Left untouched, 1990/2000 rows dominate
# reduced_data purely because they were drawn 5x denser, not because the
# underlying population differs -- e.g. pre-fix, table(reduced_data$immigr,
# reduced_data$YEAR) showed ~3.89M combined rows for 1990 and ~4.68M for
# 2000, vs ~1M for a typical ACS year.
#
# 20% of a 5% sample is 1%, matching the ACS density. (A version of the
# thinning step also appears, commented out, in
# archive/scripts/Data_clean_all_data.R -- an older script that assumed
# all_data.RData was already thinned from a prior run and only needed to
# build on top of it. Section 1 of this pipeline now rebuilds all_data from
# the raw IPUMS extract on every run, so the thinning has to happen live
# here instead.)
#
# Each surviving 1990/2000 row is rescaled to carry the weight of the ~5
# original rows it now stands in for: PERWT_new = PERWT_old * (original row
# count / thinned row count), computed separately for 1990 and 2000 since
# ceiling() rounding means the true ratio isn't exactly 5.000 in either
# year. Without this rescaling, weighted totals (e.g. sum(PERWT)) for
# 1990/2000 would understate the true population relative to years that
# weren't thinned. As of 2026-09-10, the shift-share/national-inflow
# construction in 06_shift_share_iv.R uses raw respondent counts rather
# than PERWT (see the methodology note atop that script), so this
# rescaling no longer affects that construction -- it's kept because
# PERWT is still used elsewhere in the pipeline (MSA population ranking
# comparison in 05_stable_msa_selection.R, and descriptive statistics in
# the analysis scripts), where accurate population weights still matter.
dense_years <- c(1990, 2000)

thin_counts <- all_data[YEAR %in% dense_years, .(orig_N = .N), by = YEAR]

all_data_thin <- all_data[YEAR %in% dense_years,
                          .SD[sample(.N, size = ceiling(0.2 * .N))],
                          by = YEAR]

thin_counts <- merge(thin_counts,
                     all_data_thin[, .(thinned_N = .N), by = YEAR],
                     by = "YEAR")
thin_counts[, rescale := orig_N / thinned_N]
cat("Thinning summary (1990/2000) -- PERWT rescaled by the exact inverse of\n")
cat("each year's actual thinning ratio:\n")
print(thin_counts)

all_data_thin <- merge(all_data_thin, thin_counts[, .(YEAR, rescale)], by = "YEAR")
all_data_thin[, PERWT := PERWT * rescale]
all_data_thin[, rescale := NULL]

all_data <- rbind(all_data[!YEAR %in% dense_years], all_data_thin)
rm(all_data_thin, thin_counts)

# Population sanity guard: confirms PERWT still sums to a plausible US
# population in every survey wave, right after the 1990/2000 thin-and-
# rescale above and before the native-born 25% subsampling below. This is
# the one point in the pipeline where sum(PERWT) should approximate the
# true national population -- foreign-born rows are kept in full and the
# 1990/2000 rescale above is designed to restore what thinning removed, so
# a failure here means PERWT itself has drifted, not that a downstream
# analytical sampling choice (like the native-born 25% draw just below)
# has changed the implied population, which it deliberately does. Checked
# here rather than in 06_shift_share_iv.R for exactly that reason: by the
# time reduced_data reaches Section 6, only 25% of native-born rows survive
# with unrescaled PERWT, so sum(PERWT) there is expected to run well under
# the true population and isn't a meaningful check on PERWT's own
# correctness.
pop_by_year <- all_data[, .(pop_est = sum(PERWT, na.rm = TRUE)), by = YEAR]
stopifnot(pop_by_year[, all(pop_est > 200e6 & pop_est < 400e6)])
rm(pop_by_year)

foreign_born    <- all_data[immigr == 1]
us_born_sampled <- all_data[immigr == 0,
                             .SD[sample(.N, size = ceiling(0.25 * .N))],
                             by = YEAR]

reduced_data <- rbind(foreign_born, us_born_sampled)

# Validate sample distribution
print(table(reduced_data$immigr, reduced_data$YEAR))

save(reduced_data, file = "reduced_data.RData")

rm(foreign_born, us_born_sampled, all_data)
gc()

# -----------------------------------------------------------------------------

save(reduced_data,
     file = file.path(local_checkpoint_dir, "after_sec04.RData"))
cat("Checkpoint saved:", file.path(local_checkpoint_dir, "after_sec04.RData"), "\n")
cat("04_sample_construction.R complete.\n")
