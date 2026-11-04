# =============================================================================
# 05_stable_msa_selection.R  --  Section 5: Stable top 175 MSA selection
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# =============================================================================

source("scripts/data_construction/00_setup.R")
load(file.path(local_checkpoint_dir, "after_sec04.RData"))
cat("Loaded checkpoint: after_sec04.RData\n")

# Section 5: Stable top 175 MSA selection
#
# For each year, rank MSAs by population and take the top 175. Select the
# 175 MSAs that appear most frequently across all years in the top-175 list
# (ties broken by total population and average rank). This produces a
# stable set of MSAs consistent across the panel.
#
# METHODOLOGY NOTE (2026-09-11): MSA population is ranked using PERWT-
# weighted sums, not raw respondent counts. This is a separate decision
# from 06_shift_share_iv.R's national_inflow/actual_inflow construction,
# which does use raw counts divided by eligible survey waves -- that
# choice was driven by a problem specific to pooling raw PERWT sums across
# multiple survey waves (differential reweighting; see the methodology
# note at the top of that script). Ranking MSAs by population within a
# single survey year has no multi-wave pooling step, so that problem
# doesn't apply here, and PERWT is arguably the more appropriate tool: it's
# designed precisely to produce accurate standing population estimates for
# a given survey year (via raking to known age/sex/race/geography margins),
# which is exactly what this ranking needs. The raw-count version of this
# selection is also computed and saved below (selected_msas_raw /
# top175_list_named_raw / top_year_pairs_raw), purely so the two Top-175
# lists can be compared directly -- it is not used anywhere downstream.
# The two lists differ by only 2 of 175 MSAs in practice.
# -----------------------------------------------------------------------------

setDT(reduced_data)

msa_crosswalk <- read_excel("data/crosswalks/MSA2013_PUMA2000_pop2000_crosswalk.xls")
setDT(msa_crosswalk)
setnames(msa_crosswalk,
         old = c("MSA Code", "MSA Title"),
         new = c("MSA", "MSA_Name"))
msa_crosswalk_unique <- unique(msa_crosswalk, by = "MSA")
msa_crosswalk_unique[, MSA := as.numeric(MSA)]

top_n        <- 175
use_strict_top <- TRUE

build_stable_top175 <- function(pop_by_year) {
  # pop_by_year: data.table with columns YEAR, MSA, pop

  top175_by_year <- pop_by_year[order(YEAR, -pop),
                                head(.SD, top_n),
                                by = YEAR]
  top175_by_year[, rank := seq_len(.N), by = YEAR]

  freq_table <- top175_by_year[,
    .(top175_count = .N,
      total_pop    = sum(pop),
      avg_rank     = mean(rank)),
    by = MSA]

  freq_table <- freq_table[order(-top175_count, -total_pop, avg_rank)]
  selected <- freq_table[1:top_n, MSA]

  year_pairs <- unique(top175_by_year[, .(YEAR, MSA)])

  top175_by_year_named <- merge(top175_by_year,
                                msa_crosswalk_unique[, .(MSA, MSA_Name)],
                                by = "MSA", all.x = TRUE)
  list_named <- split(
    top175_by_year_named[, .(MSA, MSA_Name)],
    top175_by_year_named$YEAR
  )

  list(selected_msas = selected, top_year_pairs = year_pairs,
       top175_list_named = list_named)
}

# --- Primary: PERWT-weighted population (exclude non-metro MSA == 0) ---
msa_pop_by_year_perwt <- reduced_data[MSA != 0,
                                      .(pop = sum(PERWT, na.rm = TRUE)),
                                      by = .(YEAR, MSA)]
perwt_result <- build_stable_top175(msa_pop_by_year_perwt)

selected_msas      <- perwt_result$selected_msas
top_year_pairs      <- perwt_result$top_year_pairs
top175_list_named   <- perwt_result$top175_list_named

cat("Selected", length(selected_msas), "stable MSAs (PERWT-weighted) across",
    uniqueN(reduced_data$YEAR), "years.\n")

# --- Comparison only: raw respondent counts (not used downstream) ---
msa_pop_by_year_raw <- reduced_data[MSA != 0, .(pop = .N), by = .(YEAR, MSA)]
raw_result <- build_stable_top175(msa_pop_by_year_raw)

selected_msas_raw      <- raw_result$selected_msas
top_year_pairs_raw      <- raw_result$top_year_pairs
top175_list_named_raw   <- raw_result$top175_list_named

cat("Selected", length(selected_msas_raw), "stable MSAs (raw count,",
    "comparison only) across", uniqueN(reduced_data$YEAR), "years.\n")
cat("Overlap between PERWT-weighted and raw-count Top-175:",
    length(intersect(selected_msas, selected_msas_raw)), "of", top_n, "\n")

# -----------------------------------------------------------------------------

save(reduced_data,
     selected_msas, top_year_pairs, top175_list_named,
     selected_msas_raw, top_year_pairs_raw, top175_list_named_raw,
     file = file.path(local_checkpoint_dir, "after_sec05.RData"))
cat("Checkpoint saved:", file.path(local_checkpoint_dir, "after_sec05.RData"), "\n")
cat("05_stable_msa_selection.R complete.\n")
