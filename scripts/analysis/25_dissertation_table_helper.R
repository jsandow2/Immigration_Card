# =============================================================================
# 25_dissertation_table_helper.R  --  Dissertation table helper (write_dissertation_table + example usage)
# Part of the modular Immigration_Card_Analysis.R pipeline.
# This is the piece you run to produce table_demographic_summary.tex and similar booktabs tables. Depends on demographic_summary_table (from Section 22) and bpld_lookup (from Section 8).
# Renumbered 2026-09-07 from 21_dissertation_table_helper.R.
# Renumbered 2026-10-21 from 24_dissertation_table_helper.R (archive-and-renumber pass).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")
load(file.path(analysis_checkpoint_dir, "after_sec22.RData"))
cat("Loaded checkpoint: after_sec22.RData\n")

# -----------------------------------------------------------------------------
# write_dissertation_table() MOVED 2026-10-04 to 00b_table_helper.R so that
# 24_results_tables.R can use the same house style. 14 runs BEFORE this script,
# so while the function lived here it was unreachable from the one other place
# in the pipeline that needed it, and 14 had been emitting kableExtra-styled
# tables that this file's own header documents as the wrong structure.
#
# The function is unchanged apart from dropping a stray cat(out) that echoed
# the entire LaTeX body to the console on every call. Nothing about the written
# file changes.
# -----------------------------------------------------------------------------
source("scripts/analysis/00b_table_helper.R")

# -----------------------------------------------------------------------------
# Example: Section 22 demographic summary table (top 5 origin groups vs.
# U.S.-born, ages 25+). Run this after Section 22 so demographic_summary_table
# is in the environment.
# -----------------------------------------------------------------------------

demo_table_fmt <- demographic_summary_table %>%
  mutate(
    N             = formatC(N, format = "d", big.mark = ","),
    Pct_female    = sprintf("%.1f", Pct_female),
    Pct_nonHS     = sprintf("%.1f", Pct_nonHS),
    Pct_HS        = sprintf("%.1f", Pct_HS),
    Pct_College   = sprintf("%.1f", Pct_College),
    Pct_Grad      = sprintf("%.1f", Pct_Grad),
    Mean_income   = paste0("\\$", formatC(round(Mean_income), format = "d", big.mark = ",")),
    Median_income = paste0("\\$", formatC(round(Median_income), format = "d", big.mark = ",")),
    Pct_of_immig  = ifelse(is.na(Pct_of_immig), "---", sprintf("%.1f", Pct_of_immig))
  )

colnames(demo_table_fmt) <- c("Group", "N", "\\% Female", "\\% Non-HS", "\\% HS",
                              "\\% College", "\\% Grad", "Mean Income",
                              "Median Income", "\\% of Immig.")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write_dissertation_table(
  demo_table_fmt,
  caption = paste(
    "Demographic Summary Statistics: Top 5 Origin Groups vs.\\",
    "U.S.-Born (Ages 25+)"
  ),
  label = "tab:demographic_summary",
  notes = paste(
    "Top 5 origin groups identified by pooled average population-weighted",
    "count of immigrant arrivals, 2000--2023. Sample restricted to ages 25",
    "and older throughout (see Section~\\ref{sec:age-sensitivity}).",
    "Income figures reflect INCEARN among positive earners only",
    "(INCEARN $>$ 0 and $<$ 999{,}999); INCTOT and INCWAGE are not used",
    "due to sentinel and topcode distortions, respectively.",
    "\\% of Immig.\\ reflects each group's share of all immigrant",
    "observations in the age-restricted sample; not applicable for the",
    "U.S.-born or Combined rows.",
    "The sampling design retains all foreign-born records and a 25\\%",
    "within-year draw of native-born records, so the Combined row averages",
    "the two groups in the design's proportions rather than the population's.",
    "Every other row is a within-group composition and is unaffected."
  ),
  file = "outputs/tables/table_demographic_summary.tex"
)

# -----------------------------------------------------------------------------
# Reuse for any other Section 21/22 table (nativity, gender_educ_nativity,
# MSA volume/share, etc.) by formatting the data.frame the same way and
# calling write_dissertation_table() with an appropriate caption/label/notes.
# -----------------------------------------------------------------------------

# No checkpoint saved here -- this is the terminal script in the pipeline
# and its only output (the .tex file above) is already written to disk.
cat("25_dissertation_table_helper.R complete.\n")
