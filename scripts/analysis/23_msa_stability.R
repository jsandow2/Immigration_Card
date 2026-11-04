# =============================================================================
# 23_msa_stability.R  --  Section 16: Robustness -- MSA stability check
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 17_msa_stability.R.
# Renumbered 2026-10-21 from 16_msa_stability.R (archive-and-renumber pass).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

# CHECKPOINT LOAD REPLACED 2026-09-30. This loaded after_sec15.RData, chaining
# through the tables and figures scripts, though all it needs is
# top175_list_named -- which comes from 05_stable_msa_selection.R and is saved
# into analysis_data.RData directly. 15 is being rewritten, so the chain link
# was both unnecessary and fragile.
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
top175_list_named <- .tmp_env$top175_list_named
all_natives_175   <- .tmp_env$all_natives_175
rm(.tmp_env); gc()
stopifnot(length(top175_list_named) > 1L,
          all(c("MSA", "MSA_Name") %in% names(top175_list_named[[1]])))
cat(sprintf("Loaded top175_list_named from analysis_data.RData (%d years)\n",
            length(top175_list_named)))

# Section 16: Robustness -- MSA stability check
#
# For each consecutive pair of years, identifies which MSAs entered or exited
# the top-175 list. Produces a LaTeX table for inclusion in the appendix.
#
# UNAFFECTED BY THE 2026-09-30 RESPECIFICATION. top175_list_named is built in
# 05_stable_msa_selection.R from PERWT-ranked metro populations by year. The
# 06/07 changes -- the undivided inflow series, the windowed endogenous
# variable, and the two window-construction fixes -- do not touch it, so this
# table's content is identical to before.
#
# One choice worth a sentence in the appendix text: the ranking reported here
# is the PERWT-weighted one. 05 also produces top175_list_named_raw, and the
# note in 06_shift_share_iv.R records that 173 of 175 metros overlap between
# the two, so the choice is immaterial to the selected sample.
# -----------------------------------------------------------------------------

years              <- names(top175_list_named)
comparison_results <- list()

for (year_index in 2:length(years)) {
  prev_year <- years[year_index - 1]
  curr_year <- years[year_index]

  prev_dt <- top175_list_named[[prev_year]]
  curr_dt <- top175_list_named[[curr_year]]

  entered <- curr_dt[!curr_dt$MSA %in% prev_dt$MSA, ]
  exited  <- prev_dt[!prev_dt$MSA %in% curr_dt$MSA, ]

  comparison_results[[paste(prev_year, curr_year, sep = "_vs_")]] <- list(
    same         = setequal(prev_dt$MSA, curr_dt$MSA),
    entered_MSAs = entered,
    exited_MSAs  = exited
  )
}

stability_table <- data.frame(
  Year_Comparison = character(),
  Same            = logical(),
  Entered_MSAs    = character(),
  Exited_MSAs     = character(),
  stringsAsFactors = FALSE
)

for (year_pair in names(comparison_results)) {
  result <- comparison_results[[year_pair]]

  # Separated with "; " rather than ", " (2026-10-19). Metropolitan area names
  # contain their own commas -- "Amarillo, TX" -- so a comma-separated list of
  # them cannot be read, by a person or by code, as a list of names. The count
  # columns below are taken from the vectors directly for the same reason.
  entered_vec <- if (nrow(result$entered_MSAs) > 0) {
    as.character(na.omit(result$entered_MSAs$MSA_Name))
  } else character(0)

  exited_vec <- if (nrow(result$exited_MSAs) > 0) {
    as.character(na.omit(result$exited_MSAs$MSA_Name))
  } else character(0)

  entered_names <- if (length(entered_vec)) {
    paste(entered_vec, collapse = "; ")
  } else "None"

  exited_names <- if (length(exited_vec)) {
    paste(exited_vec, collapse = "; ")
  } else "None"

  stability_table <- rbind(stability_table, data.frame(
    Year_Comparison = year_pair,
    Same            = result$same,
    N_Entered       = length(entered_vec),
    N_Exited        = length(exited_vec),
    Entered_MSAs    = entered_names,
    Exited_MSAs     = exited_names,
    stringsAsFactors = FALSE
  ))
}

# Written to .tex rather than left as a bare expression (2026-09-22). A
# top-level kable() call auto-prints in the console but is silently discarded
# under source(), since source() runs with print.eval = FALSE -- the same
# issue found in 24_results_tables.R. column_spec() widths are required here
# because the entered/exited cells hold up to thirty comma-separated MSA
# names and will otherwise run off the page in longtable.
# REWRITTEN 2026-10-19 to match the table convention used elsewhere in the
# chapter. Three changes:
#
#   1. No "striped". kable_styling(striped) emits \cellcolor{gray!10} on
#      alternating rows, which requires \usepackage[table]{xcolor}. Loading
#      that package is the option clash resolved on 2026-10-04 by taking
#      striping out of 24_results_tables.R; the same reasoning applies here,
#      and this table was the last source of \cellcolor in the pipeline.
#
#   2. The "Same" column is gone. It was FALSE in every row, so it separated
#      nothing -- a column with no variation is a sentence, and it is now one,
#      in the note. In its place are the counts entering and leaving, which
#      are the quantities the text's churn claim actually rests on and which
#      a reader cannot recover from the name lists by counting commas.
#
#   3. Year pairs read 1990--2000 rather than 1990_vs_2000, so no underscore
#      escaping is needed and the column fits in less width.
stability_out <- data.frame(
  `Year pair` = gsub("_vs_", "--", stability_table$Year_Comparison),
  `In / Out`  = sprintf("%d / %d", stability_table$N_Entered,
                        stability_table$N_Exited),
  `Entered`   = stability_table$Entered_MSAs,
  `Exited`    = stability_table$Exited_MSAs,
  check.names = FALSE, stringsAsFactors = FALSE)

tbl_msa_stability <- kable(
  stability_out, format = "latex", booktabs = TRUE, longtable = TRUE,
  align = c("l", "c", "l", "l"), escape = TRUE,
  caption = paste0(
    "Metropolitan Area Entry and Exit in the Top 175, Consecutive Sample ",
    "Years. The stable panel used throughout retains the 175 metropolitan ",
    "areas appearing most often among the largest 175 across sample years, ",
    "so the churn recorded here is in the year-by-year ranking rather than ",
    "in the estimation panel, whose composition is fixed by construction. ",
    "``In / Out'' counts the metropolitan areas entering and leaving the ",
    "top 175 between the two years named. The set differs between every ",
    "consecutive pair of sample years, which is why a per-year top 175 was ",
    "not used. Names are separated by semicolons because metropolitan area ",
    "names contain commas."),
  caption.short = "Metropolitan area entry and exit in the top 175",
  label = "msa_stability") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"),
                font_size = 8) %>%
  column_spec(1, width = "1.8cm") %>%
  column_spec(2, width = "1.2cm") %>%
  column_spec(3, width = "5.8cm") %>%
  column_spec(4, width = "5.8cm")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(tbl_msa_stability),
           "outputs/tables/table_msa_stability.tex")
cat("Table written: outputs/tables/table_msa_stability.tex\n")

# --- The other half of the appendix promise ----------------------------------
# The Data section says the full list of the 175 metropolitan areas is
# reported in the appendix. It was not, because nothing wrote it. It is
# written here rather than in a new script because this script already holds
# top175_list_named and is the one that owns metropolitan-panel composition.
# The stable set is taken from all_natives_175, not from top175_list_named.
# The latter is the ranking year by year, so its union is every metro that
# ever placed in a top 175 -- a larger set than the 175 actually retained.
# all_natives_175 is the estimation sample after the restriction, so the
# metros in it are the stable panel by definition. Names are resolved from
# the per-year lists, which are the only place MSA_Name lives.
name_lookup <- unique(do.call(rbind, lapply(top175_list_named,
                                            function(x) x[, c("MSA", "MSA_Name")])))
stable_msa  <- unique(all_natives_175$MSA)
msa_names   <- sort(unique(na.omit(
  name_lookup$MSA_Name[match(stable_msa, name_lookup$MSA)])))

cat(sprintf("Stable panel: %d metros in all_natives_175, %d resolved to names\n",
            length(stable_msa), length(msa_names)))
if (length(msa_names) != length(stable_msa))
  warning(sprintf("%d metro(s) in the stable panel have no name in the ",
                  length(stable_msa) - length(msa_names)),
          "per-year lists and are omitted from the appendix table.")

# Laid out in three columns down the page so 175 names do not run to many
# pages. Padded to a whole number of rows before the matrix is filled.
ncol_list <- 3L
nrow_list <- ceiling(length(msa_names) / ncol_list)
padded    <- c(msa_names, rep("", nrow_list * ncol_list - length(msa_names)))
msa_grid  <- as.data.frame(matrix(padded, nrow = nrow_list, byrow = FALSE),
                           stringsAsFactors = FALSE)
names(msa_grid) <- rep("Metropolitan area", ncol_list)

tbl_msa_list <- kable(
  msa_grid, format = "latex", booktabs = TRUE, longtable = TRUE,
  align = rep("l", ncol_list), escape = TRUE, row.names = FALSE,
  caption = paste0(
    "The 175 Metropolitan Areas in the Estimation Panel. These are the 175 ",
    "metropolitan areas appearing most frequently among the largest 175 ",
    "across the sample years, held fixed across all waves. Entry and exit ",
    "in the underlying year-by-year ranking is reported in the preceding ",
    "table."),
  caption.short = "The 175 metropolitan areas in the estimation panel",
  label = "msa_list") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"),
                font_size = 8) %>%
  column_spec(1:ncol_list, width = "4.6cm")

writeLines(as.character(tbl_msa_list),
           "outputs/tables/table_msa_list.tex")
cat("Table written: outputs/tables/table_msa_list.tex\n")

# --- Sample-year count: the chapter states two different numbers -------------
# The Data section says 18 sample years; the robustness text says 19. This
# reports what the MSA ranking was actually built over, so the discrepancy is
# settled from the data rather than by picking one.
cat(sprintf("\nTop-175 ranking computed over %d years: %s\n",
            length(top175_list_named),
            paste(names(top175_list_named), collapse = ", ")))

print(stability_table)

# -----------------------------------------------------------------------------

# CHECKPOINT SAVE REMOVED 2026-09-30. after_sec16.RData is orphaned: nothing
# loads it. 10_exposure_missing_check.R was the only consumer and now reads
# analysis_data.RData directly, as do 19, 20, 20b, 29, 31b and this script.
# Writing a checkpoint no reader needs is how the chain acquired links that
# were load-bearing in appearance only.
cat("23_msa_stability.R complete.\n")
