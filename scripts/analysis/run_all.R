# =============================================================================
# run_all.R  --  Chapter 2 analysis pipeline, in order
#
# REBUILT 2026-10-21 after the archive-and-renumber pass. Every script that
# produces something the paper uses, or that backs a claim the paper makes, is
# sourced below; nothing else remains in this directory. Eleven scripts with
# no artefact and no discussion in the text moved to archive/ on the same day.
#
# NUMBERING. Renumbered 2026-10-21 into run-order blocks:
#   00, 00b  setup and the shared table helper
#   08       data loading
#   09-14    estimation on the main design (09b and 14b are suffixed rather
#            than shifting 15-29: 09b describes the treatment variable, 14b
#            is a decomposition of 14's outcome coding)
#   15-19    shift-share identification and its diagnostics
#   20-21    the Card (2001) replications
#   22-23    descriptives and panel composition
#   24-29    rendering: LaTeX tables and figures
#
# The sequence is continuous with scripts/data_construction/, which runs 01-07,
# so the whole pipeline reads 01 through 29 in order and can be run start to
# finish by someone who did not build it.
# scripts/README.md carries the full old-to-new mapping.
#
# ORDER MATTERS in three places and nowhere else:
#   - 08 must run first. It is the fresh start of the analysis session.
#     (bpld_lookup no longer depends on it: 00_setup_analysis.R builds and
#     caches it, so any script can be run standalone -- see 2026-10-21 note
#     there.)
#   - 24 must run after 09, 10, 11, 14b, 15, 16, 19, 20 and 21, because it renders
#     their CSVs. It reads nothing else and fits nothing.
#   - 25 must run after 22, which writes the checkpoint 25 reads, and 26 must
#     run after both 08 and 22, whose CSVs it renders.
# Everything else loads its own data from analysis_data.RData and could be run
# standalone.
#
# The only surviving checkpoint dependency in Chapter 2 is 25 <- after_sec22,
# written by 22. Every other after_secNN load was removed between 2026-09-30
# and 2026-10-04.
# =============================================================================

source("scripts/analysis/08_load_analysis_subsets.R")   # subsets, extract_composition.csv

# --- Estimation on the main design ------------------------------------------
source("scripts/analysis/09_main_results.R")            # headline; main_results.csv
source("scripts/analysis/09b_exposure_distribution.R")  # treatment distribution; renders its own table
source("scripts/analysis/10_exposure_missing_check.R")  # pooled vs identifying
source("scripts/analysis/11_birth_cohort_final_spec.R") # cohort gradient
source("scripts/analysis/12_final_spec_no_age_fe.R")    # age-FE over-control check
source("scripts/analysis/13_two_margin.R")              # high-school vs college margin
source("scripts/analysis/14_gender_difference_test.R")  # difference by sex
source("scripts/analysis/14b_multinomial_check.R")      # multinomial decomposition

# --- Shift-share identification and diagnostics -----------------------------
source("scripts/analysis/15_rotemberg.R")               # Rotemberg decomposition
source("scripts/analysis/16_adao_placebo.R")            # AKM shock randomisation
source("scripts/analysis/17_asian_placebo_concentration.R")  # resolves 16's Asian failure
source("scripts/analysis/18_rate_spec_loo.R")           # leave-one-out robustness
source("scripts/analysis/19_attribution_stability_check.R")  # metro attribution error
source("scripts/analysis/19b_bhj_incomplete_shares.R")  # BHJ sum-of-shares control; renders its own table

# --- Card (2001) replications -----------------------------------------------
source("scripts/analysis/20_migration_iv.R")            # results_migration_iv.csv
source("scripts/analysis/21_wage_iv.R")                 # results_wage_iv.csv

# --- Descriptives and panel composition -------------------------------------
source("scripts/analysis/22_origin_composition.R")      # origin composition; after_sec22
source("scripts/analysis/23_msa_stability.R")           # MSA churn; both appendix tables

# --- Rendering (must come after every CSV above) ----------------------------
source("scripts/analysis/24_results_tables.R")          # main result and robustness tables
source("scripts/analysis/25_dissertation_table_helper.R")  # demographic summary
source("scripts/analysis/26_appendix_tables.R")         # three appendix tables
source("scripts/analysis/27_education_figures.R")       # Educ_Prop*, Educ_By_Origin
source("scripts/analysis/28_map_origin_world.R")        # origin heat maps
source("scripts/analysis/29_map_destination_us.R")      # destination maps

cat("\nChapter 2 analysis pipeline complete.\n")
