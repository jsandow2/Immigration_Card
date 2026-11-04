# =============================================================================
# Immigration_Card_Analysis.R
# Kept its number in the 2026-10-21 renumbering: analysis continues the
# data_construction sequence, which ends at 07.
# Joe Sandow, University of Nebraska-Lincoln
#
# Full analysis pipeline for:
#   "Does Immigrant Exposure During Formative Years Shape Native
#    Educational Attainment?"
#
# Requires: analysis_data.RData produced by
#           scripts/data_construction/07_exposure_iv.R
#
# Pipeline order (renumbered 2026-09-07 to continue data_construction/'s
# 01-07 sequence, rather than restarting at 01 -- see scripts/README.md):
#   08 -- Load data and construct analysis subsets
#   09 -- Education models: run_models() function
#          (Binary IV/LPM, Binary GLM 2SRI, Multinomial 2SRI)
#   10 -- Meta-analytic chunk loop
#          Toggle: RUN_MULTINOMIAL controls multinomial + bootstrap
#          Four Binary IV models (Educ_Binary, HighSch, College, Grad)
#          each with MSA-clustered and state-clustered HC1 SEs
#   11 -- Meta-analytic pooling via metafor::rma()
#   12 -- Individual-level migration models (Card 2001 replication
#          + diff_iv novel instrument)
#   13 -- Wage IV models by occupation group
#   14 -- Results tables (LaTeX via kableExtra)
#   15 -- Education composition figures (Educ_Prop*.png): original,
#          age at arrival >= 25, and age at arrival 18-24
#   16 -- Robustness: MSA stability check
#   17 -- Robustness: exposure_missing subsample check
#   18 -- Robustness: birth-cohort FE
#   19 -- Rotemberg weight decomposition (GPS 2020)
#   20 -- Adao et al. (2019) placebo test
#   21 -- Descriptive stats: Mexican vs Other immigrants
#          (education, gender, occupation, MSA volume/share,
#          named via official 2013 OMB CBSA delineation file)
#   22 -- Demographic summary table: top 5 origin groups
#          (by pooled average population-weighted count, 2000-2023)
#          vs. All Immigrants, U.S.-born, and Full sample
#   23 -- Estimation sample summary statistics (Table 1)
#   24 -- Dissertation table helper (write_dissertation_table())
#
# Outputs:
#   chunk_outputs/                        -- chunk .fst files (completion flags)
#   meta_outputs_fst/
#     meta_first_stage.fst               -- first-stage F per chunk
#     meta_binary_iv_tidy.fst            -- Binary IV Educ_Binary MSA-clustered
#     meta_binary_iv_state_tidy.fst      -- Binary IV Educ_Binary state-clustered
#     meta_iv_hs_msa_tidy.fst            -- Binary IV HighSch MSA-clustered
#     meta_iv_hs_state_tidy.fst          -- Binary IV HighSch state-clustered
#     meta_iv_college_msa_tidy.fst       -- Binary IV College MSA-clustered
#     meta_iv_college_state_tidy.fst     -- Binary IV College state-clustered
#     meta_iv_grad_msa_tidy.fst          -- Binary IV Grad MSA-clustered
#     meta_iv_grad_state_tidy.fst        -- Binary IV Grad state-clustered
#     meta_binary_glm_tidy.fst           -- Binary GLM (2SRI)
#     meta_multinom_tidy.fst             -- Multinomial (2SRI)
#     meta_boot_se_tidy.fst              -- Bootstrap MSA-clustered SEs
#     meta_analysis/
#       combined_meta_summary_full.fst   -- comprehensive pooled results
#   results_educ_full.csv               -- comprehensive meta-analytic results
#   results_multinom_boot.csv           -- multinomial bootstrapped SE results
#   Educ_Prop.png                       -- education composition (unrestricted)
#   Educ_Prop_Age25plus.png             -- education composition (age at arrival >= 25)
#   Educ_Prop_Age18to24.png             -- education composition (age at arrival 18-24)
#   model_summaries.txt                 -- model diagnostics
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# 08_load_analysis_subsets.R  --  Section 8: Load data and construct analysis subsets
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Fresh session start (matches original script). Requires analysis_data.RData
# produced by 07_exposure_iv.R in the data_construction/ pipeline. Renumbered
# 2026-09-07 from 09_load_analysis_subsets.R to continue data_construction/'s
# 01-07 sequence instead of restarting analysis/ at 01.
# =============================================================================

# Section 8: Load data and construct analysis subsets
# -----------------------------------------------------------------------------

rm(list = ls())
graphics.off()
gc()

library(tidyverse)
library(data.table)
library(AER)
library(nnet)
library(sandwich)
library(lmtest)
library(broom)
library(boot)
library(metafor)
library(knitr)
library(kableExtra)
library(modelsummary)
library(ggplot2)
library(fst)
library(fixest)
library(ipumsr)

# Per-section resume checkpoints for this half of the pipeline live in their
# own subfolder of the same local, non-OneDrive-synced local_checkpoint_dir
# data_construction/ uses -- set directly here (not sourced from
# data_construction/00_setup.R, which also does rm(list=ls()) and library()
# calls this fresh session doesn't want) and duplicated in 00_setup_analysis.R
# for every other script in this half.
local_checkpoint_dir    <- "C:/R_checkpoints/Immigration_Card"
analysis_checkpoint_dir <- file.path(local_checkpoint_dir, "analysis")
dir.create(analysis_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

ddi         <- read_ipums_ddi(file.path(getwd(), "Codebooks", "usa_00007.xml"))
bpld_lookup <- as.data.table(ipums_val_labels(ddi, BPLD))[
  , .(BPLD = as.integer(val), label = lbl)]
cat(sprintf("BPLD lookup: %d codes loaded from codebook\n", nrow(bpld_lookup)))

load("analysis_data.RData")
# Objects loaded:
#   reduced_data_iv            -- full analysis sample with IV variables
#   all_natives_175            -- U.S.-born in stable top-175 MSAs
#   selected_msas              -- vector of 175 stable MSA codes
#   msa_year_predicted_inflows -- MSA x YEAR predicted inflows
#   shift_share_data           -- full MSA x BPLD x YEAR shift-share data
#   top175_list_named          -- named list for MSA stability check
#   top_year_pairs             -- year-specific top-175 MSA pairs

# -----------------------------------------------------------------------------
# ESTIMATION-SAMPLE AGE FLOOR -- ADDED 2026-09-11, REMOVED 2026-10-02
#
# `ESTIMATION_MIN_AGE <- 30L` used to be applied here, cutting all_natives_175
# from 8,537,597 to 5,145,128 (60.3% retained). It is gone. The final Chapter 2
# specification carries NO age restriction, and leaving the floor in this
# script meant the pipeline's shared data-construction step silently disagreed
# with the specification every analysis script runs.
#
# It survived the 2026-09-29/30 respecification only because every rebuilt
# script defensively reloads all_natives_175 from analysis_data.RData, which
# holds the unfloored sample, overwriting whatever this script left in the
# session. So the headline results were never on the floored data -- 29 reports
# 8,537,597 -- but the console still announced a floor, and anything written
# later that trusted this script's output would have inherited a restriction
# the specification does not have. 12 and 13 both carry comments explaining
# that they reload precisely to escape it; that reason is now obsolete.
#
# WHY IT WAS THERE, AND WHY THAT REASONING NO LONGER APPLIES.
# The original argument was that 07 admitted people whose ages-15-25 window was
# still in progress at survey, so treatment and outcome were both incomplete
# for the same individuals. On the 100k matched sample, 8,427 of 32,528
# identifying observations (25.9%) were under 18 -- structurally 0 on
# Educ_Binary, y_College and y_Grad -- and their windows sat on 2017-2022,
# exactly the arrival years the eligible-waves divisor scaled up most, so the
# instrument's largest values landed on mechanically-zero outcomes.
#
# Three things dissolved that argument:
#   1. The `AGE_15_25` branch that admitted those observations was removed from
#      07 on 2026-09-29 as a bug. The in-progress windows are gone at source.
#   2. The eligible-waves divisor is not in the final specification, so the
#      scaling that concentrated the instrument on recent arrival years is not
#      being applied either.
#   3. Incomplete schooling among young cohorts is real but is handled by
#      MEASURING it rather than by filtering: 18b documents the cohort gradient
#      and 17 the pooled-versus-identifying gap. A floor hides it instead.
#
# Robustness at other floors, if it is ever wanted, belongs in a dedicated
# check that states the floor it uses -- not in the shared construction step.
# -----------------------------------------------------------------------------
setDT(all_natives_175)
cat(sprintf("all_natives_175: %s observations (no age floor)\n",
            format(nrow(all_natives_175), big.mark = ",")))
cat(sprintf("  exposure_missing share: %.1f%%\n",
            100 * mean(all_natives_175$exposure_missing)))

# -----------------------------------------------------------------------------
# Composition of the extract before any age or metro restriction.
# ADDED 2026-10-19. The Data section reports the immigrant and U.S.-born counts
# of the full analysis extract and their shares. Nothing in the pipeline
# computed them -- they had been carried in the text by hand, so no run could
# confirm or refresh them. The split follows directly from the sampling design
# described there (all immigrants retained, a 25% draw of the U.S.-born), and
# is the denominator for the age-restricted figures 21 reports.
# -----------------------------------------------------------------------------
extract_split <- reduced_data_iv[, .(obs = .N), by = .(immigr)]
extract_split[, group := fifelse(immigr == 1L, "Immigrant", "U.S.-born")]
extract_split[, share := 100 * obs / sum(obs)]
setorder(extract_split, -obs)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(extract_split[, .(group, obs, share)],
          "outputs/tables/extract_composition.csv", row.names = FALSE)

cat("\nAnalysis extract before age/metro restriction:\n")
for (i in seq_len(nrow(extract_split)))
  cat(sprintf("  %-10s %14s  (%.1f%%)\n", extract_split$group[i],
              format(extract_split$obs[i], big.mark = ","),
              extract_split$share[i]))
cat(sprintf("  %-10s %14s\n", "Total",
            format(sum(extract_split$obs), big.mark = ",")))
cat("  Written: outputs/tables/extract_composition.csv\n")

# Ensure Educ_Group reference level is non_HS throughout
all_natives_175$Educ_Group <- relevel(
  factor(all_natives_175$Educ_Group,
         levels = c("non_HS", "HighSch", "College", "Grad")),
  ref = "non_HS"
)

# Factor YEAR, MSA, and all covariates; construct binary outcomes.
# Education-level binary outcomes use non_HS as the 0 baseline throughout,
# consistent with the multinomial 2SRI reference category.
all_natives_175 <- all_natives_175 %>%
  mutate(
    YEAR          = factor(YEAR),
    MSA           = factor(MSA),
    SEX           = factor(SEX),
    Educ_Binary   = as.integer(Educ_Group %in% c("College", "Grad")),
    # non_HS vs higher: the HS-completion margin, as opposed to
    # Educ_Binary's College/Grad-vs-everything-below college margin
    y_HSplus      = as.integer(Educ_Group != "non_HS"),
    y_HS          = as.integer(Educ_Group == "HighSch"),
    y_College     = as.integer(Educ_Group == "College"),
    y_Grad        = as.integer(Educ_Group == "Grad"),
    Marital       = factor(Marital),
    ChildAgeGroup = factor(ChildAgeGroup),
    Mom_Educ      = factor(Mom_Educ),
    Dad_Educ      = factor(Dad_Educ)
  )

# Define demographic subsets from all_natives_175
native_men_175   <- all_natives_175 %>% filter(SEX == 1)
native_women_175 <- all_natives_175 %>% filter(SEX == 2)

# Random sample of 100,000 for computationally intensive diagnostic models
set.seed(123)
sample_dt <- all_natives_175[sample(nrow(all_natives_175), 100000), ]

all_natives  <- sample_dt
native_men   <- all_natives %>% filter(SEX == 1)
native_women <- all_natives %>% filter(SEX == 2)

cat("Subsets constructed:\n")
cat("  all_natives_175:", nrow(all_natives_175), "obs\n")
cat("  all_natives (sample):", nrow(all_natives), "obs\n")
cat("  native_men:", nrow(native_men), "obs\n")
cat("  native_women:", nrow(native_women), "obs\n")
# NO CHECKPOINT WRITTEN, as of 2026-10-04.
#
# This wrote after_sec08.RData, which nothing reads. The checkpoint chain was
# dismantled between 2026-09-30 and 2026-10-04 -- every surviving script
# loads analysis_data.RData directly -- and after that the save was writing
# large objects to disk on every run for no consumer. Worse, a stale
# checkpoint left on disk is precisely how 15_rotemberg.R came to decompose
# a superseded estimate in September. The only checkpoint still in the
# Chapter 2 pipeline is after_sec22, written by 21 and read by 24.

cat("\n08_load_analysis_subsets.R complete (no checkpoint).\n")
