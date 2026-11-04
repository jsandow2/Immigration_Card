# =============================================================================
# Immigration_Card_Data_Construction.R
# Joe Sandow, University of Nebraska-Lincoln
#
# Full data construction pipeline for:
#   "Does Immigrant Exposure During Formative Years Shape Native
#    Educational Attainment?"
#
# Pipeline order:
#   Section 0  -- Setup and libraries
#   Section 1  -- IPUMS panel data load (all_data.RData via ipumsr)
#   Section 2  -- MSA assignment and crosswalks
#   Section 3  -- Variable construction
#                  (occupation groups, education groups, race, wages,
#                   immigration flag, migration flags)
#   Section 4  -- Sample construction (reduced_data)
#   Section 5  -- Stable top 175 MSA selection
#   Section 6  -- Shift-share IV construction
#   Section 7  -- Individual exposure IV (ages 15-25)
#
# Note (2026-09-07): raw 1990 PUMS ingestion (hand-parsing the 51 fixed-width
# household/person files) is no longer part of this pipeline. The IPUMS panel
# extract loaded in Section 1 already covers 1990 as a native 5% Census
# sample alongside every ACS year, and IPUMS's rectangular extract format
# already attaches household-level variables (MET2013, STATEFIP, etc.) onto
# every person row -- both things the old ingestion step was doing by hand.
# Its script is archived at archive/scripts/01_pums_1990_ingestion_UNUSED.R
# for reference, including a widths_/col_names mismatch bug (the household
# layout was one field short of the true 1990 PUMS Census codebook) found
# and fixed there in case it's ever needed again.
#
# Inputs:
#   usa_00007.xml    -- IPUMS DDI codebook for panel extract (left at root:
#                        paired with a ~20.8GB usa_00007.csv data payload)
#   data/crosswalks/crosswalk_state_puma_msa.csv
#   data/crosswalks/crosswalk_state_pumasupr_msa.csv
#   data/crosswalks/crosswalk_state_countyfip_msa.csv
#   data/crosswalks/MSA2013_PUMA2000_pop2000_crosswalk.xls
#
# Outputs:
#   all_data.RData        -- IPUMS panel with all constructed variables
#                             (~80.1M person-year rows; 1990/2000 at 5%
#                             Census density, 2006-2019/2021-2023 at 1% ACS
#                             density -- 2001-2005 and 2020 are not in the
#                             extract: pre-2006 ACS demonstration years and
#                             the COVID-disrupted 2020 ACS, respectively)
#   reduced_data.RData     -- Analysis sample (all foreign-born +
#                              25% U.S.-born sample by year)
#   analysis_data.RData    -- Final handoff to the analysis/ pipeline
#
# Per-section resume checkpoints (after_secN.RData) are saved to
# local_checkpoint_dir (below) rather than into this OneDrive-synced folder,
# so they don't consume OneDrive storage quota. Only the three cross-stage
# handoff files above -- the ones later sessions or the analysis/ pipeline
# actually need to reach from another device -- stay in this project folder.
# Each section now saves only the specific objects the next section needs
# (via save(...), not save.image()) rather than the whole environment.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 0: Setup
# -----------------------------------------------------------------------------

rm(list = ls())
graphics.off()
gc()

library(tidyverse)
library(data.table)
library(ipumsr)
library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(parallel)

# Local-only (non-OneDrive-synced) folder for per-section resume checkpoints.
# These exist purely to resume an interrupted run on this machine; they don't
# need the cross-device reach that all_data.RData/reduced_data.RData/
# analysis_data.RData do, so keeping them off OneDrive avoids eating into
# OneDrive storage quota.
local_checkpoint_dir <- "C:/R_checkpoints/Immigration_Card"
dir.create(local_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

cat("00_setup.R: libraries loaded, environment cleared.\n")
