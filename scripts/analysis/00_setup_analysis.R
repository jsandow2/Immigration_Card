# =============================================================================
# 00_setup_analysis.R -- shared library loads for the analysis-half pipeline.
# Extracted verbatim from the library() block that originally opened Section 9
# of Immigration_Card_Analysis.R. Sourced by every 09-29 script EXCEPT
# 08_load_analysis_subsets.R, which keeps its own inline copy (it is the fresh
# start of the analysis session, mirroring the original script's behavior).
#
# Renumbered 2026-09-07: analysis/ scripts continue the same sequence as
# data_construction/ (which ends at 07), rather than restarting at 01 --
# old 09_load_analysis_subsets.R through 21_dissertation_table_helper.R are
# now 08 through 24. See scripts/README.md for the full mapping.
#
# analysis_checkpoint_dir holds this half's per-section resume checkpoints,
# same local_checkpoint_dir used by data_construction/00_setup.R (a local,
# non-OneDrive-synced folder) but in its own "analysis" subfolder so the two
# halves don't collide writing e.g. after_sec01.RData to the same path.
# Defined here (not sourced from data_construction/00_setup.R, since that
# script also does rm(list=ls()) and library() calls this half doesn't want)
# and duplicated inline in 08_load_analysis_subsets.R, which doesn't source
# this file.
# =============================================================================

local_checkpoint_dir   <- "C:/R_checkpoints/Immigration_Card"
analysis_checkpoint_dir <- file.path(local_checkpoint_dir, "analysis")
dir.create(analysis_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

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

# -----------------------------------------------------------------------------
# bpld_lookup -- ADDED 2026-10-21.
#
# This was built only inside 08_load_analysis_subsets.R and never saved, so it
# existed solely as session state. Any script needing BPLD labels failed in a
# fresh session with "bpld_lookup not found. Run 08 first", and 08 loads the
# full analysis data -- an expensive step to run for a 254-row lookup table.
# It also made the blanket rm(list = ls()) at the top of several scripts
# destructive, which is why those read rm(list = setdiff(ls(), "bpld_lookup")).
#
# Built here instead, so every script that sources this file has it. Cached to
# RDS because read_ipums_ddi parses the whole codebook on each call; the cache
# is derived data and regenerates itself if deleted. 08 keeps its own inline
# copy, as it does for the paths above, because it deliberately does not
# source this file.
# -----------------------------------------------------------------------------
if (!exists("bpld_lookup")) {
  .bpld_cache <- file.path(local_checkpoint_dir, "bpld_lookup.rds")
  bpld_lookup <- if (file.exists(.bpld_cache)) {
    readRDS(.bpld_cache)
  } else {
    .ddi <- read_ipums_ddi(file.path(getwd(), "Codebooks", "usa_00007.xml"))
    .bl  <- as.data.table(ipums_val_labels(.ddi, BPLD))[
              , .(BPLD = as.integer(val), label = lbl)]
    saveRDS(.bl, .bpld_cache)
    cat(sprintf("BPLD lookup: %d codes built from the DDI, cached to %s\n",
                nrow(.bl), .bpld_cache))
    .bl
  }
  rm(list = intersect(ls(), c(".bpld_cache", ".ddi", ".bl")))
}
