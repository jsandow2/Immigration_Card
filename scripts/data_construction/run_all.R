# Runs the full data-construction pipeline in order, exactly reproducing what
# the original monolithic Immigration_Card_Data_Construction.R did in one
# session. Safe to re-run from any point: comment out completed sections and
# instead load(file.path(local_checkpoint_dir, "after_secNN.RData")) for the
# last one you finished. Per-section checkpoints live in local_checkpoint_dir
# (set in 00_setup.R, currently C:/R_checkpoints/Immigration_Card) -- not in
# this OneDrive-synced folder.
source("scripts/data_construction/00_setup.R")
source("scripts/data_construction/01_ipums_panel_load.R")
source("scripts/data_construction/02_msa_crosswalk.R")
source("scripts/data_construction/03_variable_construction.R")
source("scripts/data_construction/04_sample_construction.R")
source("scripts/data_construction/05_stable_msa_selection.R")
source("scripts/data_construction/06_shift_share_iv.R")
source("scripts/data_construction/07_exposure_iv.R")
cat("\nData construction pipeline complete. analysis_data.RData is ready for the analysis/ pipeline.\n")
