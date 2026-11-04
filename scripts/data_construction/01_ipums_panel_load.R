# =============================================================================
# 01_ipums_panel_load.R  --  Section 1: IPUMS panel data load
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# =============================================================================

source("scripts/data_construction/00_setup.R")

# Section 1: IPUMS panel data load
#
# Loads the IPUMS USA panel extract covering 1990-2023 (~80.1M person-year
# rows; 1990/2000 at 5% Census density, all other years at 1% ACS density;
# 2001-2005 and 2020 are not in the extract -- see 00_setup.R header).
#
# Read via read_ipums_micro_chunked() rather than a single-shot
# read_ipums_micro() call. On this machine, a single-shot read of the full
# ~20.8GB usa_00007.csv reliably hung partway through (endpoint protection
# software appears to interfere with vroom's large memory-mapped reads);
# chunked reading avoided the problem. Confirmed 2026-09-07 that chunking
# introduces no duplication: nrow(all_data) == uniqueN(all_data[, .(YEAR,
# SERIAL, PERNUM)]), and per-year row counts match known Census (5%) and ACS
# (1%) sampling densities.
# -----------------------------------------------------------------------------

ddi <- read_ipums_ddi("usa_00007.xml")

all_data <- read_ipums_micro_chunked(
  ddi,
  callback = IpumsDataFrameCallback$new(function(x, pos) {
    cat("Chunk starting at row", pos, "-", nrow(x), "rows read\n")
    x
  }),
  chunk_size = 1000000
)

setDT(all_data)

save(all_data, file = "all_data.RData")

# -----------------------------------------------------------------------------

save(all_data,
     file = file.path(local_checkpoint_dir, "after_sec01.RData"))
cat("Checkpoint saved:", file.path(local_checkpoint_dir, "after_sec01.RData"), "\n")
cat("01_ipums_panel_load.R complete.\n")
