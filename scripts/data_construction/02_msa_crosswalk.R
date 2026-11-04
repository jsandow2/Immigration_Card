# =============================================================================
# 02_msa_crosswalk.R  --  Section 2: MSA assignment and crosswalks
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# =============================================================================

source("scripts/data_construction/00_setup.R")
load(file.path(local_checkpoint_dir, "after_sec01.RData"))
cat("Loaded checkpoint: after_sec01.RData\n")

# Section 2: MSA assignment and crosswalks
#
# Assigns each observation to an MSA using a priority hierarchy:
#   1. MET2013 (directly reported in IPUMS)
#   2. MET2013_crosswalk (derived from PUMA-to-MSA crosswalk)
# Migration origin MSAs (MIGMSA1, MIGMSA5) are assigned via three crosswalks:
#   PUMA, PUMASUPR, and county FIPS -- coalesced in that order.
#
# Each crosswalk is copied and its output column renamed BEFORE merging to
# prevent data.table duplicate column errors. The county crosswalk is reused
# for both MIGMSA1 and MIGMSA5 so it is copied separately each time.
# -----------------------------------------------------------------------------

setDT(all_data)

# -----------------------------------------------------------------------------
# Construct MET2013_crosswalk (added 2026-09-07)
#
# MET2013 (IPUMS-reported metro area) is missing for a meaningful share of
# observations. MET2013_crosswalk is a fallback derived from three external
# crosswalk files, ported from the old standalone archive/scripts/
# Data_clean_all_data.R -- which had assumed all_data.RData already had this
# column baked in from a prior run, rather than building it fresh. Now that
# Section 1 rebuilds all_data from the raw IPUMS extract every time, this
# has to be constructed here instead. (Dropped from the ported logic: a live
# fetch of 90mfips.txt from the Census Bureau's website, which the original
# script computed but never actually used in building MET2013_crosswalk.)
#
# Priority order: 2000-PUMA-based MSA crosswalk (via PUMASUPR) > county-based
# CBSA crosswalk (via COUNTYFIP) > 1990-specific PUMA-to-county fallback
# (via PUMA, since PUMASUPR is only defined post-1990) > assign rural (0)
# for any 1990 observation still unmatched after all three.
# -----------------------------------------------------------------------------

msa_crosswalk_2013 <- as.data.table(read_excel("data/crosswalks/MSA2013_PUMA2000_pop2000_crosswalk.xls"))
cbsa_crosswalk_raw  <- as.data.table(read_excel("data/crosswalks/list1.xls"))
puma_1990_5pct      <- as.data.table(read_excel("data/crosswalks/1990_PUMAs_5pct.xls"))

# Clean 2000-PUMA-based MSA crosswalk (keep the highest-population-share MSA
# per STATEFIP x PUMA2000, since a PUMA can span more than one MSA)
msa_crosswalk_clean <- msa_crosswalk_2013[
  , .(
    STATEFIP = as.integer(`State FIPS Code`),
    PUMA2000 = as.integer(`2000 PUMA Code`),
    MSA_Code = `MSA Code`,
    Percent_PUMA = `Percent PUMA Population`
  )
][
  , .SD[which.max(Percent_PUMA)], by = .(STATEFIP, PUMA2000)
]

all_data <- merge(
  all_data,
  msa_crosswalk_clean,
  by.x = c("STATEFIP", "PUMASUPR"),
  by.y = c("STATEFIP", "PUMA2000"),
  all.x = TRUE
)

# Clean CBSA (county-based) crosswalk
cbsa_crosswalk_clean <- cbsa_crosswalk_raw[
  , .(
    STATEFIP = as.integer(`FIPS State Code`),
    COUNTYFIP = as.integer(`FIPS County Code`),
    CBSA_Code = `CBSA Code`
  )
]

all_data <- merge(
  all_data,
  cbsa_crosswalk_clean,
  by = c("STATEFIP", "COUNTYFIP"),
  all.x = TRUE
)

all_data[, MET2013_crosswalk := fifelse(!is.na(MSA_Code), MSA_Code, CBSA_Code)]
all_data[, MET2013_source := fcase(
  !is.na(MSA_Code), "Super PUMA-based",
  !is.na(CBSA_Code), "County-based",
  default = "Unmatched"
)]

# Build PUMA-to-county mapping from the 1990 5% sample (PUMASUPR doesn't
# exist for 1990, so this is the fallback for that year specifically)
puma_to_county <- puma_1990_5pct[
  !is.na(`County FIPS code`),
  .(
    STATEFIP = as.integer(State),
    PUMA_for_merge = as.integer(PUMA),
    COUNTYFIP = as.integer(`County FIPS code`),
    POPULATION = `Area population`
  )
][
  , .SD[which.max(POPULATION)], by = .(STATEFIP, PUMA_for_merge)
]

puma_to_cbsa <- merge(
  puma_to_county,
  cbsa_crosswalk_clean,
  by = c("STATEFIP", "COUNTYFIP"),
  all.x = TRUE
)
setnames(puma_to_cbsa, "CBSA_Code", "CBSA_Code_fallback")

all_data[, PUMA_for_merge := fifelse(YEAR == 1990 & is.na(PUMASUPR),
                                      as.integer(PUMA), as.integer(PUMASUPR))]

all_data <- merge(
  all_data,
  puma_to_cbsa[, .(STATEFIP, PUMA_for_merge, CBSA_Code_fallback)],
  by = c("STATEFIP", "PUMA_for_merge"),
  all.x = TRUE
)

all_data[, MET2013_crosswalk := fifelse(is.na(MET2013_crosswalk), CBSA_Code_fallback, MET2013_crosswalk)]
all_data[, MET2013_source := fcase(
  !is.na(MSA_Code), "Super PUMA-based",
  is.na(MSA_Code) & !is.na(CBSA_Code), "County-based",
  is.na(MET2013_crosswalk) & !is.na(CBSA_Code_fallback), "PUMA-to-county fallback",
  default = "Unmatched"
)]

all_data[, MET2013_crosswalk := as.numeric(MET2013_crosswalk)]

# Any 1990 observation still unmatched after all three crosswalk passes is
# assigned rural (MSA = 0) rather than left NA
all_data[YEAR == 1990 & is.na(MET2013_crosswalk), MET2013_crosswalk := 0]
all_data[YEAR == 1990 & MET2013_crosswalk == 0, MET2013_source := "Assigned rural (MET2013 = 0)"]

cat("MET2013_crosswalk source distribution:\n")
print(all_data[, .N, by = MET2013_source])
cat("Unmatched MET2013_crosswalk (excluding assigned-rural 1990 cases):",
    sum(is.na(all_data$MET2013_crosswalk)), "\n")

# Drop merge-only helper columns -- MET2013_crosswalk and MET2013_source are
# the only outputs of this block that anything downstream needs
all_data[, c("MSA_Code", "Percent_PUMA", "CBSA_Code",
             "PUMA_for_merge", "CBSA_Code_fallback") := NULL]
rm(msa_crosswalk_2013, cbsa_crosswalk_raw, puma_1990_5pct,
   msa_crosswalk_clean, cbsa_crosswalk_clean, puma_to_county, puma_to_cbsa)
gc()

# -----------------------------------------------------------------------------

# Primary MSA assignment
# NOTE: MET2013 wraps both as.numeric() -- ipumsr variables often carry a
# haven_labelled class alongside their numeric value (for IPUMS's value-
# label metadata), which fifelse() treats as a different class than plain
# numeric even when the underlying values are just numbers.
all_data[, MSA := fifelse(!is.na(MET2013), as.numeric(MET2013), as.numeric(MET2013_crosswalk))]

# Load and deduplicate base crosswalks
puma_crosswalk     <- fread("data/crosswalks/crosswalk_state_puma_msa.csv")
pumasupr_crosswalk <- fread("data/crosswalks/crosswalk_state_pumasupr_msa.csv")
county_crosswalk   <- fread("data/crosswalks/crosswalk_state_countyfip_msa.csv")

setnames(puma_crosswalk,
         c("STATEFIP", "PUMA",     "MSA"),
         c("MIGPLAC",  "MIGPUMA",  "MSA_puma"))
setnames(pumasupr_crosswalk,
         c("STATEFIP", "PUMASUPR", "MSA"),
         c("MIGPLAC",  "MIGPUMA",  "MSA_pumasupr"))
setnames(county_crosswalk,
         c("STATEFIP", "COUNTYFIP", "MSA"),
         c("STATEFIP",  "COUNTYFIP", "MSA_county"))

puma_crosswalk     <- unique(puma_crosswalk,     by = c("MIGPLAC", "MIGPUMA"))
pumasupr_crosswalk <- unique(pumasupr_crosswalk, by = c("MIGPLAC", "MIGPUMA"))
county_crosswalk   <- unique(county_crosswalk,   by = c("STATEFIP", "COUNTYFIP"))

# --- MIGMSA1: 1-year migration origin ---
# Copy each crosswalk and rename its output column before merging so that
# data.table never sees two columns with the same name in the result.
puma_crosswalk_1     <- copy(puma_crosswalk)
pumasupr_crosswalk_1 <- copy(pumasupr_crosswalk)
county_crosswalk_1   <- copy(county_crosswalk)

setnames(puma_crosswalk_1,     "MSA_puma",     "MIGMSA1_puma")
setnames(pumasupr_crosswalk_1, "MSA_pumasupr", "MIGMSA1_pumasupr")
setnames(county_crosswalk_1,   "MSA_county",   "MIGMSA1_county")

all_data <- merge(all_data, puma_crosswalk_1,
                  by.x = c("MIGPLAC1", "MIGPUMA1"),
                  by.y = c("MIGPLAC",  "MIGPUMA"), all.x = TRUE)

all_data <- merge(all_data, pumasupr_crosswalk_1,
                  by.x = c("MIGPLAC1", "MIGPUMA1"),
                  by.y = c("MIGPLAC",  "MIGPUMA"), all.x = TRUE)

all_data <- merge(all_data, county_crosswalk_1,
                  by = c("STATEFIP", "COUNTYFIP"), all.x = TRUE)

# --- MIGMSA5: 5-year migration origin ---
puma_crosswalk_5     <- copy(puma_crosswalk)
pumasupr_crosswalk_5 <- copy(pumasupr_crosswalk)
county_crosswalk_5   <- copy(county_crosswalk)

setnames(puma_crosswalk_5,     "MSA_puma",     "MIGMSA5_puma")
setnames(pumasupr_crosswalk_5, "MSA_pumasupr", "MIGMSA5_pumasupr")
setnames(county_crosswalk_5,   "MSA_county",   "MIGMSA5_county")

all_data <- merge(all_data, puma_crosswalk_5,
                  by.x = c("MIGPLAC5", "MIGPUMA5"),
                  by.y = c("MIGPLAC",  "MIGPUMA"), all.x = TRUE)

all_data <- merge(all_data, pumasupr_crosswalk_5,
                  by.x = c("MIGPLAC5", "MIGPUMA5"),
                  by.y = c("MIGPLAC",  "MIGPUMA"), all.x = TRUE)

all_data <- merge(all_data, county_crosswalk_5,
                  by = c("STATEFIP", "COUNTYFIP"), all.x = TRUE)

# Coalesce to single migration MSA per direction
all_data[, MIGMSA1 := fcoalesce(MIGMSA1_puma, MIGMSA1_pumasupr, MIGMSA1_county)]
all_data[, MIGMSA5 := fcoalesce(MIGMSA5_puma, MIGMSA5_pumasupr, MIGMSA5_county)]

# Drop intermediate crosswalk columns
all_data[, c("MIGMSA1_puma", "MIGMSA1_pumasupr", "MIGMSA1_county",
             "MIGMSA5_puma", "MIGMSA5_pumasupr", "MIGMSA5_county") := NULL]

rm(puma_crosswalk,   pumasupr_crosswalk,   county_crosswalk,
   puma_crosswalk_1, pumasupr_crosswalk_1, county_crosswalk_1,
   puma_crosswalk_5, pumasupr_crosswalk_5, county_crosswalk_5)
gc()

# -----------------------------------------------------------------------------

save(all_data,
     file = file.path(local_checkpoint_dir, "after_sec02.RData"))
cat("Checkpoint saved:", file.path(local_checkpoint_dir, "after_sec02.RData"), "\n")
cat("02_msa_crosswalk.R complete.\n")
