# =============================================================================
# 03_variable_construction.R  --  Section 3: Variable construction
# Part of the modular Immigration_Card_Data_Construction.R pipeline.
# =============================================================================

source("scripts/data_construction/00_setup.R")
load(file.path(local_checkpoint_dir, "after_sec02.RData"))
cat("Loaded checkpoint: after_sec02.RData\n")

# Section 3: Variable construction
# -----------------------------------------------------------------------------

# --- 3.1 Geographic category indicators ---

northeast_states <- c(09, 23, 25, 33, 34, 36, 42, 44, 50)

# Top 30 MSAs by weighted population per year (time-varying flag)
top30_by_year <- all_data[MSA != 0,
                           .(pop = sum(PERWT, na.rm = TRUE)),
                           by = .(YEAR, MSA)]
top30_by_year <- top30_by_year[order(-pop), head(.SD, 30), by = YEAR]
top30_by_year[, IsLargest30 := TRUE]

all_data[, IsLargest30 := FALSE]
all_data[top30_by_year, on = .(YEAR, MSA), IsLargest30 := TRUE]

all_data[, MSA_category := fifelse(
  MSA == 0,        "Rural",
  fifelse(IsLargest30,           "IsLargest30",
  fifelse(STATEFIP == 06,        "IsCalifornia",
  fifelse(STATEFIP == 48,        "IsTexas",
  fifelse(STATEFIP == 12,        "IsFlorida",
  fifelse(STATEFIP %in% northeast_states, "IsNortheast",
          "Other_MSAs"))))))]

# --- 3.2 Origin country indicators ---

all_data[, IsMexico       := BPLD == 20000]
all_data[, IsAsia         := BPLD >= 50000 & BPLD <= 52140]
all_data[, IsCanadaAusEur := (BPLD >= 15000 & BPLD <= 15083) |
                              (BPLD >= 40000 & BPLD <= 49900) |
                              (BPLD >= 70010 & BPLD <= 70014)]

# --- 3.3 Occupation groups (OCC2010) ---
#
# FIXED 2026-09-12. These six substantive categories (Low-Skill Labor &
# Service, Operatives & Craft Workers, Clerical, Sales, Managers,
# Professional & Technical) are named to match Card (2001, JOLE), which
# groups workers into "I. laborers, farm workers, and low-skilled service
# workers; II. operatives and craft workers; III. clerical workers; IV.
# sales workers; V. managers; and VI. professional and technical workers"
# (p.32). Card assigns these via a multinomial logit predicting each
# person's probability of working in each group from demographics
# (including education), applied to the whole population including
# non-workers -- Occup_Group here is instead a direct classification of
# each worker's own observed OCC2010 code, used as a control alongside a
# separate YearSch term in 12_migration_iv.R's Card-following
# specification. Education is already its own control there, so a
# logit-based assignment, which uses education as a predictor, would be
# circular for that purpose, so the classification here is direct.
#
# This review corrected four issues found against the full OCC2010
# codebook (see OCC2010_crosswalk_review.csv), three of them outright
# miscodings:
#   - 9999 (NIU -- not asked the occupation question; mostly children and
#     non-workers) was claimed by the Low-Skill range before the Not in
#     Labor Force branch (AGE<16 | WORKEDYR==0) could ever fire. This was
#     responsible for a ~52% "Low-Skill" share on the full sample.
#   - 9030/9040 (Aircraft Pilots and Flight Engineers; Air Traffic
#     Controllers) were hardcoded into Operatives & Craft Workers; moved
#     to Professional & Technical.
#   - 6200 (First-Line Supervisors of Construction Trades and Extraction
#     Workers) was listed under Managers but unreachable, since Operatives
#     & Craft's 6200:6765 range was checked first in fcase() order and
#     always claimed it. Range narrowed to 6210:6765 so 6200 now reaches
#     the Managers branch it was already listed in.
#   - Farming/Fishing/Forestry (6005, 6010, 6040, 6050, 6100, 6120, 6130)
#     was split three ways across Managers/Professional/Operatives.
#     Card's Group I explicitly names "farm workers" alongside laborers
#     and low-skilled service workers, so the whole set moves to Low-Skill
#     to match his definition.
#
# SECOND PASS (2026-09-12), after checking group_counts proportions on the
# first pass above: renormalized across the six substantive groups
# (excluding Not in Labor Force and Military), Professional & Technical
# came in at ~25%, noticeably above a defensible 15-20% band, while the
# other five were all in or near that band. The cause was Healthcare
# Support (3600-3655) and all of Protective Service (3700-3955) still
# being swept in via the 3000:3955 range meant for Healthcare
# Practitioners (3000-3540).
#
# Moving those two blocks to Low-Skill Labor & Service isn't just a fix
# for the proportions -- it matches the standard 1990 Census
# major-occupation-group taxonomy Card's six groups are evidently built
# from, which places Protective Service alongside food service, cleaning,
# and personal care under one broad "Service occupations" category,
# distinct from "Managerial and Professional Specialty." That is the same
# grouping Card's Group I name ("low-skilled service workers") describes.
#
# Business Operations Specialists/Financial Specialists remain folded into
# Managers via 500:950: Managers alone already landed at ~14.6%, right at
# the edge of the 15-20% band, so pulling specialists out would only move
# it further away, and Card's paper gives no basis for splitting this one
# further.

all_data[, Occup_Group := fcase(
  OCC2010 %in% c(4000:4150, 4200:4250, 4300:4650,
                 6005, 6010, 6040, 6050, 6100, 6120, 6130,
                 3600, 3610, 3620, 3630, 3640, 3645, 3646,
                 3647, 3648, 3649, 3655,
                 3700, 3710, 3720, 3730, 3740, 3750,
                 3800, 3820, 3840, 3850,
                 3900, 3910, 3930, 3940, 3945, 3955,
                 9000, 9050, 9110, 9120, 9130, 9140, 9150,
                 9200, 9240, 9260, 9300, 9310, 9350, 9360,
                 9410, 9415, 9420, 9510, 9520, 9560, 9600,
                 9610, 9620, 9630, 9640, 9650, 9720, 9750),
  "Low-Skill Labor & Service",

  OCC2010 %in% c(6210:6765, 6800:6940, 7000:7630, 7700:8965),
  "Operatives & Craft Workers",

  OCC2010 %in% 5000:5940,
  "Clerical Workers",

  OCC2010 %in% 4700:4965,
  "Sales Workers",

  OCC2010 %in% c(10:430, 500:950, 6200),
  "Managers",

  OCC2010 %in% c(1000:1980, 2000:2150, 2160, 2200:2920,
                 3000:3540, 9030, 9040),
  "Professional & Technical Workers",

  OCC2010 %in% 9800:9830,
  "Military Occupations",

  OCC2010 == 9920 | OCC2010 == 9999 | AGE < 16 | WORKEDYR == 0,
  "Not in Labor Force",

  default = "Other/Unclassified"
)]

# Validate occupation group distribution
group_counts <- all_data[, .N, by = Occup_Group][order(-N)]
group_counts[, Proportion := round(N / sum(N), 4)]
print(group_counts)

# --- 3.4 Education groups (EDUCD) ---

all_data[, Educ_Group := fcase(
  EDUCD %in% c(000:030), "non_HS",
  EDUCD %in% c(040:050), "non_HS",
  EDUCD %in% c(060:064), "HighSch",
  EDUCD %in% c(065:100, 101), "College",
  EDUCD %in% c(110:116), "Grad",
  EDUCD == 999, NA_character_,
  default = NA_character_
)]

all_data[, Educ_Group := factor(Educ_Group,
                                 levels = c("non_HS", "HighSch",
                                            "College", "Grad"))]

# Validate education distribution (working-age only)
Educ_table <- all_data[AGE >= 16, .N, by = Educ_Group][order(-N)]
Educ_table[, Proportion := round(N / sum(N), 4)]
print(Educ_table)

# --- 3.5 Immigration and region of birth ---

all_data[, immigr := ifelse(BPL >= 120, 1, 0)]

all_data[, RegionBPL := fcase(
  BPL == 500,                             "China",
  BPL == 521,                             "India",
  BPL == 200,                             "Mexico",
  BPL %in% 453:454,                       "Germany",
  BPL %in% 520:547,                       "Middle East",
  BPL %in% 001:120,                       "USA",
  BPL %in% 500:599 & !(BPL %in% c(500, 521)), "Asia",
  BPL %in% 400:499 & !(BPL %in% 453:454), "Europe",
  BPL %in% c(150, 155, 160, 199),         "North America",
  BPL %in% 210:249,                       "Central America",
  BPL %in% c(250:299, 110, 115),          "Caribbean",
  BPL == 300,                             "South America",
  BPL == 600,                             "Africa",
  BPL %in% c(700, 710),                   "Oceania",
  BPL == 900,                             "Abroad, not specified",
  BPL == 800,                             "At Sea",
  default = "Other"
)]

# --- 3.6 Years in the U.S. and arrival before age 6 ---

all_data[, YearsInUS  := ifelse(YRIMMIG > 0, YEAR - YRIMMIG, NA_real_)]
all_data[, YearsInUS2 := YearsInUS^2]
all_data[, ArrBefore6 := ifelse(immigr == 1 & YearsInUS >= AGE - 6, 1, 0)]

# --- 3.7 Hourly wage ---
# hrwage = annual earnings / (usual hours per week * weeks worked)
# Weeks worked derived from categorical WKSWORK2

all_data[, WEAKWK := fcase(
  WKSWORK2 == 0, 0,
  WKSWORK2 == 1, 13,
  WKSWORK2 == 2, 26,
  WKSWORK2 == 3, 39,
  WKSWORK2 == 4, 47,
  WKSWORK2 == 5, 49,
  default = NA_real_
)]

all_data[, hrwage := INCEARN / (UHRSWORK * WEAKWK)]
all_data[UHRSWORK == 0 | WEAKWK == 0 |
         is.na(UHRSWORK) | is.na(WEAKWK), hrwage := NA_real_]

# --- 3.8 Race indicators ---

all_data[, Black     := as.integer(RACED %in% 200:210)]
all_data[, Aboriginal := as.integer(RACED %in% c(300:328, 350,
                                                   352:353, 355:358,
                                                   361:363, 370:399))]
all_data[, Asian     := as.integer(RACED %in% c(400:620, 640:643,
                                                  651, 656:679))]

# --- 3.9 Migration / outflow flag ---
# MOVED == 1 if the individual lived in a different PUMA or state
# one year ago (MIGPLAC1/MIGPUMA1) or five years ago (MIGPLAC5/MIGPUMA5)
#
# NOTE: STATEFIP/PUMA and MIGPLAC1/MIGPUMA1/MIGPLAC5/MIGPUMA5 are all
# IPUMS-native columns, each carrying its own codebook value labels
# (haven_labelled). Comparing two labelled columns directly (STATEFIP ==
# MIGPLAC1, etc.) triggers harmless-but-noisy "..1 and ..2 have conflicting
# value labels" warnings from haven's label-reconciliation logic when the
# two columns' codebooks disagree on what a shared numeric code means (the
# underlying numeric comparison is correct either way -- only the label
# metadata is ambiguous). as.integer() on both sides strips the labels and
# removes the warning, same pattern as the PUMA/PUMASUPR and MET2013 fixes
# in 02_msa_crosswalk.R.
STATEFIP_int <- as.integer(all_data$STATEFIP)
PUMA_int     <- as.integer(all_data$PUMA)
MIGPLAC1_int <- as.integer(all_data$MIGPLAC1)
MIGPUMA1_int <- as.integer(all_data$MIGPUMA1)
MIGPLAC5_int <- as.integer(all_data$MIGPLAC5)
MIGPUMA5_int <- as.integer(all_data$MIGPUMA5)

all_data[, MOVED := fcase(
  !is.na(MIGPLAC1_int) & !is.na(MIGPUMA1_int) & MIGPLAC1_int != 0L & MIGPUMA1_int > 3L &
    ((STATEFIP_int == MIGPLAC1_int & PUMA_int != MIGPUMA1_int) |
     (STATEFIP_int != MIGPLAC1_int & PUMA_int == MIGPUMA1_int) |
     (STATEFIP_int != MIGPLAC1_int & PUMA_int != MIGPUMA1_int)), 1L,

  !is.na(MIGPLAC5_int) & !is.na(MIGPUMA5_int) & MIGPLAC5_int != 0L & MIGPUMA5_int != 0L &
    ((STATEFIP_int == MIGPLAC5_int & PUMA_int != MIGPUMA5_int) |
     (STATEFIP_int != MIGPLAC5_int & PUMA_int == MIGPUMA5_int) |
     (STATEFIP_int != MIGPLAC5_int & PUMA_int != MIGPUMA5_int)), 1L,

  default = 0L
)]

rm(STATEFIP_int, PUMA_int, MIGPLAC1_int, MIGPUMA1_int, MIGPLAC5_int, MIGPUMA5_int)

# --- 3.10 Parental education ---

all_data[, Mom_Educ := fcase(
  is.na(EDUC_MOM),                      "Missing",
  EDUC_MOM == 6,                         "HighSch",
  EDUC_MOM %in% c(7, 8, 9, 10),         "College",
  EDUC_MOM == 11,                        "Grad",
  default = "non_HS"
)]

all_data[, Dad_Educ := fcase(
  is.na(EDUC_POP),                      "Missing",
  EDUC_POP == 6,                         "HighSch",
  EDUC_POP %in% c(7, 8, 9, 10),         "College",
  EDUC_POP == 11,                        "Grad",
  default = "non_HS"
)]

all_data[, Mom_Educ := factor(Mom_Educ,
                               levels = c("non_HS", "Missing",
                                          "HighSch", "College", "Grad"))]
all_data[, Dad_Educ := factor(Dad_Educ,
                               levels = c("non_HS", "Missing",
                                          "HighSch", "College", "Grad"))]

# --- 3.11 Binary education outcome ---

all_data[, Educ_Binary := as.integer(Educ_Group %in% c("College", "Grad"))]

# --- 3.12 Potential experience ---
# Standard Mincer definition: Age - Years of schooling - 6

all_data[, YearSch := case_when(
  EDUCD %in% 000:010 ~ 0,
  EDUCD %in% 011:030 ~ 4,
  EDUCD == 040       ~ 8,
  EDUCD == 050       ~ 9,
  EDUCD == 060       ~ 10,
  EDUCD == 061       ~ 11,
  EDUCD == 062       ~ 11,
  EDUCD == 063       ~ 12,
  EDUCD == 064       ~ 12,
  EDUCD == 065       ~ 13,
  EDUCD == 070       ~ 13,
  EDUCD == 071       ~ 13,
  EDUCD == 080       ~ 14,
  EDUCD == 081       ~ 14,
  EDUCD == 082       ~ 14,
  EDUCD == 083       ~ 14,
  EDUCD == 090       ~ 14,
  EDUCD == 100       ~ 16,
  EDUCD == 101       ~ 16,
  EDUCD == 110       ~ 17,
  EDUCD == 111       ~ 17,
  EDUCD == 112       ~ 17,
  EDUCD == 113       ~ 18,
  EDUCD == 114       ~ 18,
  EDUCD == 115       ~ 17,
  EDUCD == 116       ~ 20,
  TRUE               ~ NA_real_
)]

all_data[, potential_experience := AGE - YearSch - 6]
all_data[potential_experience < 0, potential_experience := 0]

# --- 3.13 Child age group ---
#
# FIXED 2026-09-12: this was previously built from the RESPONDENT's own AGE
# (AGE < 5 -> "Under5", etc.), which made "NoChild" mean "respondent is 18
# or older" rather than "has no children" -- i.e. it was an age-of-
# respondent proxy, not a household-composition control, and its name was
# misleading throughout the pipeline. It also happened to partially absorb
# the under-18 contamination diagnosed in the estimation-sample age-floor
# work (see 08_load_analysis_subsets.R), and it collapses to a single level
# under an AGE >= 30 estimation sample -- the immediate reason this needed
# fixing.
#
# Now built from YNGCH (age of the youngest own child, including step/
# adopted, co-resident with the respondent; 99 = no own children present).
# YNGCH can run up to 98, but only ages through 17 get their own bucket --
# a co-resident child 18 or older is not a dependent in the sense this
# control is meant to capture, so those rows are folded into "NoChild"
# alongside respondents with no co-resident children at all.
all_data[, ChildAgeGroup := fcase(
  YNGCH == 99 | YNGCH >= 18, "NoChild",
  YNGCH %in% 0:4,            "Under5",
  YNGCH %in% 5:12,           "5to12",
  YNGCH %in% 13:17,          "13to17"
)]

all_data[, Marital := as.integer(MARST %in% 1:2)]

# Save all_data with all constructed variables
save(all_data, file = "all_data.RData")


save(all_data,
     file = file.path(local_checkpoint_dir, "after_sec03.RData"))
cat("Checkpoint saved:", file.path(local_checkpoint_dir, "after_sec03.RData"), "\n")
cat("03_variable_construction.R complete.\n")
