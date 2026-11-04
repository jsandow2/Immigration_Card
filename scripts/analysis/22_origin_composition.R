# =============================================================================
# 22_origin_composition.R  --  Descriptive statistics: nativity, and the
# Renumbered 2026-10-21 from 21_origin_composition.R (archive-and-renumber pass).
#                              composition of immigration by origin
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# Supersedes AND MERGES two scripts, both archived 2026-10-02:
#   21_descriptives_mexico_other.R   nativity descriptives + a Mexico/Other split
#   22_demographic_summary.R         demographic profile of the top 5 origins
#
# WHY THE MEXICO / OTHER SPLIT IS GONE
# The old 21 existed to answer one question: why did dropping Mexico RAISE the
# estimate? Its header said the Mexican-vs-other composition "motivates the
# Rotemberg sensitivity finding in Section 19 -- Mexico's large IV weight
# reflects its distinctive low-education, low-skill composition."
#
# That question has since been answered, and not by composition. The anomalous
# sign was a construction artefact: 07_exposure_iv.R merged actual inflows on
# SURVEY year while the instrument was built over the ages-15-25 window. With
# that fixed, Mexico carries alpha = 0.189 and beta_k = +1.264e-05 -- the
# largest POSITIVE weight in the decomposition. There is no anomaly left for a
# composition table to explain.
#
# What survives is the underlying descriptive interest: origin countries differ
# sharply in the education their migrants bring, and that is worth showing in
# its own right. A Mexico/Other dummy is the wrong instrument for it, because
# it collapses every other origin into a single residual and hides exactly the
# variation that makes the fact interesting.
#
# STRUCTURE
#   A  Nativity -- immigrant vs US-born, by gender, by education, over time
#   B  Parental education missingness (for the Data section's limitations)
#   C  Education composition across the top TOP_N origins   [NEW]
#      Wide on countries, narrow on characteristics.
#   D  Full demographic profile, top 5 origins vs US-born   [was Section 22]
#      Narrow on countries, wide on characteristics.
#
# C and D are deliberately complementary rather than redundant: C shows how far
# apart origins are on the one dimension this dissertation is about, across
# enough countries for the spread to be visible; D gives the conventional
# five-row table with income and gender attached.
#
# AGE >= 25 THROUGHOUT. Educ_Group is assigned in the construction script from
# EDUCD alone with no age floor, so anyone still in school is bucketed as
# non_HS -- the 1990 US-born male non-HS share came out at 42% before this
# restriction against roughly 25% nationally for ages 25+, and roughly 20%
# after. AGE >= 25 is the standard convention for attainment statistics. This
# is a descriptive table and not the estimation sample; 23 covers the latter,
# which carries no age restriction.
#
# SELF-SUFFICIENT. Loads analysis_data.RData directly rather than chaining
# through after_sec20 / after_sec21, which will not regenerate on a clean
# rebuild. Still writes after_sec22.RData because 25_dissertation_table_helper.R
# reads it.
#
# Writes outputs/tables/{parent_educ_missingness, educ_by_origin,
#                        demographic_summary_table, bpld_codes}.csv
# =============================================================================

# Clears the session but KEEPS bpld_lookup. It is built from the IPUMS DDI in
# 08_load_analysis_subsets.R and is not carried in analysis_data.RData, so a
# blanket rm(list = ls()) would delete the one object this script cannot
# rebuild -- and would do it in the very session that just created it, making
# the guard below unsatisfiable no matter what order things were run in.
rm(list = setdiff(ls(), "bpld_lookup")); gc()
source("scripts/analysis/00_setup_analysis.R")

library(data.table)

TOP_N   <- 20L   # origins in Panel C
TOP_N_D <- 5L    # origins in Panel D (the conventional table)

cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
reduced_data_iv <- as.data.table(.tmp_env$reduced_data_iv)
all_natives_175 <- as.data.table(.tmp_env$all_natives_175)

# MSA_Name is not carried on the microdata; it lives in top175_list_named,
# the per-year ranking built in 05. Pulled here, while the environment still
# exists, for the metropolitan concentration output in Panel E. Left NULL if
# absent so an older analysis_data.RData reports codes instead of failing.
msa_name_lookup <- if (!is.null(.tmp_env$top175_list_named)) {
  unique(rbindlist(lapply(.tmp_env$top175_list_named,
                          function(x) as.data.table(x)[, .(MSA, MSA_Name)])))
} else NULL

rm(.tmp_env); gc()

# bpld_lookup is built from the IPUMS DDI in 08_load_analysis_subsets.R and is
# not carried in analysis_data.RData, so it has to already be in the session.
if (!exists("bpld_lookup"))
  stop("bpld_lookup not found. Run 08_load_analysis_subsets.R first -- it is\n",
       "built from the IPUMS DDI and is not in analysis_data.RData.")
bl <- as.data.table(bpld_lookup)

desc_data <- reduced_data_iv[AGE >= 25]
cat(sprintf("Descriptive sample restricted to AGE >= 25: %s of %s obs retained (%.1f%%)\n",
            format(nrow(desc_data), big.mark = ","),
            format(nrow(reduced_data_iv), big.mark = ","),
            100 * nrow(desc_data) / nrow(reduced_data_iv)))

# YEAR arrives as a factor on some paths and numeric on others. Resolving it
# once here avoids the silent-failure case where a factor comparison like
# YEAR >= 2000 compares level ORDER rather than year.
desc_data[, YEAR_n := as.integer(as.character(YEAR))]

# IPUMS columns arrive as haven_labelled. The class is harmless in comparisons
# and arithmetic, which is why it goes unnoticed, but it breaks every function
# that checks classes strictly -- rbind() refuses a column whose class differs
# across items, fifelse() refuses mismatched `yes`/`no`, and scales::rescale()
# has no method for it at all. Stripping the whole table once here rather than
# discovering them one function at a time, which is how the first two surfaced.
.hv <- names(desc_data)[vapply(desc_data, inherits, logical(1), "haven_labelled")]
if (length(.hv)) {
  cat(sprintf("Stripping haven_labelled from %d column(s): %s\n",
              length(.hv), paste(.hv, collapse = ", ")))
  desc_data[, (.hv) := lapply(.SD, function(x) {
    v <- unclass(x); attributes(v) <- NULL; v
  }), .SDcols = .hv]
}
rm(.hv)

rm(reduced_data_iv); gc()

# =============================================================================
# A. NATIVITY
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("A. NATIVITY (AGE >= 25)\n")
cat(strrep("=", 70), "\n", sep = "")

# THE NATIVITY SPLIT IS A DESIGN ARTEFACT, NOT AN ESTIMATE.
# 04_sample_construction.R keeps foreign-born rows IN FULL and draws 25% of
# native-born rows within each year. So the sample runs about 41% immigrant
# against a national foreign-born share near 14%, and that gap is the sampling
# design, not the population.
#
# PERWT DOES NOT FIX THIS -- it makes it worse. 04 deliberately leaves native
# PERWT unrescaled after the 25% draw (its own comment says so), so weighting
# understates natives by a factor of four rather than correcting anything. The
# 1990/2000 thinning rescale is a separate correction and is already applied.
#
# The counts are therefore reported as counts, with the implied correction
# shown beside them rather than silently applied. Nothing in this chapter needs
# the nativity share as an estimate; anyone who wants it should cite the Census.
NATIVE_DRAW <- 0.25

cat("\n--- Immigrant vs U.S.-born (OBSERVATION COUNTS, not population shares) ---\n")
nativity <- desc_data[, .(n = .N), by = immigr]
nativity[, group := ifelse(immigr == 1, "Immigrant", "U.S.-born")]
nativity[, sample_share := n / sum(n)]
# Undo the 25% native draw to show what the split would look like at the
# design's own implied population. This is arithmetic on a known sampling
# fraction, not a reweighting of the data, and it is printed for orientation
# only -- no downstream table uses it.
nativity[, n_design_adj := fifelse(immigr == 1, as.numeric(n),
                                   as.numeric(n) / NATIVE_DRAW)]
nativity[, implied_share := n_design_adj / sum(n_design_adj)]
print(nativity[order(-immigr), .(group, n, sample_share, implied_share)],
      row.names = FALSE)
cat(sprintf("  sample is %.1f%% immigrant by row count; %.1f%% after undoing the\n",
            100 * nativity[immigr == 1, sample_share],
            100 * nativity[immigr == 1, implied_share]))
cat(sprintf("  %.0f%% native draw. Neither is a population estimate -- the first is\n",
            100 * NATIVE_DRAW))
cat("  the design, the second assumes the draw is the only departure from it.\n")
cat("  Every panel below is a WITHIN-GROUP composition and is unaffected by\n")
cat("  this, because the native draw is random within year.\n")

cat("\n--- Gender by nativity ---\n")
gender_nativity <- desc_data[, .(n = .N), by = .(immigr, SEX)]
gender_nativity[, share := n / sum(n), by = immigr]
gender_nativity[, group  := ifelse(immigr == 1, "Immigrant", "U.S.-born")]
gender_nativity[, gender := ifelse(SEX == 1, "Male", "Female")]
print(gender_nativity[order(immigr, SEX), .(group, gender, n, share)],
      row.names = FALSE)

cat("\n--- Gender x education by nativity ---\n")
gender_educ_nativity <- desc_data[
  !is.na(Educ_Group)
][, .(n = .N), by = .(immigr, SEX, Educ_Group)
][, share := n / sum(n), by = .(immigr, SEX)
][, group  := ifelse(immigr == 1, "Immigrant", "U.S.-born")
][, gender := ifelse(SEX == 1, "Male", "Female")]
print(gender_educ_nativity[order(immigr, SEX, Educ_Group),
                           .(group, gender, Educ_Group, n, share)],
      row.names = FALSE)

# Immigrants are grouped by year of ARRIVAL (YRIMMIG) because education is
# largely fixed at arrival; the US-born by survey year. Same cohort convention
# as the Educ_Prop figure. Restricted to actual sample years because YRIMMIG is
# a continuous calendar year and is not tied to survey design.
cat("\n--- Gender x education by nativity, over time ---\n")
sample_years <- c(1990, 2000, 2006:2023)

immig_gender_educ_yr <- desc_data[
  immigr == 1 & YRIMMIG %in% sample_years & !is.na(Educ_Group)
][, .(n = .N), by = .(Year = as.integer(YRIMMIG), SEX, Educ_Group)
][, share := n / sum(n), by = .(Year, SEX)
][, `:=`(group = "Immigrant", gender = ifelse(SEX == 1, "Male", "Female"))]

native_gender_educ_yr <- desc_data[
  immigr == 0 & YEAR_n %in% sample_years & !is.na(Educ_Group)
][, .(n = .N), by = .(Year = YEAR_n, SEX, Educ_Group)
][, share := n / sum(n), by = .(Year, SEX)
][, `:=`(group = "U.S.-born", gender = ifelse(SEX == 1, "Male", "Female"))]

gender_educ_nativity_yr <- rbind(
  immig_gender_educ_yr[, .(group, gender, Year, Educ_Group, n, share)],
  native_gender_educ_yr[, .(group, gender, Year, Educ_Group, n, share)]
)
print(gender_educ_nativity_yr[order(group, gender, Year, Educ_Group)],
      row.names = FALSE)

# =============================================================================
# B. PARENTAL EDUCATION MISSINGNESS
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("B. PARENTAL EDUCATION MISSINGNESS (all_natives_175)\n")
cat(strrep("=", 70), "\n\n", sep = "")
parent_educ_missing <- data.frame(
  Variable    = c("Mom_Educ", "Dad_Educ"),
  Pct_Missing = c(100 * mean(all_natives_175$Mom_Educ == "Missing"),
                  100 * mean(all_natives_175$Dad_Educ == "Missing"))
)
print(parent_educ_missing, row.names = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(parent_educ_missing, "outputs/tables/parent_educ_missingness.csv",
          row.names = FALSE)

# =============================================================================
# C. EDUCATION COMPOSITION BY ORIGIN
# =============================================================================
# Origins ranked by pooled average population-weighted count: PERWT summed
# within each survey year, then averaged ACROSS years so each survey year
# counts equally. Without that averaging the denser 2006-2023 ACS coverage
# swamps the single 2000 5% sample and the ranking becomes a statement about
# survey design. Same ranking rule as Panel D, so the two tables nest.
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("C. EDUCATION COMPOSITION, TOP %d ORIGINS\n", TOP_N))
cat(strrep("=", 70), "\n", sep = "")

origin_rank <- desc_data[
  immigr == 1 & YEAR_n >= 2000 & YEAR_n <= 2023,
  .(pop_weighted = sum(PERWT, na.rm = TRUE)), by = .(BPLD, YEAR_n)
][, .(avg_pop_weighted = mean(pop_weighted)), by = BPLD
][order(-avg_pop_weighted)]
origin_rank <- merge(origin_rank, bl, by = "BPLD", all.x = TRUE)
setorder(origin_rank, -avg_pop_weighted)

top_origins <- head(origin_rank, TOP_N)
total_immig_n <- desc_data[immigr == 1, .N]

educ_by_origin <- desc_data[
  immigr == 1 & BPLD %in% top_origins$BPLD & !is.na(Educ_Group),
  .(N            = .N,
    Pct_female   = 100 * mean(SEX == 2, na.rm = TRUE),
    Pct_nonHS    = 100 * mean(Educ_Group == "non_HS",  na.rm = TRUE),
    Pct_HS       = 100 * mean(Educ_Group == "HighSch", na.rm = TRUE),
    Pct_College  = 100 * mean(Educ_Group == "College", na.rm = TRUE),
    Pct_Grad     = 100 * mean(Educ_Group == "Grad",    na.rm = TRUE)),
  by = BPLD]
educ_by_origin[, Pct_CollegePlus := Pct_College + Pct_Grad]
educ_by_origin[, Pct_of_immig    := 100 * N / total_immig_n]

# PERWT-weighted college-or-more, alongside the unweighted figures above.
# The unweighted basis is inherited from the old Section 22's make_summary;
# the weighted one is what a population composition actually calls for, and is
# what figure (d) in 27_education_figures.R draws. Both are reported so the
# table and the figure can be reconciled rather than silently disagreeing.
# If the two diverge materially for an origin, that origin's sampling weights
# are unusual and it is worth knowing before either number goes in the text.
wtd <- desc_data[immigr == 1 & BPLD %in% top_origins$BPLD & !is.na(Educ_Group),
                 .(w = sum(PERWT, na.rm = TRUE)), by = .(BPLD, Educ_Group)]
wtd[, share := w / sum(w), by = BPLD]
# Summed directly rather than via dcast: a dcast column is only created if some
# origin actually has that education level, so an empty category would silently
# remove the column and the next line would fail on a missing object.
wtd <- wtd[, .(Pct_CollegePlus_wtd = 100 * sum(share[Educ_Group %in%
                                                     c("College", "Grad")]),
               Pct_nonHS_wtd       = 100 * sum(share[Educ_Group == "non_HS"])),
           by = BPLD]
educ_by_origin <- merge(educ_by_origin, wtd, by = "BPLD", all.x = TRUE)
educ_by_origin <- merge(educ_by_origin, top_origins[, .(BPLD, label,
                                                        avg_pop_weighted)],
                        by = "BPLD")
setorder(educ_by_origin, -avg_pop_weighted)
setcolorder(educ_by_origin, c("label", "BPLD", "N", "Pct_of_immig",
                              "Pct_female", "Pct_nonHS", "Pct_HS",
                              "Pct_College", "Pct_Grad", "Pct_CollegePlus"))

# Reference rows. The US-born line is what makes Panel C readable: without it
# the origin figures are numbers without a scale.
ref <- desc_data[!is.na(Educ_Group), .(
  N               = .N,
  Pct_female      = 100 * mean(SEX == 2, na.rm = TRUE),
  Pct_nonHS       = 100 * mean(Educ_Group == "non_HS",  na.rm = TRUE),
  Pct_HS          = 100 * mean(Educ_Group == "HighSch", na.rm = TRUE),
  Pct_College     = 100 * mean(Educ_Group == "College", na.rm = TRUE),
  Pct_Grad        = 100 * mean(Educ_Group == "Grad",    na.rm = TRUE)),
  by = .(label = ifelse(immigr == 1, "ALL IMMIGRANTS", "U.S.-BORN"))]
ref[, `:=`(Pct_CollegePlus = Pct_College + Pct_Grad, BPLD = NA_integer_,
           Pct_of_immig = NA_real_, avg_pop_weighted = NA_real_)]

cat(sprintf("\n  Top %d origins account for %.1f%% of immigrant observations.\n",
            TOP_N, educ_by_origin[, sum(Pct_of_immig)]))
cat("\n--- By size (share of immigrant observations) ---\n")
print(rbind(educ_by_origin, ref, fill = TRUE)[
  , .(label, N, Pct_of_immig, Pct_nonHS, Pct_HS, Pct_College, Pct_Grad,
      Pct_CollegePlus)], row.names = FALSE, digits = 3)

# The same table sorted on the outcome dimension. Size order buries the
# finding; this ordering is the figure-in-a-table and is what the text should
# quote from.
cat("\n--- Sorted by college-or-more, which is where the spread is ---\n")
cat("    (_wtd columns are PERWT-weighted and are what figure (d) in 15 draws)\n")
print(educ_by_origin[order(-Pct_CollegePlus),
                     .(label, N, Pct_nonHS, Pct_CollegePlus,
                       Pct_nonHS_wtd, Pct_CollegePlus_wtd)],
      row.names = FALSE, digits = 3)
cat(sprintf("\n  Max |weighted - unweighted| on college-or-more: %.2f pp (%s)\n",
            educ_by_origin[, max(abs(Pct_CollegePlus_wtd - Pct_CollegePlus))],
            educ_by_origin[which.max(abs(Pct_CollegePlus_wtd -
                                         Pct_CollegePlus)), label]))

us_cp <- ref[label == "U.S.-BORN", Pct_CollegePlus]
cat(sprintf("\n  Range across these %d origins: %.1f%% to %.1f%% college-or-more",
            nrow(educ_by_origin), educ_by_origin[, min(Pct_CollegePlus)],
            educ_by_origin[, max(Pct_CollegePlus)]))
cat(sprintf("\n  (a factor of %.1f), against %.1f%% for the U.S.-born.\n",
            educ_by_origin[, max(Pct_CollegePlus) / min(Pct_CollegePlus)],
            us_cp))
cat(sprintf("  %d of %d origins sit ABOVE the U.S.-born rate.\n",
            educ_by_origin[Pct_CollegePlus > us_cp, .N], nrow(educ_by_origin)))
cat("\n  This spread is the salient point for the text. 'Immigrants' is not a\n")
cat("  composition -- it is an average over origins that differ by more on\n")
cat("  this margin than the immigrant/native gap itself.\n")

write.csv(educ_by_origin, "outputs/tables/educ_by_origin.csv", row.names = FALSE)

# =============================================================================
# D. DEMOGRAPHIC PROFILE, TOP 5 ORIGINS vs U.S.-BORN
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("D. DEMOGRAPHIC PROFILE, TOP %d ORIGINS vs U.S.-BORN\n", TOP_N_D))
cat(strrep("=", 70), "\n", sep = "")

if (!"income_clean" %in% names(desc_data))
  desc_data[, income_clean := fifelse(INCEARN > 0 & INCEARN < 999999,
                                      INCEARN, NA_real_)]

top5_bpld <- head(origin_rank, TOP_N_D)
cat(sprintf("\n--- Top %d origins (pooled average population-weighted count, 2000-2023) ---\n",
            TOP_N_D))
print(top5_bpld[, .(BPLD, label, avg_pop_weighted)], row.names = FALSE)

make_summary <- function(data, group_label, total_immig = NULL) {
  out <- data[, .(
    N             = .N,
    Pct_female    = mean(SEX == 2, na.rm = TRUE) * 100,
    Pct_nonHS     = mean(Educ_Group == "non_HS",  na.rm = TRUE) * 100,
    Pct_HS        = mean(Educ_Group == "HighSch", na.rm = TRUE) * 100,
    Pct_College   = mean(Educ_Group == "College", na.rm = TRUE) * 100,
    Pct_Grad      = mean(Educ_Group == "Grad",    na.rm = TRUE) * 100,
    Mean_income   = mean(income_clean,   na.rm = TRUE),
    Median_income = median(income_clean, na.rm = TRUE))]
  out[, Group        := group_label]
  out[, Pct_of_immig := if (is.null(total_immig)) NA_real_
                        else N / total_immig * 100]
  setcolorder(out, c("Group", "N", "Pct_female", "Pct_nonHS", "Pct_HS",
                     "Pct_College", "Pct_Grad", "Mean_income",
                     "Median_income", "Pct_of_immig"))
  out
}

summary_rows <- vector("list", nrow(top5_bpld) + 3L)
# "Combined" rather than "Full sample": every row here is AGE >= 25, so putting
# the restriction on one label would imply the others lack it, and "full"
# suggests a population the sampling design does not deliver. "Combined" says
# what the row actually is -- immigrants and natives added together, in the
# design's proportions. The age restriction belongs to the whole table and is
# stated in the caption.
names(summary_rows) <- c(top5_bpld$label, "All Immigrants", "U.S.-born",
                         "Combined")
for (i in seq_len(nrow(top5_bpld)))
  summary_rows[[top5_bpld$label[i]]] <- make_summary(
    desc_data[immigr == 1 & BPLD == top5_bpld$BPLD[i]],
    top5_bpld$label[i], total_immig_n)
summary_rows[["All Immigrants"]] <- make_summary(desc_data[immigr == 1],
                                                 "All Immigrants", total_immig_n)
summary_rows[["U.S.-born"]]      <- make_summary(desc_data[immigr == 0],
                                                 "U.S.-born")
summary_rows[["Combined"]]       <- make_summary(desc_data, "Combined")

demographic_summary_table <- rbindlist(summary_rows)
cat("\n--- Demographic summary table ---\n")
print(demographic_summary_table, row.names = FALSE, digits = 4)

# WARNING ON THE "Combined" ROW. Every other row is a within-group
# composition and is clean: the native draw is random within year, so native
# composition is unbiased, and immigrant rows are a complete enumeration.
# "Combined" is the one row that MIXES the two, so it averages them in the
# design's 41/59 proportions rather than the population's ~14/86. It is not a
# statistic about anything.
#
# Left in place rather than removed because this table feeds
# demographic_summary_table -> after_sec22 -> 24 -> the LaTeX, and dropping a
# row silently would change a table already drafted against. DECISION NEEDED:
# either drop the row or caption it as a sample average. Do not report it as
# a population figure.
cat("\n  NOTE: the 'Combined' row mixes natives and immigrants in the\n")
cat("  sampling design's proportions (foreign-born kept in full, 25% native\n")
cat("  draw), so it is not a population average. Every other row is a\n")
cat("  within-group composition and is unaffected.\n")

# =============================================================================
# PANEL E -- Mexico against all other origins
# =============================================================================
# RESTORED 2026-10-19. These breakdowns were produced by
# 21_descriptives_mexico_other.R, which was archived when 21 and 22 were
# merged. The merge dropped the outputs without dropping the claims resting on
# them: ten figures in the Data section traced to CSVs last written
# 2026-07-24 by a script no longer in the pipeline, so no current run could
# reproduce them. They are recomputed here from desc_data, on the AGE >= 25
# restriction every other panel uses.
#
# Mexico is identified through the codebook lookup rather than the IsMexico
# indicator, so this panel resolves origins exactly as the demographic summary
# above does and cannot silently disagree with it.
mex_bpld <- bl[label == "Mexico", BPLD]
stopifnot(length(mex_bpld) == 1L)

mex_base <- copy(desc_data[immigr == 1])
mex_base[, Origin := fifelse(BPLD == mex_bpld, "Mexico", "Other")]

# Shares are within-origin, so neither is affected by the sampling design's
# overrepresentation of the foreign-born -- both columns are immigrants.
comp_by <- function(d, var, dim_label) {
  out <- d[!is.na(get(var)),
           .(n = .N), by = .(Origin, Category = as.character(get(var)))]
  out[, share := 100 * n / sum(n), by = Origin]
  out[, Dimension := dim_label]
  setcolorder(out, c("Dimension", "Origin", "Category", "n", "share"))
  out[order(Dimension, Origin, -share)]
}

# Occupation is reported on two denominators because they answer different
# questions and the chapter's wording decides which it needs. "Occupation
# (all)" is over everyone aged 25+; "Occupation (employed)" drops Not in
# Labor Force. The distinction matters here in a way it did not in July:
# Occup_Group now assigns NILF on OCC2010 9920/9999, AGE < 16 OR
# WORKEDYR == 0, and that last clause puts roughly a quarter of both origin
# groups in NILF where the earlier construction put well under one per cent.
# A sentence saying someone "works in" an occupation needs the employed
# denominator; a sentence describing the population needs the other.
mexico_vs_other <- rbindlist(list(
  comp_by(mex_base, "Educ_Group",  "Education"),
  comp_by(mex_base, "Occup_Group", "Occupation (all)"),
  comp_by(mex_base[Occup_Group != "Not in Labor Force"],
          "Occup_Group", "Occupation (employed)"),
  comp_by(mex_base, "SEX",         "Sex")))

write.csv(mexico_vs_other, "outputs/tables/mexico_vs_other.csv",
          row.names = FALSE)

cat("\n--- Mexico vs all other origins (immigrants, age 25+) ---\n")
print(mexico_vs_other, row.names = FALSE, digits = 4)

# Metropolitan concentration: the volume claim (which metro holds the most
# Mexican-born observations) and the share claim (which metros are most
# overwhelmingly Mexican among their immigrants) are different questions and
# are reported as separate orderings rather than one table sorted two ways.
# MSA == 0 is the residual "not in an identified metropolitan area" bucket,
# not a metro. Left in, it ranks second by volume and would be read as one.
mex_msa <- mex_base[!is.na(MSA) & MSA != 0, .(
  mexican    = sum(Origin == "Mexico"),
  immigrants = .N), by = MSA]
mex_msa[, share := 100 * mexican / immigrants]

if (!is.null(msa_name_lookup)) {
  mex_msa <- merge(mex_msa, msa_name_lookup, by = "MSA", all.x = TRUE)
} else {
  mex_msa[, MSA_Name := NA_character_]
  warning("top175_list_named not in analysis_data.RData; ",
          "metro concentration reported by MSA code without names.")
}

write.csv(mex_msa[order(-mexican)], "outputs/tables/mexico_msa_concentration.csv",
          row.names = FALSE)
cat("\n--- Top 5 metros by Mexican-born volume ---\n")
print(head(mex_msa[order(-mexican)], 5), row.names = FALSE)
cat("\n--- Top 5 metros by Mexican share of immigrants (min 500 immigrants) ---\n")
print(head(mex_msa[immigrants >= 500][order(-share)], 5), row.names = FALSE)

write.csv(demographic_summary_table,
          "outputs/tables/demographic_summary_table.csv", row.names = FALSE)
write.csv(bl, "outputs/tables/bpld_codes.csv", row.names = FALSE)

# 25_dissertation_table_helper.R loads this checkpoint for
# demographic_summary_table and bpld_lookup. Kept for that reason only --
# nothing in this script reads a checkpoint.
save(bpld_lookup, demographic_summary_table,
     file = file.path(analysis_checkpoint_dir, "after_sec22.RData"))

cat("\n", strrep("-", 70), "\n", sep = "")
cat("Saved: outputs/tables/parent_educ_missingness.csv\n")
cat("       outputs/tables/educ_by_origin.csv\n")
cat("       outputs/tables/demographic_summary_table.csv\n")
cat("       outputs/tables/bpld_codes.csv\n")
cat(sprintf("       %s (for 24)\n",
            file.path(analysis_checkpoint_dir, "after_sec22.RData")))
cat("22_origin_composition.R complete.\n")
