# =============================================================================
# 24_results_tables.R  --  All results and robustness tables (LaTeX)
# Renumbered 2026-10-21 from 14_results_tables.R (archive-and-renumber pass).
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# REWRITTEN 2026-10-05 to the Chapter 1 table convention, after comparing
# against the guns chapter's tables directly.
#
# WHY THE PREVIOUS VERSION WAS WRONG
# The 2026-10-04 rewrite got the WRAPPER right -- \caption and \label before
# the table, \resizebox, the \minipage notes block -- but not the table itself.
# Chapter 1 follows the standard economics convention and the previous version
# did not:
#
#   specifications are NUMBERED COLUMNS, coefficients run down the rows
#   standard errors sit on their OWN ROW beneath the estimate, not stacked
#     inside the cell with \makecell
#   stars are RAISED math superscripts, $^{***}$, not bare ***
#   a \midrule separates the coefficient block from a DIAGNOSTICS PANEL --
#     fixed effects, clustering, first stage, observations
#
# Set beside Chapter 1, the old layout read as a different author's work. The
# \makecell dependency is gone with it; the preamble no longer needs it for
# these tables.
#
# A PURE RENDERING STEP
# Reads CSVs from outputs/tables/ and writes .tex. Loads no checkpoint, fits
# nothing. If a number here is wrong, it is wrong in the script that wrote the
# CSV.
#
#   table_education_main.tex    <- main_results.csv            (29)
#   table_migration_iv.tex      <- results_migration_iv.csv    (12)
#   table_wage_iv.tex           <- results_wage_iv.csv         (13)
#   table_exposure_missing.tex  <- exposure_missing_check.csv  (17)   [NEW]
#   table_birth_cohort.tex      <- birth_cohort_final_spec.csv (18b)  [NEW]
#   table_adao_placebo.tex      <- adao_placebo_results.csv    (20)   [NEW]
#
# THE FOUR NEW TABLES EXIST BECAUSE OF A STANDING RULE: the writeup never
# sends the reader to the code. Any claim drawn from the data must be supported
# by a table or figure in the paper. The robustness discussion makes claims
# about the identifying subsample, the cohort gradient and the placebo
# diagnostics; those claims need tables, and these are they.
#
# Descriptive tables are NOT here. Table 1 (estimation sample) comes from 23
# and the origin-composition tables from 21; both currently emit CSV only and
# need a home before the Data section can cite them.
#
# PREAMBLE: booktabs, graphicx. No xcolor (this style does not stripe) and no
# makecell (no stacked cells any more).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")
source("scripts/analysis/00b_table_helper.R")

library(data.table)

OUT <- "outputs/tables"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

need <- function(f) {
  p <- file.path(OUT, f)
  if (!file.exists(p))
    stop("Missing ", p, "\nRun the script that produces it first -- see this ",
         "file's header for the mapping.")
  as.data.table(read.csv(p, check.names = FALSE))
}

# --- formatting -------------------------------------------------------------
# Threshold fix retained from 2026-09-21: *** is p < 0.01, matching every notes
# block below. Raised math superscripts, per Chapter 1.
stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.01, "$^{***}$",
      ifelse(p < 0.05, "$^{**}$",
        ifelse(p < 0.10, "$^{*}$", ""))))
}
# write_dissertation_table() passes escape = FALSE, which is required for the
# math and the raised stars -- but it also means NOTHING is escaped, including
# text lifted from a CSV. Origin labels and group names go through this.
tex_txt <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  gsub("([&%$#_{}])", "\\\\\\1", x)
}
sci  <- function(x, d = 2) formatC(x, format = "e", digits = d)
fx   <- function(x, d = 3) formatC(x, format = "f", digits = d)
cm   <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
est  <- function(x, p = NULL, d = 2, f = sci)
  paste0(f(x, d), if (is.null(p)) "" else stars(p))
se   <- function(x, d = 2, f = sci) paste0("(", f(x, d), ")")

# Row constructors. A row is a character vector: label, then one entry per
# specification column. Empty-label rows carry standard errors, as in Ch1.
rw   <- function(label, ...) c(label, ...)
blank <- function(ncol) c("", rep("", ncol))
MID  <- function(ncol) c("@MIDRULE@", rep("", ncol))

build <- function(rows, colnames_) {
  m <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(m) <- colnames_
  m
}

STARS_NOTE <- "$^{***}$ $p<0.01$, $^{**}$ $p<0.05$, $^{*}$ $p<0.10$."
CLUSTER_NOTE <- paste(
  "Standard errors in parentheses are clustered on metropolitan area.")

# =============================================================================
# TABLE 1 -- Education, the headline specification
# =============================================================================
m <- need("main_results.csv")
setkey(m, outcome)
h  <- m[spec %like% "^HEADLINE"][1]     # y_HSplus
s2 <- m[spec %like% "^SECONDARY"][1]    # y_CollPlus
r3 <- m[spec %like% "^ROBUSTNESS"][1]   # y_HSplus + race controls

# The identifying count comes from 17's restricted row rather than being
# asserted: it is the same quantity, computed once, and this keeps the two
# tables arithmetically consistent.
xm <- need("exposure_missing_check.csv")
n_identifying <- xm[sample %like% "^restricted", n][1]


edu <- build(list(
  rw("Reduced form",  est(h$reduced_form, h$rf_p), est(s2$reduced_form, s2$rf_p),
                      est(r3$reduced_form, r3$rf_p)),
  rw("\\quad MSA-clustered SE",   se(h$rf_se),  se(s2$rf_se),  se(r3$rf_se)),
  rw("\\quad State-clustered SE", se(h$rf_se_state), se(s2$rf_se_state),
                                   se(r3$rf_se_state)),
  rw("2SLS",          est(h$tsls, h$tsls_p), est(s2$tsls, s2$tsls_p),
                      est(r3$tsls, r3$tsls_p)),
  rw("\\quad MSA-clustered SE",   se(h$tsls_se), se(s2$tsls_se), se(r3$tsls_se)),
  rw("\\quad State-clustered SE", se(h$tsls_se_state), se(s2$tsls_se_state),
                                   se(r3$tsls_se_state)),
  rw("Effect per SD (pp)", fx(h$pp_per_sd, 2), fx(s2$pp_per_sd, 2),
                           fx(r3$pp_per_sd, 2)),
  MID(3),
  rw("First-stage coefficient",   fx(h$first_stage, 3), fx(s2$first_stage, 3),
                      fx(r3$first_stage, 3)),
  rw("\\quad MSA-clustered SE", se(h$fs_se, 3, fx), se(s2$fs_se, 3, fx),
                                 se(r3$fs_se, 3, fx)),
  rw("Year FE",       "Yes", "Yes", "Yes"),
  rw("MSA FE",        "Yes", "Yes", "Yes"),
  rw("Race controls", "Base", "Base", "Extended"),
  rw("Observations",  cm(h$n), cm(s2$n), cm(r3$n)),
  rw("Identifying observations", cm(n_identifying), cm(n_identifying),
                                 cm(n_identifying)),
  rw("MSA clusters",   "175", "175", "175"),
  rw("State clusters", "45",  "45",  "45")
), c(" ", "(1)", "(2)", "(3)"))

write_dissertation_table(
  edu,
  caption = paste("Immigrant Exposure During Formative Years and U.S.-Born",
                  "Educational Attainment"),
  label = "tab:education_main",
  align = "lccc",
    notes = paste(
      "Column (1) reports the high-school margin: completing high school or more",
      "against failing to complete high school. Column (2) reports the college",
      "margin: completing some college or above against all lower attainment.",
      "Column (3) repeats column (1) with indicators for Asian and Aboriginal",
      "identity added to the control set. The instrument is the 1990",
      "origin-group settlement share interacted with national arrival counts,",
      "averaged over each individual's ages 15--25 window; the endogenous",
      "variable is the realized inflow averaged over the \\emph{same} window.",
      "Controls are race, sex, marital status, child age group, poverty status",
      "and parental education. Standard errors in parentheses are clustered on",
      "metropolitan area; state-clustered standard errors are reported alongside",
      "for every coefficient so that their ratio can be formed where it is",
      "wanted. ``Effect per SD'' rescales the reduced form to a",
      "one-standard-deviation increase in exposure among identifying",
      "observations and is the interpretable magnitude; the raw coefficients are",
      "per unit of inflow. Identifying",
      "observations are those with non-zero exposure; the remainder are retained",
      "and contribute to the controls and fixed effects but not to the",
      "instrument's variation.",
      STARS_NOTE),
  file = file.path(OUT, "table_education_main.tex")
)

# =============================================================================
# TABLE 2 -- Internal migration
# =============================================================================
g <- need("results_migration_iv.csv")
g[, grp := sub("\\s*\\(Spec1:.*$", "", Group)]
a <- g[grp %like% "All"][1]; mn <- g[grp %like% "Men"][1]; wm <- g[grp %like% "Women"][1]

mig <- build(list(
  rw("Levels (Card specification)", est(a$Est_Card, a$p_Card),
     est(mn$Est_Card, mn$p_Card), est(wm$Est_Card, wm$p_Card)),
  rw("", se(a$SE_Card_MSA), se(mn$SE_Card_MSA), se(wm$SE_Card_MSA)),
  rw("Differences", est(a$Est_diff, a$p_diff), est(mn$Est_diff, mn$p_diff),
     est(wm$Est_diff, wm$p_diff)),
  rw("", se(a$SE_diff_MSA), se(mn$SE_diff_MSA), se(wm$SE_diff_MSA)),
  MID(3),
  rw("Mean P(moved)", fx(a$Mean_P_moved, 4), fx(mn$Mean_P_moved, 4),
     fx(wm$Mean_P_moved, 4)),
  rw("First-stage Wald, levels", fx(a$FS_F_Card, 1), fx(mn$FS_F_Card, 1),
     fx(wm$FS_F_Card, 1)),
  rw("First-stage Wald, differences", fx(a$FS_F_diff, 1), fx(mn$FS_F_diff, 1),
     fx(wm$FS_F_diff, 1)),
  rw("SE clustering", "MSA", "MSA", "MSA"),
  rw("Observations, levels", cm(a$N_Spec1), cm(mn$N_Spec1), cm(wm$N_Spec1)),
  rw("Observations, differences", cm(a$N_Spec2), cm(mn$N_Spec2), cm(wm$N_Spec2))
), c(" ", "All", "Men", "Women"))

write_dissertation_table(
  mig,
  caption = "IV Estimates of Immigrant Inflows on U.S.-Born Internal Migration",
  label = "tab:migration_iv",
  align = "lccc",
  notes = paste(
    "The two specifications use different samples, hence the separate",
    "observation counts. The levels specification follows",
    "\\citet{card2001immigrant} on a seeded random draw from U.S.-born",
    "respondents aged 16--68. The differences specification uses the full sample",
    "of U.S.-born respondents aged 25--68; its instrument differences a baseline",
    "window at ages 14--19 against the preceding five years, and those windows",
    "overlap below age 24, which would drive the instrument mechanically toward",
    "zero in a younger sample. Standard errors in parentheses are clustered on",
    "metropolitan area; state-clustered alternatives do not change any",
    "conclusion. The first-stage Wald statistic is computed against the same",
    "metropolitan-area-clustered covariance matrix as the standard errors rather",
    "than as a homoskedastic $F$: the instrument varies at the",
    "metropolitan-area-by-year level, so a statistic assuming independence",
    "across individual records would largely restate the number of records",
    "rather than measure instrument strength.",
    STARS_NOTE),
  file = file.path(OUT, "table_migration_iv.tex")
)

# =============================================================================
# TABLE 3 -- Wages by occupation group
# =============================================================================
w <- need("results_wage_iv.csv")
w[, Occup_Group := gsub("&", "\\\\&", Occup_Group)]
grps <- unique(w$Occup_Group)
pick <- function(occ, who, col) w[Occup_Group == occ & Group == who][[col]][1]

wage_rows <- list()
for (i in seq_along(grps)) {
  o <- grps[i]
  lab <- if (grepl("^Military", o)) paste0(o, "$^{\\dagger}$") else o
  wage_rows[[length(wage_rows) + 1L]] <- rw(lab,
    est(pick(o, "All", "Estimate"), NULL), est(pick(o, "Men", "Estimate"), NULL),
    est(pick(o, "Women", "Estimate"), NULL))
  # Stars come from 13's own Stars column rather than being recomputed, so the
  # table cannot disagree with the script that produced the estimates.
  wage_rows[[length(wage_rows)]][2:4] <- paste0(
    wage_rows[[length(wage_rows)]][2:4],
    vapply(c("All", "Men", "Women"), function(k) {
      s <- pick(o, k, "Stars"); if (is.na(s) || !nzchar(trimws(s))) "" else
        paste0("$^{", trimws(s), "}$") }, character(1)))
  wage_rows[[length(wage_rows) + 1L]] <- rw("",
    se(pick(o, "All", "SE")), se(pick(o, "Men", "SE")), se(pick(o, "Women", "SE")))
  wage_rows[[length(wage_rows) + 1L]] <- rw("\\quad Effect per SD (\\$/hour)",
    fx(pick(o, "All", "Est_per_SD"), 3), fx(pick(o, "Men", "Est_per_SD"), 3),
    fx(pick(o, "Women", "Est_per_SD"), 3))
  if (i < length(grps)) wage_rows[[length(wage_rows) + 1L]] <- blank(3)
}

wage <- build(wage_rows, c("Occupation Group", "All", "Men", "Women"))

write_dissertation_table(
  wage,
  caption = paste("IV Estimates of Immigrant Inflows on Mean Hourly Wages by",
                  "Occupation Group"),
  label = "tab:wage_iv",
  align = "lccc",
  notes = paste(
    "Coefficients are IV estimates within",
    "metropolitan-area-by-occupation-by-year cells, with metropolitan-area and",
    "survey-year fixed effects; cell means are unweighted. Occupation-specific",
    "inflows are divided by the number of survey waves that could have observed",
    "each arrival year, as elsewhere in the chapter, and arrival years with fewer",
    "than two eligible waves are dropped.",
    "Standard errors in parentheses are clustered on metropolitan area. The",
    "per-standard-deviation column is the comparable one, since inflow volume",
    "differs by an order of magnitude across occupation groups and the raw",
    "coefficients are therefore measured against different-sized shocks.",
    "$^{\\dagger}$~Military occupations are reported for completeness but are",
    "not interpretable; the reasons are given in the text.",
    STARS_NOTE),
  file = file.path(OUT, "table_wage_iv.tex")
)

# =============================================================================
# TABLE 4 -- Pooled versus identifying sample  [NEW]
# =============================================================================
p_hs <- xm[outcome == "y_HSplus"   & sample %like% "^pooled"][1]
r_hs <- xm[outcome == "y_HSplus"   & sample %like% "^restricted"][1]
p_cl <- xm[outcome == "y_CollPlus" & sample %like% "^pooled"][1]
r_cl <- xm[outcome == "y_CollPlus" & sample %like% "^restricted"][1]

emx <- build(list(
  rw("Reduced form", est(p_hs$reduced_form, p_hs$rf_p),
     est(r_hs$reduced_form, r_hs$rf_p), est(p_cl$reduced_form, p_cl$rf_p),
     est(r_cl$reduced_form, r_cl$rf_p)),
  rw("", se(p_hs$rf_se), se(r_hs$rf_se), se(p_cl$rf_se), se(r_cl$rf_se)),
  rw("2SLS", est(p_hs$tsls, p_hs$tsls_p), est(r_hs$tsls, r_hs$tsls_p),
     est(p_cl$tsls, p_cl$tsls_p), est(r_cl$tsls, r_cl$tsls_p)),
  rw("", se(p_hs$tsls_se), se(r_hs$tsls_se), se(p_cl$tsls_se), se(r_cl$tsls_se)),
  rw("Effect per SD (pp)", fx(p_hs$pp_per_sd, 2), fx(r_hs$pp_per_sd, 2),
     fx(p_cl$pp_per_sd, 2), fx(r_cl$pp_per_sd, 2)),
  MID(4),
  rw("First-stage coefficient", fx(p_hs$first_stage, 3), fx(r_hs$first_stage, 3),
     fx(p_cl$first_stage, 3), fx(r_cl$first_stage, 3)),
  rw("Outcome mean (per cent)",
     if (is.null(p_hs$outcome_mean)) "---" else fx(p_hs$outcome_mean, 1),
     if (is.null(r_hs$outcome_mean)) "---" else fx(r_hs$outcome_mean, 1),
     if (is.null(p_cl$outcome_mean)) "---" else fx(p_cl$outcome_mean, 1),
     if (is.null(r_cl$outcome_mean)) "---" else fx(r_cl$outcome_mean, 1)),
  rw("Observations", cm(p_hs$n), cm(r_hs$n), cm(p_cl$n), cm(r_cl$n))
), c(" ", "Pooled", "Identifying", "Pooled", "Identifying"))

write_dissertation_table(
  emx,
  caption = paste("Pooled Sample versus Identifying Subsample,",
                  "Both Education Margins"),
  label = "tab:exposure_missing",
  align = "lcccc",
  notes = paste(
    "Columns (1) and (2) report the high-school margin; columns (3) and (4)",
    "the college margin. ``Pooled'' is the full estimation sample.",
    "``Identifying'' restricts to observations with non-zero exposure -- those",
    "whose ages 15--25 window opens in 1990 or later, the base year of the",
    "settlement shares.",
    CLUSTER_NOTE,
    "The outcome mean is the share of each sample with the outcome, which is",
    "what makes a percentage-point effect interpretable.",
    "The restricted estimate is roughly three times the pooled estimate on",
    "both margins. This is not a validation of the pooled specification but a",
    "statement about who identifies it: the condition that creates non-zero",
    "exposure also makes the identifying sample young, and younger respondents",
    "are more likely to be recorded as not having completed a given level of",
    "schooling because they have not yet finished it. The restricted magnitude",
    "is correspondingly less plausible as a causal effect, not more. The",
    "cohort structure behind this is reported in",
    "Table~\\ref{tab:birth_cohort}.",
    STARS_NOTE),
  file = file.path(OUT, "table_exposure_missing.tex")
)

# =============================================================================
# TABLE 5 -- Birth-cohort gradient  [NEW]
# =============================================================================
# Panel B (cumulative restriction) rather than Panel A (disjoint bands): each
# row retains the full older sample, so the trend is read off a nested sequence
# rather than off bands whose standard errors widen mechanically at the edges.
bc <- need("birth_cohort_final_spec.csv")
ref <- bc[panel == "reference"][1]
cum <- bc[panel == "B: cumulative"][order(-n)]

cols <- c(ref$band, cum$band)
getv <- function(f) c(ref[[f]], cum[[f]])

bcx <- build(list(
  rw("Reduced form", mapply(est, getv("reduced_form"), getv("rf_p"))),
  rw("", vapply(getv("rf_se"), se, character(1))),
  rw("Effect per SD (pp)", vapply(getv("pp_per_sd"), fx, character(1), d = 2)),
  MID(length(cols)),
  rw("First-stage coefficient", vapply(getv("first_stage"), fx, character(1), d = 3)),
  rw("Mean age", vapply(getv("mean_age"), fx, character(1), d = 1)),
  rw("Per cent non-high-school", vapply(getv("pct_non_HS"), fx, character(1), d = 1)),
  rw("Observations", vapply(getv("n"), cm, character(1)))
), c(" ", gsub("born ", "", cols)))

write_dissertation_table(
  bcx,
  caption = "Estimates by Birth Cohort: Cumulative Restriction to Older Cohorts",
  label = "tab:birth_cohort",
  align = paste0("l", strrep("c", length(cols))),
  notes = paste(
    "The first column repeats the headline specification. Each subsequent column",
    "restricts the sample to birth cohorts at or before the stated year, so the",
    "columns are nested and each retains the full older sample. Cohorts born",
    "before 1975 carry no exposure variation by construction -- their ages",
    "15--25 window opens before the 1990 base year -- and cannot be estimated.",
    "Read each estimate against the mean age and the non-high-school rate in the",
    "same column. Standard errors in parentheses are clustered on metropolitan",
    "area.",
    STARS_NOTE),
  file = file.path(OUT, "table_birth_cohort.tex")
)

# =============================================================================
# TABLE 6 -- Rotemberg weights and the sensitivity to dropping origins  [NEW]
# =============================================================================
# Two tables rather than one. The weight decomposition is indexed by origin and
# the sensitivity by how many origins are removed, so they do not share a
# column structure and forcing them together would produce a header that fits
# neither.
#
# BOTH RUN ON 19's 500,000-row draw, not the full estimation sample, so the
# reference estimate here (1.107e-05) is close to but not identical with the
# headline 2SLS in Table~\ref{tab:education_main}. The notes say so; the
# sensitivity is a percentage change against its own reference, which is the
# comparison that matters.
rw_top <- need("rotemberg_top_groups.csv")
rw_all <- need("rotemberg_results.csv")
rw_top <- merge(rw_top, rw_all[, .(BPLD, beta_k)], by = "BPLD", all.x = TRUE)
setorder(rw_top, -alpha_k)

N_SHOW <- 10L
tw <- head(rw_top, N_SHOW)
tw[, cum := cumsum(alpha_k)]

# Weight counts, computed rather than asserted. 82 of the 230 origins carry
# exactly zero weight -- they have no 1990 settlement to interact with a shock
# -- so a claim about "230 origins" overstates how many are doing any work.
n_all  <- nrow(rw_all)
n_nz   <- rw_all[alpha_k != 0, .N]
n_pos  <- rw_all[alpha_k  > 0, .N]
n_neg  <- rw_all[alpha_k  < 0, .N]

rot_rows <- lapply(seq_len(nrow(tw)), function(i)
  rw(tex_txt(tw$label[i]), fx(100 * tw$alpha_k[i], 1), fx(100 * tw$cum[i], 1),
     sci(tw$beta_k[i])))
rot_rows <- c(rot_rows, list(
  MID(3),
  rw("Origins in the decomposition", cm(n_all), "", ""),
  rw("\\quad with non-zero weight",  cm(n_nz),  "", ""),
  rw("\\quad positive",              cm(n_pos), "", ""),
  rw("\\quad negative",              cm(n_neg), "", "")))

write_dissertation_table(
  build(rot_rows, c("Origin", "Weight (per cent)", "Cumulative (per cent)",
                    "$\\beta_k$")),
  caption = "Rotemberg Weights and Origin-Specific Estimates, Ten Largest Origins",
  label = "tab:rotemberg",
  align = "lccc",
  notes = paste(
    "\\citet{goldsmith2020bartik} show that a shift-share estimator is a",
    "weighted average of just-identified estimators, one per exposure share,",
    "with weights $\\alpha_k$ summing to one. The design is only as credible",
    "as the exogeneity of the shares carrying substantial weight, so the",
    "distribution of those weights is reported directly.",
    "$\\beta_k$ is the just-identified estimate using origin $k$'s share alone.",
    "Weights are computed on a 500{,}000-observation random draw from the",
    "estimation sample.",
    "Most origins contribute nothing: those with no 1990 settlement in any",
    "sample metropolitan area have a share of zero and therefore a weight of",
    "zero, which is why the count with non-zero weight is far below the total.",
    "Of those that do contribute, the two with negative weight carry",
    "$-3.0 \\times 10^{-6}$ and $-3.3 \\times 10^{-7}$ between them, which is",
    "negligible against a total of one.",
    "Table~\\ref{tab:rotemberg_sensitivity} reports what happens when the",
    "highest-weighted origins are removed."),
  file = file.path(OUT, "table_rotemberg.tex")
)

# --- Sensitivity ------------------------------------------------------------
rs <- need("rotemberg_sensitivity.csv")
setorder(rs, drop_n)
lab <- c("0" = "None (full instrument)", "1" = "Largest origin",
         "3" = "Three largest", "5" = "Five largest", "10" = "Ten largest")

sens_rows <- list(
  rw("Estimate", vapply(rs$beta_reduced, sci, character(1))),
  rw("\\quad MSA-clustered SE",
     vapply(rs$se_reduced, function(x) if (is.na(x)) "---" else se(x),
            character(1))),
  rw("Change from full instrument (per cent)",
     vapply(rs$pct_change, function(x) if (x == 0) "---" else
            sprintf("%+.2f", x), character(1))),
  MID(nrow(rs)),
  rw("Rotemberg weight removed (per cent)",
     vapply(100 * rs$weight_dropped, fx, character(1), d = 1)))

write_dissertation_table(
  build(sens_rows, c("Origins dropped", unname(lab[as.character(rs$drop_n)]))),
  caption = paste("Sensitivity of the Estimate to Removing the",
                  "Highest-Weighted Origins"),
  label = "tab:rotemberg_sensitivity",
  align = paste0("l", strrep("c", nrow(rs))),
  notes = paste(
    "Each column re-estimates the specification with the stated origins",
    "excluded from the construction of the instrument, on the same",
    "500{,}000-observation draw used for Table~\\ref{tab:rotemberg}. The",
    "reference column repeats the estimate on the full instrument, so the",
    "percentage changes are measured against a like-for-like baseline rather",
    "than against the full-sample headline.",
    "Removing the largest origin alone moves the estimate by 3.3 per cent, and",
    "removing the ten largest -- which together carry 56.6 per cent of the",
    "identifying weight -- moves it by 4.5 per cent. The estimate is therefore",
    "not an artefact of any one origin's exposure shares, which is the",
    "specific fragility the Rotemberg decomposition is designed to detect.",
    "The changes are not monotone in the number of origins removed. With",
    "roughly half the identifying weight gone the estimate is determined by a",
    "different and much smaller set of shares, and small non-monotonicities",
    "across those columns should not be read as structure."),
  file = file.path(OUT, "table_rotemberg_sensitivity.tex")
)

# =============================================================================
# TABLE 7 -- Metropolitan-area attribution in the actual-inflow series  [NEW]
# =============================================================================
# IPUMS reports residence AT SURVEY, not at arrival, so the actual-inflow
# series credits an arrival to the metro where the respondent was eventually
# enumerated. The instrument is unaffected -- national arrivals by origin do
# not depend on where anyone later moves -- so only the endogenous variable
# carries the error, and measurement error there is absorbed by the first stage
# UNLESS it correlates with the instrument.
#
# The three rows are the three things that would have to be true for the error
# to matter, each measured on the arrival years where the survey calendar
# permits a near-arrival comparison. If 34 has not been re-run the table is
# skipped rather than built from stale input.
att_path <- file.path(OUT, "attribution_summary.csv")
if (!file.exists(att_path)) {
  cat("\n  SKIPPED table_attribution.tex -- run 19_attribution_stability_check.R\n")
} else {
  at <- as.data.table(read.csv(att_path, check.names = FALSE))
  getv <- function(pat) at[grepl(pat, quantity), value][1]

  att <- build(list(
    rw("Correlation of metro shares, near-arrival against all",
       fx(getv("^Correlation of within-year"), 4)),
    rw("Mean absolute discrepancy per cell (pp)",
       fx(getv("^Mean absolute"), 4)),
    MID(1),
    rw("Gap against 1990 share index, across metros",
       fx(getv("signed relative gap"), 4)),
    rw("Gap against instrument, within metro",
       fx(getv("within metro"), 4)),
    rw("\\quad mean of per-metro correlations",
       fx(getv("^Mean of per-metro"), 4)),
    rw("\\quad metros beyond $+0.3$ (per cent)",
       fx(getv("^Metros with correlation"), 1)),
    rw("\\quad expected beyond $+0.3$ by chance (per cent)",
       fx(getv("^Expected beyond"), 1))
  ), c(" ", "Value"))

  write_dissertation_table(
    att,
    caption = paste("Metropolitan-Area Attribution in the Actual-Inflow Series,",
                    "Arrivals", at$arrival_years[1]),
    label = "tab:attribution",
    align = "lc",
    notes = paste(
      "The rows are ordered by increasing relevance to the design. The first two",
      "establish magnitude, comparing the distribution of arrivals across",
      "metropolitan areas measured from all respondents against the same",
      "distribution measured only from respondents observed within two years of",
      "arrival. The third asks whether the discrepancy is larger in metropolitan",
      "areas with large 1990 co-national communities. The fourth is the one that",
      "binds, because metropolitan-area fixed effects absorb anything constant",
      "within a metropolitan area: it asks whether the discrepancy moves with the",
      "instrument across arrival years \\emph{within} a metropolitan area. The",
      "per-metropolitan-area tail count is reported against the rate chance alone",
      "produces at this number of arrival years, since a correlation of $0.3$ on",
      "roughly a dozen observations sits well inside sampling noise. The",
      "comparison is restricted to arrivals in 2006--2019, the years for which the",
      "survey calendar supplies a wave within two years of every arrival year."),
    file = file.path(OUT, "table_attribution.tex")
  )
}

# =============================================================================
# TABLE 8 -- Shock-randomization placebo  [NEW]
# =============================================================================
ad <- need("adao_placebo_results.csv")
setkey(ad, outcome)
ord <- c("Black", "Asian", "Aboriginal")
ad  <- ad[match(ord, outcome)]

adx <- build(list(
  rw("Estimate", mapply(est, ad$estimate, ad$msa_p)),
  rw("\\quad MSA-clustered SE", vapply(ad$msa_se, se, character(1))),
  rw("\\quad State-clustered SE", vapply(ad$state_se, se, character(1))),
  MID(nrow(ad)),
  rw("$p$-value, MSA clustering", vapply(ad$msa_p, fx, character(1), d = 4)),
  rw("$p$-value, state clustering", vapply(ad$state_p, fx, character(1), d = 4)),
  rw("Magnitude, per cent of main estimate",
     vapply(ad$pct_of_reference, fx, character(1), d = 1))
), c(" ", ad$outcome))

write_dissertation_table(
  adx,
  caption = paste("Shock-Randomization Placebo: Outcomes Determined at Birth"),
  label = "tab:adao_placebo",
  align = paste0("l", strrep("c", nrow(ad))),
  notes = paste(
    "Each column regresses an indicator determined at birth on the instrument,",
    "with the full control set and fixed effects of the main specification.",
    "Because these outcomes cannot be caused by immigrant exposure, a",
    "significant coefficient indicates the instrument is correlated with",
    "something other than the channel of interest. The final row expresses each",
    "placebo estimate as a percentage of the main educational-attainment",
    "estimate, which is the relevant scale, since statistical significance on a",
    "sample this size is not itself informative about magnitude. Both clustering",
    "schemes are reported so that their ratio can be formed directly.",
    STARS_NOTE),
  file = file.path(OUT, "table_adao_placebo.tex")
)

# =============================================================================
# TABLE 11 -- Multinomial decomposition of the binary estimate
# =============================================================================
# ADDED 2026-10-21. Produced by 14b_multinomial_check.R on a 100,000-observation
# random draw, not the full sample: the full-sample multinomial is a ten-day
# run and is not worth it for a decomposition that, as the table shows, carries
# no size information the binary does not.
#
# Bootstrap standard errors are deliberately NOT reported. 14b uses 99
# MSA-resampled replications, which it describes as adequate for a keep-or-cut
# decision and not for reporting. The table therefore gives point estimates and
# the proportionality comparison, which is what the result turns on.
mm <- need("multinomial_final_spec_test.csv")
if (!all(c("pop_share", "effect_share") %in% names(mm))) {
  cat("\nSKIPPED table_multinomial: re-run 14b_multinomial_check.R --",
      "the CSV predates the pop_share / effect_share columns.\n")
} else {
  lab <- c(non_HS = "No high school", HighSch = "High school",
           College = "Some college",   Grad = "Graduate degree")
  pick <- function(cat, col) {
    v <- mm[[col]][match(cat, mm$category)]
    if (length(v) && !is.na(v)) v else NA_real_
  }
  rows <- list()
  rows[[1]] <- rw(lab[["non_HS"]], fx(pick("non_HS", "pp_per_sd"), 2), "---", "---")
  rows[[2]] <- MID(3)
  for (k in c("HighSch", "College", "Grad"))
    rows[[length(rows) + 1L]] <- rw(lab[[k]], fx(pick(k, "pp_per_sd"), 2),
                                    fx(pick(k, "effect_share"), 1),
                                    fx(pick(k, "pop_share"), 1))
  rows[[length(rows) + 1L]] <- MID(3)
  implied <- sum(mm$pp_per_sd[mm$category %in% c("HighSch", "College", "Grad")])
  rows[[length(rows) + 1L]] <- rw("Implied high school or above",
                                  fx(implied, 2), "100.0", "100.0")
  rows[[length(rows) + 1L]] <- rw("Binary high school or above",
                                  fx(pick("binary_HSplus", "pp_per_sd"), 2),
                                  "---", "---")

  mnt <- build(rows, c(" ", "Effect (pp per SD)", "Share of effect",
                       "Share of population"))

  write_dissertation_table(
    mnt,
    caption = paste("Multinomial Decomposition of the Binary Estimate,",
                    "100,000-Observation Draw"),
    label = "tab:multinomial",
    notes = paste(
      "Two-stage residual inclusion with a four-category outcome. Multinomial",
      "coefficients across log-odds equations are not comparable to a binary",
      "model's, so both models are converted to average marginal effects on",
      "predicted probabilities and rescaled to a one-standard-deviation increase",
      "in exposure. The four marginal effects sum to zero by construction.",
      "``Share of effect'' is each completing category's share of the total",
      "movement out of no-high-school; ``share of population'' is its share of",
      "the completing population in the same draw. Standard errors are not",
      "reported: the bootstrap behind the estimate uses 99 metropolitan-area",
      "resamples, which is too coarse to report."),
    file = file.path(OUT, "table_multinomial.tex")
  )
}

# =============================================================================
cat("\nTen tables written in the Chapter 1 convention:\n")
for (f in c("table_education_main", "table_migration_iv", "table_wage_iv",
            "table_exposure_missing", "table_birth_cohort", "table_rotemberg",
            "table_rotemberg_sensitivity", "table_attribution",
            "table_adao_placebo", "table_multinomial"))
  cat("  ", file.path(OUT, paste0(f, ".tex")), "\n", sep = "")
cat("\nPreamble: booktabs, graphicx. No xcolor, no makecell.\n")
cat("Descriptive tables (Table 1 from 23, origin composition from 21) are not\n")
cat("rendered here and still emit CSV only.\n")
cat("24_results_tables.R complete.\n")
