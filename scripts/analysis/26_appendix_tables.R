# =============================================================================
# 26_appendix_tables.R  --  Appendix tables for Chapter 2
# Renumbered 2026-10-21 from 24b_appendix_tables.R (archive-and-renumber pass).
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renders the three documentation tables the Data section references:
#   table_extract_composition.tex   from 08's extract_composition.csv
#   table_mexico_vs_other.tex       from 21's mexico_vs_other.csv
#   table_parent_educ_missing.tex   from 21's parent_educ_missingness.csv
# =============================================================================
#
# Reads the CSVs rather than the microdata, so this is a pure rendering step
# like 14: it needs no checkpoint, holds no estimation, and can be re-run on
# its own after 08 and 21 without reloading analysis_data.RData. Run order is
# 08 -> 21 -> 24 -> 24b.
#
# The .tex files were first written by hand on 2026-10-19 to avoid another
# pipeline run; this script reproduces them byte-for-byte from the same CSVs
# so they are not orphaned the way the July Mexico breakdowns were.

source("scripts/analysis/00_setup_analysis.R")
source("scripts/analysis/00b_table_helper.R")

stopifnot(exists("write_dissertation_table"))
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

req <- c("outputs/tables/extract_composition.csv",
         "outputs/tables/demographic_summary_table.csv",
         "outputs/tables/mexico_vs_other.csv",
         "outputs/tables/parent_educ_missingness.csv")
missing <- req[!file.exists(req)]
if (length(missing))
  stop("Missing input(s):\n  ", paste(missing, collapse = "\n  "),
       "\nRun 08 and 21 first.")

fmt_n   <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_pct <- function(x) sprintf("%.1f", x)

# --- A. Composition of the analysis extract ----------------------------------
xc <- as.data.table(read.csv(req[1], stringsAsFactors = FALSE))
ds <- as.data.table(read.csv(req[2], stringsAsFactors = FALSE))

# Panel labels. Without them the table is six rows reading
# Immigrant/U.S.-born/Total twice, and only the notes say which is which.
extract_tbl <- rbind(
  data.frame(Group = "\\emph{Full analysis extract}", Observations = "", Pct = ""),
  data.frame(Group = c("\\quad Immigrant", "\\quad U.S.-born"),
             Observations = fmt_n(c(xc[group == "Immigrant", obs],
                                    xc[group == "U.S.-born", obs])),
             Pct = fmt_pct(c(xc[group == "Immigrant", share],
                             xc[group == "U.S.-born", share]))),
  data.frame(Group = "\\quad Total", Observations = fmt_n(sum(xc$obs)), Pct = "100.0"),
  data.frame(Group = "@MIDRULE@", Observations = "", Pct = ""),
  data.frame(Group = "\\emph{Ages 25 and over}", Observations = "", Pct = ""),
  data.frame(Group = c("\\quad Immigrant", "\\quad U.S.-born", "\\quad Total"),
             Observations = fmt_n(ds[match(c("All Immigrants", "U.S.-born",
                                             "Combined"), Group), N]),
             Pct = fmt_pct(100 * ds[match(c("All Immigrants", "U.S.-born",
                                            "Combined"), Group), N] /
                           ds[Group == "Combined", N])))
colnames(extract_tbl) <- c("Group", "Observations", "Per cent")

write_dissertation_table(
  extract_tbl,
  caption = paste("Composition of the Analysis Extract, Before and After the",
                  "Age Restriction"),
  label = "tab:extract_composition",
  notes = paste(
    "The upper panel is the full analysis extract; the lower panel applies the",
    "age 25 and over restriction used for the descriptive statistics in",
    "Section~\\ref{sec:immigration_data}. Shares are of each panel's own total.",
    "Immigrants are overrepresented by design: the extract retains every",
    "foreign-born record and a 25 per cent within-year draw of the U.S.-born,",
    "following \\citet{card2001immigrant}, so these shares describe the sample",
    "and not the population."),
  file = "outputs/tables/table_extract_composition.tex")

# --- B. Mexico against all other origins -------------------------------------
mv <- as.data.table(read.csv(req[3], stringsAsFactors = FALSE,
                             check.names = FALSE))
pick <- function(dim, origin, cat) {
  v <- mv[Dimension == dim & Origin == origin & Category == cat, share]
  if (length(v) == 1L) v else NA_real_
}
pair <- function(dim, cat, lab) data.frame(
  Row = lab, Mexico = fmt_pct(pick(dim, "Mexico", cat)),
  Other = fmt_pct(pick(dim, "Other", cat)))

EMP <- "Occupation (employed)"
if (!EMP %in% mv$Dimension)
  stop("mexico_vs_other.csv has no '", EMP, "' rows. Re-run 21.")

educ <- rbindlist(Map(function(c, l) pair("Education", c, l),
  c("non_HS", "HighSch", "College", "Grad"),
  c("No high school", "High school", "Some college", "Graduate degree")))

occ_cats <- c("Low-Skill Labor & Service", "Operatives & Craft Workers",
              "Clerical Workers", "Sales Workers", "Managers",
              "Professional & Technical Workers")
occ <- rbindlist(lapply(occ_cats, function(c)
  pair(EMP, c, gsub("&", "\\\\&", c))))

mex_tbl <- rbind(educ,
                 data.frame(Row = "@MIDRULE@", Mexico = "", Other = ""),
                 occ,
                 data.frame(Row = "@MIDRULE@", Mexico = "", Other = ""),
                 pair("Sex", "1", "Male"))
colnames(mex_tbl) <- c("", "Mexico", "All other origins")

write_dissertation_table(
  mex_tbl,
  caption = paste("Mexican-Born Immigrants Against All Other Origins,",
                  "Ages 25 and Over"),
  label = "tab:mexico_vs_other",
  notes = paste(
    "Columns are per cent within origin group, so neither is affected by the",
    "sampling design's overrepresentation of the foreign-born. The upper panel",
    "is educational attainment, the middle panel occupation and the lower panel",
    "sex. Occupation shares are over those who worked in the year observed:",
    "\\texttt{Occup\\_Group} assigns Not in Labor Force on the basis of",
    "\\texttt{WORKEDYR}, which places roughly a quarter of each origin group",
    "outside the occupational distribution, and they are excluded here so the",
    "shares describe workers. Military occupations are omitted from the panel;",
    "they account for under 0.2 per cent of either group."),
  file = "outputs/tables/table_mexico_vs_other.tex")

# --- C. Parental education missingness ---------------------------------------
pe <- as.data.table(read.csv(req[4], stringsAsFactors = FALSE))
nice <- c(Mom_Educ = "Mother's education (\\texttt{EDUC\\_MOM})",
          Dad_Educ = "Father's education (\\texttt{EDUC\\_POP})")
par_tbl <- data.frame(Variable = unname(nice[pe$Variable]),
                      Pct = fmt_pct(pe$Pct_Missing))
colnames(par_tbl) <- c("Variable", "Per cent missing")

write_dissertation_table(
  par_tbl,
  caption = "Missingness in Parental Education",
  label = "tab:parent_educ_missing",
  notes = paste(
    "Share of the estimation sample for which the variable is unavailable.",
    "Both are coded as an explicit ``Missing'' category rather than imputed or",
    "dropped, and both are retained as controls throughout. The missingness is",
    "structural rather than nonresponse: IPUMS constructs these variables by",
    "linking each individual to a coresident parent's record, so they are",
    "populated only where that parent is in the same household at the time of",
    "the survey."),
  file = "outputs/tables/table_parent_educ_missing.tex")

cat("\n26_appendix_tables.R complete. Written:\n")
cat("  outputs/tables/table_extract_composition.tex\n")
cat("  outputs/tables/table_mexico_vs_other.tex\n")
cat("  outputs/tables/table_parent_educ_missing.tex\n")
