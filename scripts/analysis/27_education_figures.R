# =============================================================================
# 27_education_figures.R  --  Section 15: Education composition figures (Educ_Prop*.png)
# Part of the modular Immigration_Card_Analysis.R pipeline.
# Renumbered 2026-09-07 from 16_education_figures.R.
# Renumbered 2026-10-21 from 15_education_figures.R (archive-and-renumber pass).
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

# UNCHAINED 2026-10-02. This loaded after_sec14.RData, a RESULTS checkpoint,
# for content that is entirely descriptive -- it needed nothing from 14 except
# the objects it then re-saved into after_sec15. Nothing reads after_sec15 any
# more either: 23_msa_stability.R stopped chaining it on 2026-09-30 and now
# takes top175_list_named and all_natives_175 from analysis_data.RData
# directly. Both the load and the save are gone, so this script runs standalone
# and does not depend on the retired chunk loop.

# Section 15: Education composition figures
#
# FIGURE 4 ADDED 2026-10-02 -- education composition by origin country, the
# figure counterpart to Panel C of 22_origin_composition.R. It is built here
# rather than in 21 because 21 owns the tables and this script owns the
# figures, and it recomputes from analysis_data.RData rather than reading 21's
# CSV so neither script depends on the other having run.
#
# Plots the share of each education group over time separately for:
#   - Immigrant arrivals: grouped by year of immigration (YRIMMIG)
#   - U.S.-born: grouped by survey year (YEAR)
# Solid lines = immigrant arrivals; dashed lines = U.S.-born.
#
# Three versions are produced to test the sensitivity of the apparent 2005
# inflection to the age composition of arrivals (see data_section_draft.tex,
# Section "Sensitivity to Age at Arrival"):
#   (a) Educ_Prop.png            -- unrestricted (original figure; immigrant
#                                    side includes arrivals of any age,
#                                    including school-age children who have
#                                    not yet completed their schooling)
#   (b) Educ_Prop_Age25plus.png  -- immigrants restricted to age at arrival
#                                    >= 25; natives restricted to current
#                                    AGE >= 25 for a like-for-like adult
#                                    comparison (same convention as the
#                                    Section 21/22 desc_data restriction)
#   (c) Educ_Prop_Age18to24.png -- immigrants restricted to age at arrival
#                                    18-24 (young adult arrivals); natives
#                                    shown with the same AGE >= 25 reference
#                                    line as (b), since natives have no
#                                    "age at arrival" analog and AGE >= 25
#                                    is the stable adult benchmark against
#                                    which the young-adult-arrival series is
#                                    compared
#
# Age at arrival is computed as AGE - (survey YEAR - YRIMMIG), i.e. current
# age minus years since arrival -- NOT current AGE, which would introduce a
# recency-selection artifact (recent arrival cohorts can only satisfy a
# current-age filter if they were already older on arrival).
# -----------------------------------------------------------------------------

educ_colors <- c(
  "non_HS"         = "#E07B7B",
  "HighSch"        = "#7BB87B",
  "College"        = "#5BB8C8",
  "Grad"           = "#C890D4",
  "College + Grad" = "#333333"
)

# Display labels for the legends. ADDED 2026-10-21. The figures were showing
# the raw Educ_Group factor values -- non_HS, HighSch, Grad -- which are
# codebook names, not English, and are the only place a variable name reached
# a figure in the paper. The underlying factor keeps the original levels so
# nothing downstream changes; only what the legend prints does.
EDUC_LABELS <- c(
  "non_HS"         = "No high school",
  "HighSch"        = "High school",
  "College"        = "Some college",
  "Grad"           = "Graduate degree",
  "College + Grad" = "Some college or more"
)

# Helper: build the immigrant/native education-composition data frame shared
# by all three figure variants. immig_data/native_data should already be
# filtered to the desired age restriction before calling this.
build_educ_prop_data <- function(immig_data, native_data) {
  immig_prop <- immig_data %>%
    filter(immigr == 1, YRIMMIG > 0, !is.na(Educ_Group)) %>%
    group_by(YEAR = YRIMMIG, Educ_Group) %>%
    summarize(n = sum(PERWT, na.rm = TRUE), .groups = "drop") %>%
    group_by(YEAR) %>%
    mutate(share = n / sum(n)) %>%
    ungroup() %>%
    mutate(Group = "Immigrant arrivals")

  native_prop <- native_data %>%
    filter(immigr == 0, !is.na(Educ_Group)) %>%
    group_by(YEAR, Educ_Group) %>%
    summarize(n = sum(PERWT, na.rm = TRUE), .groups = "drop") %>%
    group_by(YEAR) %>%
    mutate(share = n / sum(n)) %>%
    ungroup() %>%
    mutate(Group = "Natives (survey)")

  immig_college_grad <- immig_prop %>%
    filter(Educ_Group %in% c("College", "Grad")) %>%
    group_by(YEAR, Group) %>%
    summarize(share = sum(share), .groups = "drop") %>%
    mutate(Educ_Group = "College + Grad")

  native_college_grad <- native_prop %>%
    filter(Educ_Group %in% c("College", "Grad")) %>%
    group_by(YEAR, Group) %>%
    summarize(share = sum(share), .groups = "drop") %>%
    mutate(Educ_Group = "College + Grad")

  bind_rows(immig_prop, native_prop, immig_college_grad, native_college_grad) %>%
    filter(YEAR >= 1990) %>%
    mutate(
      linetype_var = Group,
      Educ_Group   = factor(Educ_Group,
                            levels = c("non_HS", "HighSch", "College",
                                       "Grad", "College + Grad"))
    )
}

# Helper: render and save one education-composition figure
plot_educ_prop <- function(plot_data, title, subtitle, filename, y_max = 0.55) {
  p <- ggplot(
    plot_data,
    aes(x = YEAR, y = share,
        color    = Educ_Group,
        linetype = linetype_var,
        shape    = linetype_var)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_color_manual(values = educ_colors, labels = EDUC_LABELS,
                       name = "Education") +
    scale_linetype_manual(
      values = c("Immigrant arrivals" = "solid", "Natives (survey)" = "dashed"),
      name   = NULL
    ) +
    scale_shape_manual(
      values = c("Immigrant arrivals" = 16, "Natives (survey)" = 1),
      name   = NULL
    ) +
    # coord_cartesian rather than scale_y_continuous(limits=) -- fixed
    # 2026-09-21. Setting limits on the scale DROPS observations outside the
    # range before plotting, which broke the series wherever a share exceeded
    # y_max and produced a gap implying missing data ("Removed 3 rows
    # containing missing values" on two of the three figures). coord_cartesian
    # zooms the view instead, leaving the data intact and running the line to
    # the panel edge.
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = title, subtitle = subtitle, x = "Year",
         y = "Share of population",
         # ADDED 2026-10-21. Arrivals collapsed during the pandemic, so the
         # 2020-2023 arrival cohorts are thin and their year-to-year swings
         # are sampling noise rather than compositional change. 40_map_origin
         # handles the same problem by stopping at 2019; these figures keep
         # the full range because the collapse is real history, so the
         # caveat is stated instead.
         caption = paste0(
           "Arrival cohorts from 2020 onward are small: immigration fell ",
           "sharply during the pandemic and had not recovered by the end of ",
           "the sample.\nYear-to-year movement at the right-hand end of the ",
           "immigrant series is therefore dominated by sampling variation ",
           "and should not be read as a trend.")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))

  dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, p, width = 10, height = 6, dpi = 300)
  cat(sprintf("%s saved.\n", filename))
  p
}

# reduced_data_iv is not carried through the checkpoint chain past Section 8
# (see 08_load_analysis_subsets.R) -- loaded fresh here directly from
# analysis_data.RData, via a throwaway environment so nothing else already
# in memory gets overwritten. This script now has no other data source.
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
reduced_data_iv <- .tmp_env$reduced_data_iv
rm(.tmp_env)
setDT(reduced_data_iv)

# --- (a) Unrestricted (original figure) ---
plot_data      <- build_educ_prop_data(reduced_data_iv, reduced_data_iv)
educ_prop_plot <- plot_educ_prop(
  plot_data,
  title    = "Education Composition: Immigrant Arrivals vs. U.S. Born",
  subtitle = paste("Immigrant proportions by arrival year (YRIMMIG);",
                   "U.S. born proportions by survey year"),
  filename = "outputs/figures/Educ_Prop.png",
  y_max    = 0.55
)

# --- (b) Age at arrival >= 25 (immigrants); AGE >= 25 (natives) ---
reduced_data_iv[, age_at_arrival := AGE - (as.integer(as.character(YEAR)) - YRIMMIG)]

immig_age25plus  <- reduced_data_iv[immigr == 1 & age_at_arrival >= 25]
native_age25plus <- reduced_data_iv[immigr == 0 & AGE >= 25]

plot_data_age25plus      <- build_educ_prop_data(immig_age25plus, native_age25plus)
educ_prop_plot_age25plus <- plot_educ_prop(
  plot_data_age25plus,
  title    = "Education Composition: Age 25+ at Arrival vs. U.S. Born (Age 25+)",
  subtitle = paste("Immigrant proportions by arrival year (YRIMMIG),",
                   "restricted to age at arrival >= 25"),
  filename = "outputs/figures/Educ_Prop_Age25plus.png",
  y_max    = 0.70
)

# --- (c) Age at arrival 18-24 (immigrants); natives shown with the same
#         AGE >= 25 reference line as (b) ---
immig_age18_24 <- reduced_data_iv[immigr == 1 & age_at_arrival >= 18 & age_at_arrival <= 24]

plot_data_age18_24      <- build_educ_prop_data(immig_age18_24, native_age25plus)
educ_prop_plot_age18_24 <- plot_educ_prop(
  plot_data_age18_24,
  title    = "Education Composition: Age 18-24 at Arrival vs. U.S. Born (Age 25+)",
  subtitle = paste("Immigrant proportions by arrival year (YRIMMIG),",
                   "restricted to age at arrival 18-24"),
  filename = "outputs/figures/Educ_Prop_Age18to24.png",
  y_max    = 0.70
)

# age_at_arrival is a scratch column on reduced_data_iv used only to build
# the two restricted figures above -- drop it so it doesn't linger on this
# large, long-lived object.
reduced_data_iv[, age_at_arrival := NULL]

# --- (d) Education composition by origin country -----------------------------
# The figure counterpart to Panel C of 22_origin_composition.R. Horizontal
# stacked bars, one per origin, SORTED BY COLLEGE-OR-MORE rather than by size:
# size order buries the finding, and the spread across origins is the point.
# A U.S.-born reference bar is drawn at the bottom so the origin figures have
# a scale.
#
# AGE >= 25, matching 21's desc_data and for the same reason -- Educ_Group is
# assigned from EDUCD with no age floor, so anyone still in school is bucketed
# as non_HS and the non-HS share is inflated in any table that does not add
# its own restriction.
#
# PERWT-WEIGHTED, consistent with figures (a)-(c) above, which all weight.
# NOTE: Panel C of 21 is currently UNWEIGHTED (it inherited `.N` counts from
# the old Section 22), so the two differ in basis. Weighted is the defensible
# choice for a population composition and is what this figure uses; 21 now
# reports both so the figure and the table can be reconciled.
cat("\nBuilding figure (d): education composition by origin...\n")

if (!exists("bpld_lookup"))
  stop("bpld_lookup not found. Run 08_load_analysis_subsets.R first -- it is\n",
       "built from the IPUMS DDI and is not in analysis_data.RData.")
bl <- as.data.table(bpld_lookup)

TOP_N_ORIGINS <- 20L
desc <- reduced_data_iv[AGE >= 25]
desc[, YEAR_n := as.integer(as.character(YEAR))]
# IPUMS columns arrive as haven_labelled. Harmless in comparisons and
# arithmetic, but rbind() refuses a column whose class differs across items and
# fifelse() refuses mismatched `yes`/`no`, so the U.S.-born reference row would
# fail to bind against the origin rows. Stripped once for the whole table.
.hv <- names(desc)[vapply(desc, inherits, logical(1), "haven_labelled")]
if (length(.hv))
  desc[, (.hv) := lapply(.SD, function(x) {
    v <- unclass(x); attributes(v) <- NULL; v
  }), .SDcols = .hv]
rm(.hv)

# Same ranking rule as 21: PERWT summed within survey year, then averaged
# across years so each survey year counts equally and the denser 2006-2023 ACS
# coverage does not swamp the single 2000 5% sample.
origin_rank <- desc[immigr == 1 & YEAR_n >= 2000 & YEAR_n <= 2023,
                    .(pop_weighted = sum(PERWT, na.rm = TRUE)),
                    by = .(BPLD, YEAR_n)
][, .(avg_pop_weighted = mean(pop_weighted)), by = BPLD
][order(-avg_pop_weighted)]
top_origins <- head(merge(origin_rank, bl, by = "BPLD",
                          all.x = TRUE)[order(-avg_pop_weighted)],
                    TOP_N_ORIGINS)

org <- desc[immigr == 1 & BPLD %in% top_origins$BPLD & !is.na(Educ_Group),
            .(w = sum(PERWT, na.rm = TRUE)), by = .(BPLD, Educ_Group)]
org[, share := w / sum(w), by = BPLD]
org <- merge(org, top_origins[, .(BPLD, label)], by = "BPLD")

usb <- desc[immigr == 0 & !is.na(Educ_Group),
            .(w = sum(PERWT, na.rm = TRUE)), by = Educ_Group]
usb[, `:=`(share = w / sum(w), label = "U.S.-born", BPLD = NA_integer_)]

org_plot <- rbind(org, usb, fill = TRUE)
org_plot[, Educ_Group := factor(Educ_Group,
                                levels = c("non_HS", "HighSch", "College", "Grad"))]

# Order the bars by college-or-more. The U.S.-born bar is forced to the bottom
# so it reads as a baseline rather than as one more origin in the ranking.
ord <- org_plot[Educ_Group %in% c("College", "Grad"),
                .(cp = sum(share)), by = label][order(cp)]
ord <- rbind(ord[label == "U.S.-born"], ord[label != "U.S.-born"])
org_plot[, label := factor(label, levels = ord$label)]

p_org <- ggplot(org_plot, aes(x = label, y = share, fill = Educ_Group)) +
  geom_col(width = 0.78) +
  coord_flip() +
  # reverse = TRUE because the stack renders Grad leftmost while the factor
  # orders non_HS first, so an unreversed key reads backwards against the bars
  scale_fill_manual(values = educ_colors[c("non_HS", "HighSch",
                                           "College", "Grad")],
                    labels = EDUC_LABELS[c("non_HS", "HighSch",
                                           "College", "Grad")],
                    name = "Education",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(title = "Education composition of immigrants, by country of birth",
       subtitle = paste("Top", TOP_N_ORIGINS,
                        "origins by population-weighted count; age 25+;",
                        "ordered by share with college or more"),
       x = NULL, y = "Share of the origin's immigrant population",
       caption = paste0(
         "IPUMS BPLD, person-weighted (PERWT), restricted to age 25 and over ",
         "because Educ_Group is assigned from EDUCD with no age floor.\n",
         "The U.S.-born bar is a reference, not one of the origins. Origins ",
         "are ranked by PERWT summed within each survey year and then ",
         "averaged across\nyears, so each survey year counts equally.")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7, colour = "grey30"))

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/figures/Educ_By_Origin.png", p_org, width = 10, height = 7,
       dpi = 300)
ggsave("outputs/figures/Educ_By_Origin.pdf", p_org, width = 10, height = 7)
cat("outputs/figures/Educ_By_Origin.{png,pdf} saved.\n")

cp <- org_plot[Educ_Group %in% c("College", "Grad"),
               .(college_plus = 100 * sum(share)), by = label][order(-college_plus)]
cat("\n  College-or-more by origin (weighted), as drawn:\n")
print(cp, row.names = FALSE, digits = 3)
cat(sprintf("\n  Spread across the %d origins: %.1f%% to %.1f%%, against %.1f%%\n",
            TOP_N_ORIGINS, cp[label != "U.S.-born", min(college_plus)],
            cp[label != "U.S.-born", max(college_plus)],
            cp[label == "U.S.-born", college_plus]))
cat("  for the U.S.-born. If this does not match Panel C of 21, check the\n")
cat("  weighting basis first -- this figure is PERWT-weighted.\n")

rm(desc, org, org_plot, reduced_data_iv); gc()

cat("\nSection 15 complete: four education-composition figures saved --\n")
cat("  Educ_Prop.png            (unrestricted)\n")
cat("  Educ_Prop_Age25plus.png  (age at arrival >= 25)\n")
cat("  Educ_Prop_Age18to24.png  (age at arrival 18-24)\n")
cat("  Educ_By_Origin.{png,pdf} (composition by country of birth)\n")

# No checkpoint written. after_sec15.RData existed only to pass
# all_natives_175 / bpld_lookup / top175_list_named to 23_msa_stability.R,
# which has taken them from analysis_data.RData directly since 2026-09-30.
cat("27_education_figures.R complete (no checkpoint -- see outputs/figures/).\n")
