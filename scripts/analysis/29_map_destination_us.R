# =============================================================================
# 29_map_destination_us.R  --  US maps: when each metro was the strongest
# Renumbered 2026-10-21 from 41_map_destination_us.R (archive-and-renumber pass).
#                              destination, by origin region
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# PURPOSE
# Four maps, one per origin region, each metro coloured by the five-year wave
# in which it drew its largest share of that region's US arrivals. Where the
# world maps in 40 show origin composition shifting, these show destination
# patterns shifting -- the other half of what makes a shift-share instrument
# do work.
#
# WHY "SHARE OF THE REGION'S ARRIVALS" AND NOT RAW COUNTS
# A dominant-wave map built on counts would mostly restate national flow
# volume: if one wave had the largest inflow, nearly every metro would show
# that wave as dominant and the map would come out close to monochrome. So for
# each metro and wave the measure is
#
#     share_{m,w} = arrivals to metro m in wave w / all arrivals in wave w
#
# within the region, and the dominant wave is the w where share_{m,w} peaks.
# That answers "when did this metro matter most as a destination," which is the
# timing question, rather than "when were flows largest," which is a property
# of the national series.
#
# WHAT THIS DISCARDS
# Magnitude. A metro receiving 400,000 arrivals and one receiving 400 are
# coloured identically if their timing matches. The MIN_ARRIVALS floor below
# keeps the smallest cells out, and the script reports how decisively each
# metro's dominant wave wins, so weak cases are visible. If the map reads as
# noise, switch to small multiples -- four panels per region, wave by panel,
# showing shares on a common scale.
#
# REGION DEFINITIONS (revised 2026-10-01 to five purely geographic groups):
#   Americas             Canada through South America
#   Asia                 East, Southeast and South Asia
#   Europe / Australia   Europe plus Australia, New Zealand and the Pacific
#   Africa               all of Africa, including South Africa
#   Middle East          its own region rather than folded into Africa or Asia
#
# The earlier version placed South Africa with Europe/Australia and the Middle
# East with Africa. Both were defensible groupings but neither was geographic,
# and each needed a caption sentence to justify it. Separating the Middle East
# and returning South Africa to Africa removes the need for both.
#
# GEOGRAPHY: usmap, county-based (switched from tigris/CBSA 2026-10-01)
# A metro is drawn as the union of its constituent counties, joined through
# data/crosswalks/crosswalk_state_countyfip_msa.csv. Outside New England CBSAs
# are defined as whole counties, so the shapes are close to the CBSA polygons
# the earlier version used. usmap positions Alaska and Hawaii natively, which
# removes the shift_geometry() handling the tigris version needed. State
# borders are overlaid on top of the counties.
#
# The trade-off: Puerto Rico is absent from usmap's county map, so its metros
# (San Juan among them) cannot be shaded. Because usmap's STATE map does carry
# PR, leaving it in would draw an outline that is permanently empty for a
# reason the caption does not give -- so PR is dropped from the state layer and
# the maps show the 50 states only. Metros are also lost where a CBSA code has
# no row in the county crosswalk. Both losses, plus metros drawn only in part,
# are audited per region below against the geometry usmap will actually render.
#
# REQUIRES: usmap, ggplot2
# Writes outputs/figures/fig_destination_<region>.{png,pdf} and the CSV.
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(ggplot2)
suppressPackageStartupMessages(library(usmap))

WAVES <- list("2000-2004" = 2000:2004, "2005-2009" = 2005:2009,
              "2010-2014" = 2010:2014, "2015-2019" = 2015:2019)
MIN_ARRIVALS <- 50L     # per metro-region, pooled across waves
OUT_DIR <- "outputs/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
reduced_data_iv <- as.data.table(.tmp_env$reduced_data_iv)
rm(.tmp_env); gc()

if (!exists("bpld_lookup"))
  stop("bpld_lookup not found. Run 08_load_analysis_subsets.R first -- it is\n",
       "built from the IPUMS DDI and is not in analysis_data.RData.")
bl <- as.data.table(bpld_lookup)
setnames(bl, names(bl)[1:2], c("BPLD", "origin_label"), skip_absent = TRUE)

# --- Region assignment -------------------------------------------------------
# Continents come from the IPUMS BPL code ranges (BPL = BPLD %/% 100), which
# are stable. The Middle East is carved out by LABEL rather than by code,
# because that block's exact code range is easy to get wrong and a silent
# misassignment would be invisible in the output. Matching on names makes the
# rule auditable -- the full list of what matched is printed below.
MIDDLE_EAST <- c("Afghanistan", "Iran", "Iraq", "Israel", "Israel/Palestine",
                 "Palestine", "Jordan", "Kuwait", "Lebanon", "Saudi Arabia",
                 "Syria", "Turkey", "Yemen", "Yemen Arab Republic",
                 "Yemen, PDR", "United Arab Emirates", "Oman", "Qatar",
                 "Bahrain", "Cyprus", "Gaza Strip", "West Bank")

bl[, BPL := as.integer(BPLD) %/% 100L]
bl[, region := fcase(
  BPL >= 150L & BPL <= 399L, "Americas",
  BPL >= 400L & BPL <= 499L, "Europe / Australia",
  BPL >= 500L & BPL <= 549L, "Asia",
  BPL >= 600L & BPL <= 699L, "Africa",
  BPL >= 700L & BPL <= 799L, "Europe / Australia",
  default = NA_character_
)]
# Middle East carved out of the Asia range into its own region. South Africa is
# no longer moved -- it stays in Africa with the rest of the continent.
bl[origin_label %in% MIDDLE_EAST, region := "Middle East"]

cat("\n", strrep("-", 70), "\n", sep = "")
cat("REGION ASSIGNMENT AUDIT -- check this before trusting the maps\n\n")
print(bl[!is.na(region), .(codes = .N), by = region][order(region)],
      row.names = FALSE)
una <- bl[is.na(region)]
cat(sprintf("\n  unassigned codes: %d\n", nrow(una)))
if (nrow(una)) {
  cat("  (US states, territories, 'at sea' and unknown codes belong here;\n")
  cat("   a real country in this list is a bug in the BPL ranges above)\n")
  print(head(una[, .(BPLD, origin_label, BPL)], 20), row.names = FALSE)
}
cat("\n  Middle East codes matched by label:\n")
print(bl[origin_label %in% MIDDLE_EAST, .(BPLD, origin_label)], row.names = FALSE)

# --- Arrival flows by metro, region and wave --------------------------------
arr <- reduced_data_iv[immigr == 1 & YRIMMIG > 0 & MSA != 0 &
                         YRIMMIG >= 2000L & YRIMMIG <= 2019L,
                       .(MSA, BPLD, YRIMMIG)]
arr <- merge(arr, bl[, .(BPLD, region)], by = "BPLD")
arr <- arr[!is.na(region)]
arr[, wave := NA_character_]
for (w in names(WAVES)) arr[YRIMMIG %in% WAVES[[w]], wave := w]
arr <- arr[!is.na(wave)]

cell <- arr[, .(arrivals = .N), by = .(region, wave, MSA)]
# Share of the region's national arrivals in that wave
cell[, wave_total := sum(arrivals), by = .(region, wave)]
cell[, share := arrivals / wave_total]

# --- Dominant wave per metro-region -----------------------------------------
tot <- cell[, .(total = sum(arrivals)), by = .(region, MSA)]
cell <- merge(cell, tot, by = c("region", "MSA"))
cell <- cell[total >= MIN_ARRIVALS]

dom <- cell[order(region, MSA, -share),
            .(dominant_wave = wave[1],
              top_share     = share[1],
              runner_share  = if (.N > 1) share[2] else NA_real_,
              total         = total[1]),
            by = .(region, MSA)]
dom[, decisiveness := top_share / runner_share]
dom[, single_wave := is.na(runner_share)]
dom[, dominant_wave := factor(dominant_wave, levels = names(WAVES))]

cat("\n", strrep("-", 70), "\n", sep = "")
cat("DOMINANT-WAVE SUMMARY\n\n")
print(dcast(dom, region ~ dominant_wave, value.var = "MSA",
            fun.aggregate = length), row.names = FALSE)

# Three distinct reasons a metro can go unshaded, which the map cannot
# distinguish and the reader should not have to guess at:
#   (a) no arrivals from the region at all -- absent from `tot` entirely
#   (b) arrivals below MIN_ARRIVALS -- in `tot` but filtered out
#   (c) shaded, but its dominant wave barely beats the runner-up
# Reported separately so the figure's silence is accounted for.
cat(sprintf("\n  metro-region cells with any arrivals : %d\n", nrow(tot)))
cat(sprintf("  of those, shaded (>= %d arrivals)      : %d\n",
            MIN_ARRIVALS, nrow(dom)))
cat(sprintf("  of those, filtered out as too thin    : %d\n",
            nrow(tot) - nrow(dom)))
cat("  Metros with NO arrivals from a region are absent from the count above\n")
cat("  and are also unshaded, so unshaded carries both meanings on the map.\n")

cat("\n  Decisiveness = top wave share / runner-up share. Values near 1 mean\n")
cat("  the colour is close to arbitrary for that metro. Metros with arrivals\n")
cat("  in only ONE wave are trivially dominant there and are counted\n")
cat("  separately -- they are the weakest evidence on the map, not the\n")
cat("  strongest, despite taking a solid colour.\n")
print(dom[, .(median_decisiveness = round(median(decisiveness, na.rm = TRUE), 2),
              pct_weak       = round(100 * mean(decisiveness < 1.2, na.rm = TRUE), 1),
              n_single_wave  = sum(single_wave)),
          by = region][order(region)], row.names = FALSE)

# --- Geography: MSA -> county FIPS -------------------------------------------
# The crosswalk mixes vintages: it carries both the old 4-digit MSA/PMSA codes
# and the 5-digit 2013 CBSA codes the analysis uses. Filtering to >= 10000
# keeps the CBSA-era rows, which is what MET2013 reports. 929 CBSA codes are
# present, New York resolving to 25 counties and Chicago to 14.
cw <- as.data.table(read.csv("data/crosswalks/crosswalk_state_countyfip_msa.csv"))
cw <- cw[MSA >= 10000L]
cw[, fips := sprintf("%02d%03d", STATEFIP, COUNTYFIP)]
cat(sprintf("\nCounty crosswalk: %d CBSA codes, %d county rows\n",
            uniqueN(cw$MSA), nrow(cw)))

# --- What usmap can actually draw --------------------------------------------
# usmap's STATE map includes Puerto Rico; its COUNTY map does not. The two
# disagree, and the disagreement is silent: a PR metro carries a valid 5-digit
# CBSA code, so it survives the crosswalk join above and is only dropped at
# render time, leaving an empty PR outline that the caption cannot explain.
#
# Two corrections follow. The state layer excludes PR so nothing is drawn that
# cannot be filled, and the audit below compares shaded county FIPS against the
# FIPS usmap will actually render, so any future mismatch of this kind is
# reported rather than absorbed into the white space.
drawable_fips <- {
  .cm <- usmap::us_map(regions = "counties")
  .f  <- unique(as.character(.cm$fips))
  # Puerto Rico is excluded from these maps by design, and is dropped from
  # both drawn layers below. It is dropped here as well so "drawable" means
  # what the audit reports it to mean: a PR metro is LOST AT RENDER, which is
  # true and intended, rather than counted as shaded and then silently absent.
  .f[substr(.f, 1, 2) != "72"]
}

.sl <- usmap::us_map(regions = "states")
.sl <- .sl[substr(as.character(.sl$fips), 1, 2) != "72", ]   # drop Puerto Rico
state_layer <- if (inherits(.sl, "sf")) {
  geom_sf(data = .sl, fill = NA, colour = "grey35", linewidth = 0.28)
} else {
  geom_polygon(data = .sl, aes(x = x, y = y, group = group),
               fill = NA, colour = "grey35", linewidth = 0.28)
}
cat(sprintf("State layer: %s | %d states drawn (PR excluded)\n",
            if (inherits(.sl, "sf")) "sf (usmap >= 0.7)" else "data.frame (legacy usmap)",
            uniqueN(.sl$fips)))
cat(sprintf("usmap county map carries %d county FIPS\n", length(drawable_fips)))

# --- Connecticut: legacy counties -> planning regions ------------------------
# Connecticut's eight counties stopped functioning as governmental units in
# 1960, and at the State's request the Census Bureau replaced them with the
# nine council-of-government planning regions as Connecticut's county
# equivalent (87 FR 34235, 6 June 2022). The legacy codes 09001-09015 were
# retired; the planning regions carry 09110-09190, and the two sets are
# disjoint. The MSA crosswalk predates the change, so every Connecticut metro
# joins on a retired code, survives the join carrying a valid CBSA, and then
# disappears at render against modern usmap geometry -- silently, which is
# what the render audit below was added to catch.
#
# The geographies do not nest: planning regions are built from towns and cut
# across the old county lines, so no exact correspondence exists and none is
# invented here. The Census Bureau publishes its own, in the same notice
# (Table 2, "Counties-to-Planning Regions Approximation"), and that published
# table is what is applied below. At the level this map draws, the
# approximation does no work: each of the nine regions falls inside exactly
# one Connecticut metro, so no region is contested and no metro's dominant
# wave depends on the correspondence. What shifts is where a boundary sits,
# not which metro is shaded or what colour it takes. The caption says so.
CT_LEGACY <- sprintf("09%03d", c(1L, 3L, 5L, 7L, 9L, 11L, 13L, 15L))
CT_XWALK <- data.table(
  fips_old = c("09001", "09001", "09003", "09013", "09007",
               "09009", "09009", "09011", "09015", "09005"),
  fips_new = c("09120", "09190", "09110", "09110", "09130",
               "09140", "09170", "09180", "09150", "09160"))

ct_rows <- cw[fips %in% CT_LEGACY]
if (nrow(ct_rows) && !any(CT_LEGACY %in% drawable_fips) &&
    all(CT_XWALK$fips_new %in% drawable_fips)) {

  ct_new <- merge(ct_rows, CT_XWALK, by.x = "fips", by.y = "fips_old",
                  allow.cartesian = TRUE)
  ct_new <- unique(ct_new[, .(MSA, fips = fips_new)])
  ct_new[, `:=`(STATEFIP = 9L, COUNTYFIP = NA_integer_)]
  setcolorder(ct_new, names(cw))
  cw <- rbind(cw[!fips %in% CT_LEGACY], ct_new)

  cat(sprintf(paste0("\nConnecticut: %d legacy county row(s) -> %d planning-",
                     "region row(s), %d metro(s): %s\n  Source: 87 FR 34235 ",
                     "Table 2 (Census Bureau's published approximation).\n"),
              nrow(ct_rows), nrow(ct_new), uniqueN(ct_new$MSA),
              paste(sort(unique(ct_new$MSA)), collapse = ", ")))

} else if (nrow(ct_rows) && any(CT_LEGACY %in% drawable_fips)) {
  cat("\nConnecticut: usmap still carries legacy CT counties; no remap needed.\n")
} else if (nrow(ct_rows)) {
  cat("\nConnecticut: usmap carries neither the legacy counties nor the nine\n",
      "  expected planning regions. CT will not render -- check usmap's CT\n",
      "  geometry before trusting these figures.\n", sep = "")
}

# How much of the crosswalk can usmap render at all? This is a property of the
# two sources, not of any one region, so it is reported once.
cw_undrawable <- cw[!fips %in% drawable_fips]
if (nrow(cw_undrawable)) {
  cat(sprintf("\n  %d of %d crosswalk counties (%d metros) are absent from usmap's\n",
              nrow(cw_undrawable), nrow(cw), uniqueN(cw_undrawable$MSA)))
  cat("  county map and can never be shaded. By state FIPS:\n")
  print(cw_undrawable[, .(counties = .N, metros = uniqueN(MSA)),
                      by = .(state_fips = substr(fips, 1, 2))][order(-counties)],
        row.names = FALSE)
}

# --- One map per region ------------------------------------------------------
WAVE_COLS <- c("2000-2004" = "#1b6ca8", "2005-2009" = "#2e9e5b",
               "2010-2014" = "#e8a33d", "2015-2019" = "#c7472f")
COUNTY_BORDER <- "grey75"   # NA suppresses county lines inside metros

# -----------------------------------------------------------------------------
# LEGEND AND CONTRAST SETTINGS -- added 2026-10-19, three changes, all revertible
# from this block alone. To return to the first gradient version:
#
#   (1) NA key in the hue legend. Was a visible "NA" swatch below 2015-2019,
#       reading as a fifth wave. Now absent because the fill scale is an
#       identity scale over explicit breaks. REVERT: restore
#       scale_fill_manual(..., na.value = UNSHADED) + scale_alpha_continuous().
#
#   (2) Contrast floor. ALPHA_FLOOR was 0.35 and UNSHADED was "grey96", which
#       put the palest metros close enough to the unshaded counties to be
#       misread as not-in-sample. REVERT: set ALPHA_FLOOR back to 0.35 and
#       UNSHADED back to "grey96".
#
#   (3) Bivariate legend. Hue and opacity were two separate legends, the
#       second a grey ramp that never showed what a pale blue meant. Now one
#       4x3 grid of the actual map colours: rows are waves, columns run pale
#       to dark. REVERT: drop KEY_* and the scale_fill_identity() call.
#
# The map itself still varies continuously; the legend samples that ramp at
# three points. Nothing about which metro takes which hue has changed in any
# of the three.
# -----------------------------------------------------------------------------
ALPHA_FLOOR <- 0.50          # was 0.35
UNSHADED    <- "grey88"      # was "grey96"

# NOT BINNED. Two intermediate versions were tried and rejected:
#
#   scale_fill_identity() over 12 blended hexes, to get one bivariate legend
#   while keeping a continuous map. This fails silently -- an identity scale
#   only draws breaks that occur in the data, and a continuous ramp produces
#   hundreds of distinct fills, so eleven of the twelve keys matched nothing
#   and the legend came out as a scattered half-grid.
#
#   Terciles of intensity with a 12-level factor scale, which fixed the
#   legend but quantised the map into three steps per hue.
#
# The continuous ramp is kept because the gradient is the point: it reads as
# a wave washing across the metros, which three discrete steps do not. The
# cost is two legends rather than one, since ggplot cannot build a single key
# from two scales without quantising one of them.

slug <- function(x) tolower(gsub("[^A-Za-z]+", "_", x))

for (rg in sort(unique(dom$region))) {
  d  <- dom[region == rg]
  pd <- merge(cw[, .(MSA, fips)],
              d[, .(MSA, dominant_wave, total)], by = "MSA", allow.cartesian = TRUE)

  # Volume as a second channel, added 2026-10-19. Colour alone said only WHEN
  # a metro peaked, and the caption had to concede that a metro receiving many
  # arrivals and one receiving few took the same colour. Opacity now carries
  # the count, so timing and magnitude are both legible and the concession is
  # no longer needed.
  #
  # The scale is log10 because the distribution is not remotely even: within a
  # region the largest metro runs 65 to 205 times the median, so on a linear
  # ramp every metro but the top two or three would sit at the floor. It is
  # normalised WITHIN region -- Africa's largest destination draws 16,768
  # arrivals against the Americas' 132,170 -- so each map uses its own range
  # and opacity is comparable across metros within a map, never between maps.
  # The floor of 0.35 keeps the smallest metros visible as colour rather than
  # fading them into the unshaded counties, which would confuse "few arrivals"
  # with "none".
  pd[, lg := log10(pmax(total, 1))]
  .rng <- range(pd$lg, finite = TRUE)
  pd[, intensity := if (diff(.rng) > 0) (lg - .rng[1]) / diff(.rng) else 1]

  # Two distinct ways a metro disappears, reported separately because they
  # point at different sources. A metro LOST AT THE JOIN has no row in the
  # county crosswalk. A metro LOST AT RENDER has crosswalk rows, but usmap
  # holds no geometry for any of them -- the PR case. Only the first was
  # visible before; the second was silently absorbed into the white space.
  lost_join <- setdiff(d$MSA, cw$MSA)
  pd <- pd[, drawable := fips %in% drawable_fips][]
  lost_render <- pd[, .(any_drawn = any(drawable)), by = MSA][any_drawn == FALSE, MSA]
  partial     <- pd[, .(all_drawn = all(drawable)), by = MSA][all_drawn == FALSE, MSA]
  partial     <- setdiff(partial, lost_render)
  pd          <- pd[drawable == TRUE]

  cat(sprintf("\n  %s: %d metros shaded, %d counties",
              rg, uniqueN(pd$MSA), nrow(pd)))
  if (length(lost_join))
    cat(sprintf("\n    %d metro(s) LOST AT THE JOIN (no crosswalk row): %s",
                length(lost_join), paste(head(lost_join, 6), collapse = ", ")))
  if (length(lost_render))
    cat(sprintf("\n    %d metro(s) LOST AT RENDER (no usmap geometry): %s",
                length(lost_render), paste(head(lost_render, 8), collapse = ", ")))
  if (length(partial))
    cat(sprintf("\n    %d metro(s) PARTIALLY drawn (some counties missing): %s",
                length(partial), paste(head(partial, 6), collapse = ", ")))
  cat("\n")

  # Built from the county sf directly rather than through plot_usmap(), which
  # maps one column to fill and gives no way to pass a second aesthetic. The
  # state layer was already a geom_sf, so this makes both layers consistent.
  .cty <- usmap::us_map(regions = "counties")
  if (!inherits(.cty, "sf"))
    stop("usmap returned a legacy data.frame for counties; the opacity channel ",
         "needs the sf geometry of usmap >= 0.7.")
  .cty$fips <- as.character(.cty$fips)

  # Drop Puerto Rico, matching the state layer. The original code only needed
  # this filter on the state layer, because usmap's county map carried no PR
  # geometry while its state map did -- that asymmetry is what the render
  # audit below was written to catch. Recent usmap ships PR counties too, so
  # with the county layer now drawn directly the filter is needed in both
  # places. Applied on the geometry rather than on a usmap version test, so it
  # holds whichever way the package goes.
  .cty <- .cty[substr(.cty$fips, 1, 2) != "72", ]

  .cty <- merge(.cty, as.data.frame(pd[, .(fips, dominant_wave, intensity)]),
                by = "fips", all.x = TRUE)

  p <- ggplot() +
    geom_sf(data = .cty, aes(fill = dominant_wave, alpha = intensity),
            color = COUNTY_BORDER, linewidth = 0.05) +
    state_layer +
    # na.translate = FALSE drops the "NA" key that otherwise appeared below
    # 2015-2019 and read as a fifth wave; unshaded counties still render in
    # UNSHADED via na.value, and the caption accounts for them in words.
    scale_fill_manual(values = WAVE_COLS, drop = FALSE,
                      na.value = UNSHADED, na.translate = FALSE,
                      name = "Wave of peak\ndestination share",
                      guide = guide_legend(order = 1,
                                           override.aes = list(alpha = 1))) +
    # The grey ramp is shown at the same floor the map uses, so the lightest
    # key matches the palest metro rather than implying the scale reaches 0.
    scale_alpha_continuous(range = c(ALPHA_FLOOR, 1), limits = c(0, 1),
                           breaks = c(0, 0.5, 1),
                           labels = c("fewest", "", "most"),
                           na.value = 1, name = "Arrivals\n(log scale)",
                           guide = guide_legend(order = 2,
                                                override.aes = list(
                                                  fill = "grey25", colour = NA))) +
    labs(
      title = sprintf("When each metro drew its largest share of %s arrivals", rg),
      subtitle = paste("Hue marks the five-year wave in which the metro's share of",
                       "the region's US arrivals peaked;\nopacity marks how many",
                       "arrivals it received"),
      caption = paste0(
        "Arrivals identified from year of immigration (YRIMMIG). For each metro ",
        "and wave, the measure is the metro's share of all arrivals from this ",
        "region in that wave,\nso hue reflects timing. Opacity is the metro's ",
        "total arrivals from the region over 2000-2019 on a log scale, ",
        "normalised within this map:\nthe largest metro in a region runs 65 to ",
        "205 times the median, so a linear ramp would leave all but a handful ",
        "at the floor. Opacity is\ntherefore comparable between metros on this ",
        "map and not between maps. Unshaded counties are either ",
        "non-metropolitan or belong to a metro\nthat recorded no arrivals from ",
        "the region, or fewer than ",
        MIN_ARRIVALS, " across 2000-2019. Metros are drawn as their ",
        "constituent counties.\nRegion definitions are stated in the text. ",
        "Connecticut is ",
        "drawn on the nine planning regions that replaced its counties as ",
        "county equivalents in 2022, following the Census Bureau's published ",
        "correspondence\n(87 FR 34235, Table 2). The two geographies do not ",
        "nest, so Connecticut's metro boundaries are approximate; which metro ",
        "takes which colour is not affected.")
    ) +
    # plot_usmap() supplied the blank map theme; building the layers directly
    # means supplying it here. theme_void() drops the panel, grid and axes
    # that would otherwise frame the map in grey with degree ticks.
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(hjust = 0),
          plot.caption = element_text(hjust = 0, size = 7, colour = "grey30"),
          plot.margin = margin(8, 8, 8, 8))

  f <- file.path(OUT_DIR, sprintf("fig_destination_%s", slug(rg)))
  ggsave(paste0(f, ".png"), p, width = 10, height = 6.5, dpi = 300)
  ggsave(paste0(f, ".pdf"), p, width = 10, height = 6.5)
  cat(sprintf("  saved: %s.{png,pdf}  (%d metros shaded)\n", f, uniqueN(pd$MSA)))
}

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(dom, "outputs/tables/destination_dominant_wave.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/destination_dominant_wave.csv\n")

cat("\n", strrep("-", 70), "\n", sep = "")
cat("IF THESE READ AS NOISE\n\n")
cat("Check pct_weak and n_single_wave together. A high pct_weak means the top\n")
cat("wave barely beats the runner-up for much of the map, so the colouring is\n")
cat("close to arbitrary. A high n_single_wave means many metros are coloured\n")
cat("on the basis of appearing in one wave only -- solid colour, thinnest\n")
cat("possible evidence. Either would mean the dominant-wave design is not\n")
cat("carrying information.\n\n")
cat("The fallback is small multiples: four panels per region, one per wave,\n")
cat("each showing share on a common scale. That keeps magnitude and loses the\n")
cat("single-glance timing read.\n")
cat("29_map_destination_us.R complete.\n")
