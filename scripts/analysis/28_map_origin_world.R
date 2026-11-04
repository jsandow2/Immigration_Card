# =============================================================================
# 28_map_origin_world.R  --  World maps: origin composition of US immigrant
# Renumbered 2026-10-21 from 40_map_origin_world.R (archive-and-renumber pass).
#                            arrivals, 1990 / 2005 / 2019
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# PURPOSE
# The shift-share instrument weights national arrival flows by each metro's
# 1990 origin composition. Whether that is a credible instrument depends on
# origin composition shifting over time for reasons external to any one metro,
# and at present that shift exists in the chapter only as a table of Rotemberg
# weights. These three maps show it directly.
#
# WHAT IS PLOTTED
# FLOW, not stock: arrivals in the given year, identified by YRIMMIG, rather
# than the foreign-born population resident that year. Flow is what
# national_inflow measures in 06_shift_share_iv.R, so this is the object the
# instrument is built from.
#
# SHARES, not counts. A share answers the question the figure is asking --
# what fraction of a year's arrivals came from each country -- and it is the
# quantity the shift-share is built on, since the instrument weights national
# flows by share_{m,k,1990}. It is also scale-free, so all three panels read
# against one legend and the composition shift is visible directly rather than
# having to be inferred from differing colour ranges.
#
# A SHARED COLOUR SCALE across all three panels. Without it each map's colours
# mean something different and the comparison the figure exists to support is
# not available to the reader.
#
#
# REQUIRES: sf, rnaturalearth, rnaturalearthdata, ggplot2, patchwork
# Writes outputs/figures/fig_origin_world.{png,pdf} and the underlying CSV.
# =============================================================================

source("scripts/analysis/00_setup_analysis.R")

library(data.table)
library(ggplot2)
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
  library(RColorBrewer)
})

MAP_YEARS <- c(1990L, 2005L, 2019L)   # 2019 avoids the 2020-21 COVID collapse
OUT_DIR   <- "outputs/figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading analysis_data.RData...\n")
.tmp_env <- new.env()
load("analysis_data.RData", envir = .tmp_env)
reduced_data_iv <- as.data.table(.tmp_env$reduced_data_iv)
rm(.tmp_env); gc()

# --- Arrival flows by origin and year ---------------------------------------
# reduced_data_iv retains 100% of foreign-born (04_sample_construction.R), so
# no weighting correction is needed on the immigrant side.
flows <- reduced_data_iv[immigr == 1 & YRIMMIG > 0 & YRIMMIG %in% MAP_YEARS,
                         .(arrivals = .N), by = .(YRIMMIG, BPLD)]
flows[, share := arrivals / sum(arrivals), by = YRIMMIG]
setorder(flows, YRIMMIG, -share)

cat("\nArrivals by panel year:\n")
print(flows[, .(origins = .N, arrivals = sum(arrivals)), by = YRIMMIG],
      row.names = FALSE)

# --- Origin labels from the IPUMS codebook ----------------------------------
# bpld_lookup is built in 08_load_analysis_subsets.R from ipums_val_labels().
# Rebuilt here if absent so this script stands alone.
if (!exists("bpld_lookup")) {
  stop("bpld_lookup not found. Run 08_load_analysis_subsets.R first, or load\n",
       "a checkpoint containing it -- it is built from the IPUMS DDI and\n",
       "cannot be reconstructed from analysis_data.RData alone.")
}
bl <- as.data.table(bpld_lookup)
setnames(bl, names(bl)[1:2], c("BPLD", "origin_label"), skip_absent = TRUE)
flows <- merge(flows, bl[, .(BPLD, origin_label)], by = "BPLD", all.x = TRUE)

# --- BPLD -> country join for mapping ---------------------------------------
# The join is by NAME against Natural Earth, because IPUMS BPLD codes do not
# carry ISO identifiers. Most match directly; the overrides below cover the
# cases that cannot. Aggregated codes ("Other USSR/Russia") are assigned to
# their largest constituent, which is a simplification the caption should note.
#
# AUDIT THIS. The script reports every unmatched origin with its arrival share,
# so the cost of any failure is visible rather than silent. If a large origin
# appears in that list, add it here.
ne_overrides <- c(
  "Other USSR/Russia"            = "Russia",
  "USSR, ns"                     = "Russia",
  "Russia"                       = "Russia",
  "Guyana/British Guiana"        = "Guyana",
  "Korea"                        = "South Korea",
  "South Korea"                  = "South Korea",
  "China"                        = "China",
  "Hong Kong"                    = "Hong Kong S.A.R.",
  "Taiwan"                       = "Taiwan",
  "Burma (Myanmar)"              = "Myanmar",
  "Myanmar (Burma)"              = "Myanmar",
  "Czechoslovakia"               = "Czechia",
  "Yugoslavia"                   = "Serbia",
  "Bosnia and Herzegovina"       = "Bosnia and Herz.",
  "Macedonia"                    = "North Macedonia",
  "Dominican Republic"           = "Dominican Rep.",
  "Central African Republic"     = "Central African Rep.",
  "Congo"                        = "Dem. Rep. Congo",
  "Ivory Coast"                  = "Côte d'Ivoire",
  "Cape Verde"                   = "Cabo Verde",
  "Laos"                         = "Lao PDR",
  "Trinidad and Tobago"          = "Trinidad and Tobago",
  "England"                      = "United Kingdom",
  "Scotland"                     = "United Kingdom",
  "Wales"                        = "United Kingdom",
  "Northern Ireland"             = "United Kingdom",
  "United Kingdom, ns"           = "United Kingdom",
  "Palestine"                    = "Palestine",
  "Syria"                        = "Syria",
  "South Sudan"                  = "S. Sudan",
  "Swaziland"                    = "eSwatini",
  "Equatorial Guinea"            = "Eq. Guinea",
  "Solomon Islands"              = "Solomon Is.",
  "Marshall Islands"             = "Marshall Is."
)

flows[, map_name := fifelse(origin_label %in% names(ne_overrides),
                            ne_overrides[origin_label], origin_label)]

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world[world$admin != "Antarctica", ]

# --- Match report ------------------------------------------------------------
cat("\n", strrep("-", 70), "\n", sep = "")
cat("ORIGIN-TO-MAP MATCH AUDIT\n")
matched <- unique(flows$map_name) %in% world$admin
cat(sprintf("  distinct origins: %d | matched to a country polygon: %d\n",
            length(matched), sum(matched)))
unm <- flows[!(map_name %in% world$admin),
             .(arrivals = sum(arrivals), share = sum(share) / length(MAP_YEARS)),
             by = .(origin_label, map_name)][order(-arrivals)]
if (nrow(unm)) {
  cat(sprintf("  UNMATCHED: %d origins, %.2f%% of pooled arrivals\n",
              nrow(unm), 100 * sum(unm$share)))
  cat("  Largest unmatched (add to ne_overrides if any of these matter):\n")
  print(head(unm, 15), row.names = FALSE)
} else {
  cat("  all origins matched\n")
}
cat("  Note: aggregates with no single country (e.g. 'Americas, ns') cannot\n")
cat("  match by construction and are expected here.\n")

# --- Build the plotting frame ------------------------------------------------
# A COMPLETE country x year grid, not a left join on the flow data.
#
# The earlier version merged the flow table straight onto the country polygons.
# Because that table holds one row per country-year, the merge duplicated
# geometries per year and produced three distinct failure modes: a country with
# 1990 arrivals but none in 2019 appeared only in the 1990 facet and left a
# HOLE in 2019; countries with no arrivals in any year carried YRIMMIG = NA and
# generated a spurious fourth facet; and both rendered in na.value grey, which
# was visually indistinguishable from "this origin could not be mapped."
#
# Building the full grid first and zero-filling fixes all three. Every country
# now appears in every panel, and zero arrivals is a value rather than an
# absence.
plot_dt <- flows[, .(share = sum(share)), by = .(YRIMMIG, map_name)]
grid <- CJ(admin = unique(world$admin), YRIMMIG = MAP_YEARS)
grid <- merge(grid, plot_dt,
              by.x = c("admin", "YRIMMIG"), by.y = c("map_name", "YRIMMIG"),
              all.x = TRUE)
grid[is.na(share), share := 0]

# Shared breaks across panels. Log-spaced because the distribution is extremely
# skewed -- Mexico alone is a large multiple of the median origin -- and a
# linear scale renders every panel as one dark country on a pale world.
#
# "none" is its own level rather than being folded into the smallest bin, so a
# country that sent nobody is distinguishable from one that sent a handful.
brks <- c(0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 1)
labs <- c("<0.1%", "0.1-0.5%", "0.5-1%", "1-2.5%", "2.5-5%", "5-10%", ">10%")
grid[, share_bin := as.character(cut(share, breaks = brks, labels = labs,
                                     include.lowest = TRUE))]
grid[share == 0, share_bin := "none"]
grid[, share_bin := factor(share_bin, levels = c("none", labs))]

map_sf <- merge(world[, c("admin", "geometry")], as.data.frame(grid),
                by = "admin")
map_sf <- st_transform(map_sf, crs = "+proj=robin")

# --- Share of arrivals that cannot be placed on the map ---------------------
# Reported per panel in the facet label, so the figure states how much of each
# year's flow is missing from it. Aggregated codes with no single country
# ("Americas, ns", "Asia, ns") are the bulk of this and cannot be mapped by
# construction; dropping them silently would overstate what the figure shows.
unmapped <- flows[!(map_name %in% world$admin),
                  .(unmapped = sum(share)), by = YRIMMIG]
unmapped <- merge(data.table(YRIMMIG = MAP_YEARS), unmapped,
                  by = "YRIMMIG", all.x = TRUE)
unmapped[is.na(unmapped), unmapped := 0]
cat("\nShare of arrivals not placeable on the map, by panel:\n")
print(unmapped[, .(YRIMMIG, pct_unmapped = round(100 * unmapped, 2))],
      row.names = FALSE)

# Figures for the caption. Built from the computed values rather than typed in,
# so the sentence cannot drift away from the data if anything upstream changes.
setorder(unmapped, YRIMMIG)
unmapped_phrase <- paste0(
  paste(sprintf("%.1f%%", 100 * unmapped$unmapped[-nrow(unmapped)]),
        collapse = ", "),
  " and ", sprintf("%.1f%%", 100 * unmapped$unmapped[nrow(unmapped)]))

fill_cols <- setNames(
  c("grey88", RColorBrewer::brewer.pal(7, "YlOrRd")),
  c("none", labs))

p <- ggplot(map_sf) +
  geom_sf(aes(fill = share_bin), colour = "grey60", linewidth = 0.08) +
  facet_wrap(~ YRIMMIG, ncol = 1) +
  scale_fill_manual(values = fill_cols, na.value = "white", drop = FALSE,
                    name = "Share of US\narrivals that year") +
  labs(
    title = "Origin composition of US immigrant arrivals",
    subtitle = "Share of all arrivals in the year, by country of birth",
    caption = paste0(
      "Arrivals identified from year of immigration (YRIMMIG), pooled across ",
      "survey waves. Each country is shaded by its share of that year's total ",
      "arrivals.\n'None' means no recorded arrivals from that country in that ",
      "year. Origin codes that aggregate several countries ('Americas' and ",
      "similar) cannot\nbe placed on a map; the share of each year's arrivals ",
      "they account for are ", unmapped_phrase, " respectively.")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_line(colour = "grey95", linewidth = 0.2),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    plot.caption = element_text(hjust = 0, size = 7, colour = "grey30")
  )

ggsave(file.path(OUT_DIR, "fig_origin_world.png"), p,
       width = 9, height = 12, dpi = 300)
ggsave(file.path(OUT_DIR, "fig_origin_world.pdf"), p,
       width = 9, height = 12)
cat(sprintf("\nSaved: %s/fig_origin_world.{png,pdf}\n", OUT_DIR))

# --- Top origins per panel, for the text ------------------------------------
cat("\n", strrep("-", 70), "\n", sep = "")
cat("TOP 10 ORIGINS BY PANEL (for the accompanying text)\n")
for (y in MAP_YEARS) {
  cat(sprintf("\n  %d:\n", y))
  print(flows[YRIMMIG == y][order(-share)][1:10,
        .(origin_label, arrivals, pct = round(100 * share, 2))],
        row.names = FALSE)
}

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(flows[, .(YRIMMIG, BPLD, origin_label, map_name, arrivals, share)],
          "outputs/tables/origin_shares_by_year.csv", row.names = FALSE)
cat("\nSaved: outputs/tables/origin_shares_by_year.csv\n")
cat("28_map_origin_world.R complete.\n")
