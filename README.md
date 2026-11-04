# Immigration_Card

## Mother of Exiles, Maker of Scholars: Immigrant Exposure and U.S.-born Educational Attainment

**Author:** Joe Sandow  
**Institution:** University of Nebraska-Lincoln  
**Advisor:** Ed  
**Status:** Dissertation Chapter 2 — job market paper (in progress)

---

## Overview

This repository holds the R pipeline and LaTeX source for Chapter 2 of the
dissertation, the job market paper: a causal analysis of how exposure to
immigrant inflows during the formative years (ages 15--25) affects the
educational attainment of U.S.-born adults.

The data are individual-level records from IPUMS USA covering 1990--2023,
restricted to the 175 metropolitan areas that appear most consistently among
the largest 175 across the sample years. That yields an estimation sample of
8,537,597 U.S.-born observations across 45 states, of which 3,484,287 have
non-zero exposure and identify the coefficient.

Identification follows Card (2001). The instrument is the 1990 distribution
of each origin group's U.S. population across metropolitan areas, interacted
with that origin's national arrival counts, averaged over each individual's
ages 15--25 window. The endogenous variable is the realized inflow averaged
over the *same* window, so the first stage asks how well historical settlement
predicts the inflows a person actually experienced while deciding how much
schooling to acquire.

Because exposure is a deterministic function of metropolitan area and birth
year, treatment is defined at the **metropolitan-area by birth-cohort** level
and takes 5,775 distinct values across 175 metros and 33 cohorts. Standard
errors are clustered on metropolitan area, above the level of assignment, so
precision comes from variation across the 175 clusters rather than from the
number of individual records.

---

## Main Findings

Two-stage least squares on the high-school margin gives a coefficient of
1.12e-05 (MSA-clustered SE 1.50e-06), significant at the 1 percent level. The
interpretable magnitude is **4.74 percentage points per standard deviation** of
exposure. At the college margin the effect is 3.47 percentage points per
standard deviation. Adding Asian and Aboriginal identity indicators leaves the
high-school estimate unchanged at 4.74.

The first-stage coefficient is 0.647 (SE 0.105): historical settlement
patterns predict roughly two-thirds of realized inflow over the exposure
window.

**The sign is the result.** The estimates weigh against the proposition that
immigrant inflows displace U.S.-born students from educational opportunity.
They do not support a claim about any particular magnitude — restricting to
the identifying subsample roughly triples the estimate, to 14.33 and 10.33
points, not because the pooled specification was wrong but because the
condition creating non-zero exposure also makes that subsample young, and
younger respondents are more often recorded as not having completed a level of
schooling they have not yet finished.

Supporting results replicate Card (2001) on two outcomes using longer data: no
detectable effect of immigrant inflows on U.S.-born internal migration (no
estimate distinguishable from zero), and no negative wage effect at the
low-skill end.

Identification is checked four ways: a Rotemberg decomposition, where no single
origin carries more than 11.2 percent of the weight; the Adão--Kolesár--Morales
placebo; a leave-one-out specification; and a test of whether metropolitan
attribution error in the endogenous variable correlates with the instrument
within metro, which it does not.

---

## Repository Structure

```
Immigration_Card/
├── Immigration_Card.tex          # the chapter; compiles to ~82pp
├── bibliography.bib              # references
├── nuthesis.cls, aea.bst         # UNL dissertation class, AEA bib style
├── scripts/
│   ├── data_construction/        # 01-07: ingest, crosswalk, variables,
│   │   └── run_all.R             #        sample, stable MSAs, shift-share IV
│   ├── analysis/                 # 08-29: estimation, diagnostics, rendering
│   │   ├── 00_setup_analysis.R   #        libraries, paths, BPLD lookup
│   │   ├── 00b_table_helper.R    #        shared LaTeX table writer
│   │   ├── archive/              #        retired scripts, kept for provenance
│   │   └── run_all.R
│   └── README.md                 # per-script detail and the renumbering history
├── data/crosswalks/              # published PUMA / CBSA / county reference tables
├── Codebooks/                    # where your extract's DDI goes (not committed)
├── outputs/
│   ├── tables/                   # .tex tables the chapter \input, plus their source .csv
│   └── figures/                  # .pdf / .png figures the chapter includes
└── archive/                      # superseded scripts and drafts
```

The analysis scripts are numbered **08--29**, continuing the 01--07 sequence in
`data_construction/`, so the whole pipeline reads as one ordered run.

The monolithic originals the pipeline was split from --
`Immigration_Card_Data_Construction.R`, `Immigration_Card_Analysis.R` and
`Uploaded Code.R` -- are in `archive/monolithic_originals/`. Script headers
still say "Part of the modular Immigration_Card_Analysis.R pipeline"; that is
a provenance note, not a dependency. Nothing sources them.

---

## Pipeline

```r
source("scripts/data_construction/run_all.R")   # 01-07 -> analysis_data.RData
source("scripts/analysis/run_all.R")            # 08-29 -> tables and figures
```

| Scripts | Stage |
|---|---|
| `01`--`07` | IPUMS ingest, MSA crosswalk, variable and sample construction, stable-MSA selection, shift-share and exposure IV |
| `08` | load analysis subsets |
| `09`--`14` | estimation on the main design: headline IV, pooled-vs-identifying, birth cohort, age-FE check, two margins, by sex, and `14b`'s multinomial decomposition |
| `15`--`19` | identification: Rotemberg decomposition, Adão--Kolesár--Morales placebo, Asian-placebo follow-up, leave-one-out, metro attribution |
| `20`--`21` | Card (2001) replications: internal migration, occupation wages |
| `22`--`23` | descriptives and metropolitan panel composition |
| `24`--`29` | rendering: result tables, demographic summary, appendix tables, education figures, origin and destination maps |

Order matters in three places only, documented at the top of
`scripts/analysis/run_all.R`. Everything else loads `analysis_data.RData`
itself and can be run standalone -- `00_setup_analysis.R` builds and caches the
BPLD lookup, so no script depends on another having run in the same session.

Data construction is the expensive stage. The analysis half runs in minutes.

---

## Data Sources and Reproduction

Everything the pipeline needs is in this repository **except the IPUMS USA
microdata extract**, which IPUMS requires each user to register for and
generate themselves. It cannot be redistributed.

### Already in the repository

| Path | What it is |
|---|---|
| `data/crosswalks/` | Published PUMA, CBSA and county reference tables, 5.2 MB. Read by `02_msa_crosswalk.R`. |
| `outputs/tables/`, `outputs/figures/` | The rendered tables and figures the chapter uses. Regenerated by the pipeline; committed so the document compiles without a full run. |
| `nuthesis.cls`, `aea.bst`, `bibliography.bib` | LaTeX class, bibliography style and references for `Immigration_Card.tex`. |

### You must obtain

**An IPUMS USA extract built to the specification below.** Register at
[ipums.org](https://usa.ipums.org), build the extract, and download it as a
fixed-width `.dat.gz`. IPUMS does not permit redistribution of its microdata,
so neither the data nor the DDI codebook that documents it is committed here;
you generate both when you build your own extract.

The specification below is read from the DDI of the extract this analysis
used, so an extract built to it is equivalent. **Nineteen samples:**

- 1990 Census 5 percent (IPUMS sample `199001`)
- 2000 Census 5 percent (`200001`)
- ACS 1-year for every year 2006 through 2019, then 2021, 2022 and 2023

There is no 2010 decennial sample. The 2010 Census was short-form only, so the
2010 observations are the 2010 ACS, as are every year from 2006 onward. 2020 is
absent: pandemic-year collection departed from standard sampling procedures.

**Eighty-three variables.** The pipeline reads more of these than the headline
specification uses — `OCC2010` builds `Occup_Group` and `YNGCH` builds
`ChildAgeGroup`, for instance — so request the full list rather than a subset:

```
YEAR, SAMPLE, SERIAL, CBSERIAL, NUMPREC, HHWT, CLUSTER, STATEFIP,
COUNTYICP, COUNTYFIP, PUMA, PUMASUPR, METAREA, METAREAD, MET2013,
MET2013ERR, CITY, CITYPOP, STRATA, GQ, PERNUM, PERWT, NCHILD, NCHLT5,
ELDCH, YNGCH, RELATE, RELATED, SEX, AGE, MARST, BIRTHYR, RACE, RACED,
BPL, BPLD, ANCESTR1, ANCESTR1D, ANCESTR2, ANCESTR2D, CITIZEN, YRIMMIG,
YRSUSA2, LANGUAGE, LANGUAGED, SPEAKENG, SCHOOL, EDUC, EDUCD, GRADEATT,
GRADEATTD, CLASSWKR, CLASSWKRD, OCC, OCC2010, IND, WKSWORK2, UHRSWORK,
AVAILBLE, WORKEDYR, INCTOT, INCWAGE, INCWELFR, INCINVST, INCEARN,
POVERTY, MIGPLAC5, MIGPLAC1, MIGPUMA5, MIGPUMA1, MOVEDIN, VETSTAT,
VETSTATD, VETVIETN, VETKOREA, VETWWII, VETOTHER, VETOTHERD, DEPARTS,
EDUC_MOM, EDUC_POP, EDUC_MOM2, EDUC_POP2
```

No case selections are applied at the IPUMS end; every sample restriction
happens in `scripts/data_construction/04_sample_construction.R` and
`05_stable_msa_selection.R`.

**Where the two files go.** `ipumsr` reads the `.dat.gz` and the `.xml` DDI
together, and the pipeline reads the DDI from two locations: put the `.dat.gz`
and a copy of the `.xml` in the project root, and a second copy of the `.xml`
in `Codebooks/`. `01_ipums_panel_load.R` reads the root copy, while
`00_setup_analysis.R` and `08_load_analysis_subsets.R` read the `Codebooks/`
copy to build the `BPLD` origin-label lookup.

Your extract number will differ from the one these scripts were written
against. Update the filename in three places:
`scripts/data_construction/01_ipums_panel_load.R`,
`scripts/analysis/00_setup_analysis.R` and
`scripts/analysis/08_load_analysis_subsets.R`.

### Running it

Set the working directory to the project root, then:

```r
source("scripts/data_construction/run_all.R")   # 01-07: builds analysis_data.RData
source("scripts/analysis/run_all.R")            # 08-29: estimation, tables, figures
```

The two halves are a continuous 01--29 sequence. Data construction is the
expensive part; the analysis half runs in minutes. Intermediate checkpoints go
to `local_checkpoint_dir` (set at the top of both `00_setup` files), which
points outside the repository by default so large `.RData` images are never
synced or committed.

Then compile `Immigration_Card.tex`. Its `\graphicspath` and `\input@path`
search `outputs/`, so the tables and figures resolve wherever the document is
built.

---

## Models Estimated

| Model | Outcome | Sample | SE Type |
|-------|---------|--------|---------|
| Binary IV (LPM) | Educ_Binary (College+Grad vs non_HS+HighSch) | Full | MSA-clustered HC1 |
| Binary IV (LPM) | Educ_Binary | Full | State-clustered HC1 |
| Binary GLM (2SRI) | Educ_Binary | Full | Model-based |
| Binary IV (LPM) | HighSch vs non_HS | Restricted | MSA + State clustered HC1 |
| Binary IV (LPM) | College vs non_HS | Restricted | MSA + State clustered HC1 |
| Binary IV (LPM) | Grad vs non_HS | Restricted | MSA + State clustered HC1 |
| Multinomial (2SRI) | Educ_Group (4 levels) | Full | Bootstrap MSA-clustered (B=100, trimmed SD) |

---

## Outputs

`outputs/tables/` and `outputs/figures/` **are committed**, so the chapter
compiles from a fresh clone without running the pipeline. The `.tex` and image
files are what the document uses; the `.csv` files beside them are the
estimation results the rendering scripts read, committed so each table can be
checked against its source rather than only reproduced by re-running the
estimation.

Large intermediates are excluded: `analysis_data.RData` and the IPUMS extract
(both regenerable), and the per-section checkpoints, which are written to
`local_checkpoint_dir` outside the repository.

| Output | Written by | Used by |
|---|---|---|
| `main_results.csv` | `09` | `24` |
| `exposure_missing_check.csv` | `10` | `24` |
| `birth_cohort_final_spec.csv` | `11` | `24` |
| `rotemberg_results.csv`, `rotemberg_sensitivity.csv` | `15` | `24` |
| `adao_placebo_results.csv` | `16` | `24` |
| `attribution_summary.csv` | `19` | `24` |
| `results_migration_iv.csv` | `20` | `24` |
| `results_wage_iv.csv` | `21` | `24` |
| `demographic_summary_table.csv`, `mexico_vs_other.csv`, `parent_educ_missingness.csv` | `22` | `25`, `26` |
| `extract_composition.csv` | `08` | `26` |
| `multinomial_final_spec_test.csv` | `14b` | `24` |
| 16 `table_*.tex` | `23`, `24`, `25`, `26` | the chapter |
| `Educ_Prop*.png`, `Educ_By_Origin.png` | `27` | the chapter |
| `fig_origin_world.pdf` | `28` | the chapter |
| `fig_destination_*.pdf` | `29` | the chapter |

The chunk-loop outputs this section used to list (`chunk_outputs/`,
`meta_outputs_fst/`, `results_multinom_boot.csv`) belong to the retired
design and are no longer produced.

---

## Requirements

### R Packages

```r
install.packages(c(
  "tidyverse", "data.table", "AER", "nnet", "sandwich", "lmtest",
  "broom", "boot", "metafor", "knitr", "kableExtra", "modelsummary",
  "ggplot2", "fst", "fixest", "ipumsr", "scales", "readxl", "readr",
  "stringr", "usmap", "sf", "rnaturalearth", "RColorBrewer"
))
```

`fixest` carries the IV estimation, `ipumsr` reads the DDI and extract, and
`usmap`, `sf` and `rnaturalearth` are needed only by the three map scripts
(`28`, `29`). `usmap` must be 0.7 or later: the destination maps use its `sf`
county geometry, and the Connecticut handling in `29` depends on that version
carrying the 2022 planning regions.

### R Version

Developed under R 4.x on Windows 11 (64-bit), 32 GB RAM. The binding
constraint is data construction, which holds the full IPUMS extract in memory;
16 GB is a practical minimum and 32 GB is comfortable. The analysis half works
on `analysis_data.RData` and is far lighter.

The chunk-loop and meta-analytic pooling design this project used earlier --
and the multinomial bootstrap that made it slow -- was retired in October 2026
in favour of full-sample estimation. Those scripts are in
`scripts/analysis/archive/`. Nothing in the current pipeline takes hours.

---

## Identification Strategy

The instrument is the individual-level average shift-share predicted immigrant
inflow over each native's ages 15--25:

```
avg_exposure_i = (1/11) * sum_{a=15}^{25} predicted_inflow_{m, t(a)}
```

where `predicted_inflow_{m,t}` is the Card (2001) shift-share prediction:

```
predicted_inflow_{m,t} = sum_k s_{m,k,1990} * Delta_M_{k,t}
```

`s_{m,k,1990}` is the share of all U.S. immigrants from origin `k` in 1990
located in MSA `m` (national share definition, group_by BPLD not MSA), and
`Delta_M_{k,t}` is the national inflow of immigrants from origin `k` in year
`t`. All inflow variables are rescaled to units of 100 persons before
modeling.

---

## Key Literature

- Card (2001) -- shift-share IV construction and national share definition
- Goldsmith-Pinkham, Sorkin, and Swift (2020) -- share-based IV validity
- Borusyak, Hull, and Jaravel (2022) -- shock-based IV validity
- Adão, Kolesár, and Morales (2019) -- shift-share inference and clustering
- Hunt (2017) -- closest predecessor on immigration and native education
- McHenry (2015) -- immigration and native human capital
- Peri and Sparber (2009) -- task specialization and complementarity mechanism

---

## Stable MSA Sample

The analysis restricts to natives residing in the 175 MSAs that appear
consistently across all survey years (1990--2023). MSA codes use the
MET2013 definition from IPUMS, crosswalked to CBSA codes where necessary.
The 175 stable MSAs span 45 states and the District of Columbia. Six states
(Montana, North Dakota, South Dakota, Vermont, West Virginia, Wyoming) have
no qualifying MSAs.

---

## Contact

Joe Sandow  
jsandow2@huskers.unl.edu  
Department of Economics  
University of Nebraska-Lincoln
