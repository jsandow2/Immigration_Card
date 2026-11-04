# Immigration_Card modular pipeline

This is a section-by-section split of `Immigration_Card_Data_Construction.R`
and `Immigration_Card_Analysis.R` into the same numbered-script style as the
Tinuiti project's `scripts/` folder, adapted for this project's much larger
data (a ~2GB `analysis_data.RData` extract and multi-hour model runs, versus
small CSVs).

**Nothing in your live `Immigration_Card_Data_Construction.R` or
`Immigration_Card_Analysis.R` was changed.** These are new files, built by
extracting each script's own `# Section N` blocks verbatim -- the original
two files remain the source of truth until you decide to switch over.

## Layout

```
data_construction/
  00_setup.R                     libraries + rm(list=ls()) -- sourced by every 0X script;
                                  also sets local_checkpoint_dir (see below)
  01_ipums_panel_load.R          Section 1
  02_msa_crosswalk.R             Section 2
  03_variable_construction.R     Section 3
  04_sample_construction.R       Section 4
  05_stable_msa_selection.R      Section 5
  06_shift_share_iv.R            Section 6
  07_exposure_iv.R               Section 7 -- ends by saving analysis_data.RData
  run_all.R                      sources 00 through 07 in order

  (Raw 1990 PUMS ingestion, formerly Section 1, was removed 2026-09-07 --
  unused: the IPUMS panel extract loaded in Section 1 above already covers
  1990 as a native 5% Census sample in rectangular format. Archived at
  archive/scripts/01_pums_1990_ingestion_UNUSED.R.)

analysis/
  00_setup_analysis.R            shared library() calls + analysis_checkpoint_dir,
                                  sourced by 09 onward
  08_load_analysis_subsets.R     Section 8 -- fresh session, loads analysis_data.RData
  09_run_models_function.R       Section 9
  10_meta_chunk_loop.R           Sections 10-11 (the long-running chunk loop; no
                                  checkpoint saved -- real output is the .fst files)
  11_meta_pooling.R              Section 11 (originally 12b)
  20_migration_iv.R              Section 12
  21_wage_iv.R                   Section 13
  24_results_tables.R            Section 14
  27_education_figures.R         Section 15
  23_msa_stability.R             Section 16
  10_exposure_missing_check.R    Section 17
  18_birth_cohort_check.R        Section 18 -- birth-cohort FE on the exposure_missing==0
                                  subsample, via Section 17's feols() approach (no chunk-loop re-run)
  15_rotemberg.R                 Section 19
  16_adao_placebo.R              Section 20
  21_descriptives_mexico_other.R Section 21
  22_demographic_summary.R       Section 22
  23_estimation_sample_summary.R Section 23 -- Table 1-style summary stats for the
                                  all_natives_175 estimation sample, loaded straight from
                                  after_sec08 (no chunk-loop dependency, same reasoning as
                                  Section 17); no checkpoint saved (standalone leaf)
  25_dissertation_table_helper.R write_dissertation_table() + the table_demographic_summary
                                  example; terminal script, no checkpoint saved
  run_all.R                      sources 08 through 24 in order

Renumbered 2026-09-07: analysis/ continues data_construction/'s 01-07
sequence (old 09-21 are now 08-24) rather than restarting at 01, and its
save.image()-per-section pattern was replaced with the same explicit-save()
approach data_construction/ already uses (see "How checkpointing works"
below). Two sections that deliberately skip most of the chain and load
straight from Section 8 (10_exposure_missing_check.R and
23_estimation_sample_summary.R) only see a live meta_iv_msa (Section 11's
pooled IV estimate) if Section 11 already ran earlier in the same R
session -- load() is additive and never clears the environment, so this was
true before the renumbering too and is not something this rewrite changed.
REMINDER: once the Section 10 chunk loop is re-run under the current
raw-count inflow-series construction (see 06_shift_share_iv.R), update the
hardcoded fallback value in 10_exposure_missing_check.R (currently
3.054701e-05) to the new meta_iv_msa.
```

## How checkpointing works

As of 2026-09-07 both halves use the same explicit-save() pattern; neither
uses `save.image()` anymore.

Each script does `source("00_setup_analysis.R")` (or, for
`08_load_analysis_subsets.R`/`00_setup.R`, its own inline setup) then
`load(file.path(analysis_checkpoint_dir, "after_secNN.RData"))` -- but
unlike a `save.image()`-based chain, this restores only the specific objects
the previous section explicitly chose to hand off, not everything that was
in memory at that point. At the end, each script calls
`save(obj1, obj2, ..., file = file.path(analysis_checkpoint_dir,
"after_secXX.RData"))` naming only the objects a later section actually
reads (see the per-section object list in each script's own header/footer
comments). Some scripts don't chain from their immediate predecessor at all
-- e.g. `10_meta_chunk_loop.R` and `11_meta_pooling.R` both load straight
from `08`'s checkpoint, since nothing Section 9 or 10 produces is needed
further downstream -- and a few (`10_meta_chunk_loop.R`,
`23_estimation_sample_summary.R`, `25_dissertation_table_helper.R`) save no
checkpoint at all, either because their real output is entirely `.fst`/CSV
files on disk or because they're a standalone leaf nothing else reads from.

`analysis_checkpoint_dir` (a subfolder of `data_construction/`'s
`local_checkpoint_dir`, set in `00_setup_analysis.R` and duplicated inline in
`08_load_analysis_subsets.R`, currently
`C:/R_checkpoints/Immigration_Card/analysis`) is local-only, outside
OneDrive, same rationale as `data_construction/`'s checkpoints: these exist
only to resume an interrupted run or skip straight to a later section on the
same machine, not for cross-device reach. It's a separate subfolder from
`data_construction/`'s checkpoints (which sit directly in
`local_checkpoint_dir`) so the two halves don't collide writing e.g.
`after_sec01.RData` to the same path under two different meanings.

Two scripts (`10_exposure_missing_check.R`, `23_estimation_sample_summary.R`)
deliberately skip most of the chain and load straight from `08`'s
checkpoint, to avoid depending on the multi-hour Section 10 chunk loop
having run first. Because `load()` is additive (it never clears the
environment), these two only see a live `meta_iv_msa` (Section 11's pooled
IV estimate) if Section 11 already ran earlier in the same R session --
otherwise `10_exposure_missing_check.R` falls back to a hardcoded snapshot
value. This was true before the 2026-09-07 renumbering too and isn't a
regression from it. See the REMINDER above the file listing: that hardcoded
value needs updating once the chunk loop is re-run under the current
inflow-series construction.

Three named outputs are untouched by any of this and still get written
directly into this OneDrive-synced project folder, since other tooling (and
the analysis pipeline itself) references them by exact filename:
`all_data.RData` (written by `03_variable_construction.R`),
`reduced_data.RData` (written by `04_sample_construction.R`), and
`analysis_data.RData` (written by `07_exposure_iv.R`).

## Running it

**Run everything from the project root** (`Immigration_Card/`), not from
inside `scripts/` -- exactly where `Immigration_Card.Rproj` already puts you
when you open it. All the original data paths (`analysis_data.RData`,
`Codebooks/usa_00007.xml`, output CSVs, etc.) are unchanged from the
monolithic scripts and assume the project root as the working directory;
every `source()`/`load()`/`save.image()` call these modular scripts add for
their own chaining is written as a full `scripts/...` path for exactly this
reason, so `setwd()`-ing into a subfolder will break the data paths even
though the `source()` calls would still resolve.

Either source `run_all.R` in each folder for a full reproduction, or run
individual numbered scripts once the checkpoint before them exists -- e.g.
to redo just the wage IV models:

```r
# from the project root
source("scripts/analysis/00_setup_analysis.R")
load(file.path(analysis_checkpoint_dir, "after_sec12.RData"))
source("scripts/analysis/21_wage_iv.R")
```

`11_meta_pooling.R` is the one exception worth knowing about: per its own
header comment, it only reads `meta_outputs_fst/*.fst` and can be run in a
brand-new session once the chunk loop has finished, without needing
`10_meta_chunk_loop.R`'s checkpoint at all (there isn't one -- it loads
straight from `08`'s checkpoint instead).

## Renumbering, 2026-10-21

Second renumbering pass. Eleven scripts that produced no artefact the paper
uses, and backed no claim it makes, moved to `analysis/archive/`:
`20h_effective_clusters`, `20i_inflow_rate_spec`, `20l_wild_cluster_bootstrap`,
`20m_windowed_endogenous_test`, `25_level_difference_test`,
`32_no_wave_divisor`, `41c_test_americas_annual`,
`23_estimation_sample_summary` (wrote nothing), and the three `test_*` scratch
scripts.

The 24 that remain run as a continuous sequence with `data_construction/`,
which ends at 07. The whole pipeline therefore reads **01 through 29** in
order and can be run start to finish by someone who did not build it.

| block | contents |
| --- | --- |
| `00`, `00b` | setup; shared table helper |
| `08` | data loading |
| `09`-`14` | estimation on the main design (`14b` decomposes the binary outcome) |
| `15`-`19` | shift-share identification and diagnostics |
| `20`-`21` | Card (2001) replications |
| `22`-`23` | descriptives and panel composition |
| `24`-`29` | rendering: LaTeX tables and figures |

Old to new:

| 2026-09-07 name | 2026-10-21 name |
| --- | --- |
| `08_load_analysis_subsets.R` | `08_load_analysis_subsets.R` |
| `29_main_results.R` | `09_main_results.R` |
| `17_exposure_missing_check.R` | `10_exposure_missing_check.R` |
| `18b_birth_cohort_final_spec.R` | `11_birth_cohort_final_spec.R` |
| `28_final_spec_no_age_fe.R` | `12_final_spec_no_age_fe.R` |
| `26_two_margin.R` | `13_two_margin.R` |
| `31b_gender_difference_test.R` | `14_gender_difference_test.R` |
| `19_rotemberg.R` | `15_rotemberg.R` |
| `20_adao_placebo.R` | `16_adao_placebo.R` |
| `20b_asian_placebo_concentration.R` | `17_asian_placebo_concentration.R` |
| `20j_rate_spec_loo.R` | `18_rate_spec_loo.R` |
| `34_attribution_stability_check.R` | `19_attribution_stability_check.R` |
| `12_migration_iv.R` | `20_migration_iv.R` |
| `13_wage_iv.R` | `21_wage_iv.R` |
| `21_origin_composition.R` | `22_origin_composition.R` |
| `16_msa_stability.R` | `23_msa_stability.R` |
| `14_results_tables.R` | `24_results_tables.R` |
| `24_dissertation_table_helper.R` | `25_dissertation_table_helper.R` |
| `24b_appendix_tables.R` | `26_appendix_tables.R` |
| `15_education_figures.R` | `27_education_figures.R` |
| `40_map_origin_world.R` | `28_map_origin_world.R` |
| `41_map_destination_us.R` | `29_map_destination_us.R` |

`00_setup_analysis.R` and `00b_table_helper.R` keep their numbers; they are
sourced by name throughout, and `08` keeps its own because the sequence
continues from `data_construction/`. The criterion for keeping a script was
that it either produces a table or figure the paper uses, feeds one that does,
or supports a claim the text makes in prose -- every claim has to be
defensible from the data or from a citation.
