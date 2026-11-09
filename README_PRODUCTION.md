# SignalAdaptiveStudy: production run

Replicates signal-adaptive ridge (SAR) against LASSO and ridge on the original prediction tasks, with
PatientLevelPrediction on an OMOP CDM. Two steps: feasibility counts (3b), then models (3c). Only aggregate CSVs leave
the site.

## Install

```r
# In a clone of this repository. R >= 4.4, a C++ compiler (Cyclops is built from GitHub at a pinned commit), Java for
# DatabaseConnector's JDBC drivers
install.packages("renv")
renv::restore()   # runtime packages only, binaries from Posit Package Manager; only Cyclops compiles (about a minute)
install.packages(".", repos = NULL, type = "source")   # the study package itself
```

Tests and development use the `dev` renv profile (adds testthat, Eunomia and the like): `Sys.setenv(RENV_PROFILE = "dev")`
before `renv::restore()`.

**No internet on the server:** on a connected machine with the same OS and R version, run `renv::restore()` and
then copy the project folder together with the renv cache (`renv::paths$cache()`), or build source tarballs with
`renv::snapshot()` + `renv::install()` from a local repository. JDBC drivers: `DatabaseConnector::downloadJdbcDrivers()`.

## Connection

Credentials come from environment variables (`.Renviron` or the system), never code: `DATABASE_SERVER` (the
Databricks host, or a full `jdbc:databricks://...` string), `DATABASE_PATH` (the SQL warehouse's HTTP path),
`DATABASE_PASSWORD` (a personal access token) and `DATABASE_DRIVER` (the folder with the JDBC driver from
`DatabaseConnector::downloadJdbcDrivers("spark", pathToDriver = ...)`). You need read access to the
CDM schema and write access to one schema for the cohort tables.

## Smoke test (no database, about a minute)

```sh
Rscript extras/production/smokeTest.R
```

This should print `smoke test PASSED`. Send back the printed table.

## Steps

Copy `extras/localSettings_template.R` to `extras/localSettings.R` (git-ignored, so `git pull` never touches it) and fill
in the scratch schema, the CDM schema of each database and the database of this run. Then run `extras/CodeToRun.R`
section by section, starting in a fresh R session (its first line sets a Java option the Databricks driver needs, which
only takes effect before Java starts).

1. **3b, counts:** generates 26 cohorts (Circe SQL shipped in `inst/cohorts`) and counts the training pool and test
   set of 30 tasks: the 26 original tasks and 4 new ones (PLP benchmark definitions). Cohorts and outcomes only, no
   covariates and no models. Send the zip.
2. **3c, models:** SAR, SAL, LASSO and ridge on the feasible original tasks (>= 500 training events, >= 80 test events), budgets
   100/250/500/750/1000 (as far as each task's training pool allows), 3 orderings. Runs in `workers` R processes
   (`extras/localSettings.R`, default 6; [mirai](https://mirai.r-lib.org)), one thread each, with the Andromeda (duckdb)
   memory of all workers held to half the RAM; the data are extracted in the main session. A worker that dies is
   replaced and its runs tried once more. Restartable: rerun after an interruption and finished runs are skipped.
   Send the zip.

**Runtime:** the model runs dominate (SAR and SAL about 5x LASSO and ridge), plus data extraction of minutes per task.
With workers the wall-clock time falls about in proportion to the physical cores (8 vCPUs are often 4 cores with
hyperthreading). Disk: the plpData cache is a few GB, plus the workers' temporary copies of the task data.

## Screening candidate tasks

`Rscript extras/production/runTaskScreening.R` checks the candidate tasks of the follow-up study: cohorts and counts only, no
covariates and no models. With no arguments it goes through every database in `cdmSchemas`, one after the other; name databases to run
only those (`Rscript extras/production/runTaskScreening.R MDCD MDCR`). A database that fails does not stop the others, and the last
message lists each database's result. Little CPU or memory is needed, as the work runs in the database. It is safe next to a running
study, because it uses its own cohort table,
`sar_screen_cohort_<database>`. For each database it writes
`SignalAdaptiveStudy_screening/<database>/SignalAdaptiveStudy_screening_<database>_<date>.zip`; send those zips. Rerun after an
interruption: cohorts and tasks already done are skipped, so a database screened before only runs what is new.

The tasks marked `feasibility` in `inst/settings/taskScreening/tasks.csv` use rough cohorts, not study phenotypes: they show whether a task
has enough events, not whether its definition is right.

## PLP benchmark tasks

`Rscript extras/production/runBenchmarks.R` runs the 13 PLP benchmark tasks with the final study's settings, on every database in
`cdmSchemas` (or those you name: `Rscript extras/production/runBenchmarks.R MDCD`). Start with one small database.

**Cyclops:** the study runs on the latest `develop` of Cyclops, which it also helps to test. Before a run, update it and reinstall the
package:

```r
renv::install("OHDSI/Cyclops@develop")
install.packages(".", repos = NULL, type = "source")
```

The script prints the Cyclops commit, and each result row records it in `versions`, with each fit's algorithm, iterations and status. Finished runs are not redone when Cyclops changes, so
use a new output folder if results should all come from one commit.

**Settings:**

- **Tasks:** PLPBenchmarks' cohorts, risk windows and population rules, as defined there. A task runs in a database with at least 1,000
  training-pool and 333 test events. Tasks with repeat entries per person (the COVID tasks) are split by person.
- **Feature sets**, all binary and all penalised:
  - `defaults`: FeatureExtraction's defaults without its continuous scores (Charlson, DCSI, CHADS2, CHA2DS2-VASc);
  - `phenotypes`: gender, 5-year age groups and the 64 phenotypes of the PLP vignette "Constrained Predictors", with the 14
    antibiotic classes combined into one (51 indicators).
- **Covariate window:** ends before the risk window starts. All benchmark tasks end at day -1, as PLPBenchmarks does.
- **Budgets and orderings:** 100, 250, 500 and 1000 events; 3 orderings.
- **Methods:** SAR, LASSO and ridge (primary); PLP's default LASSO, elastic net (mixing 0.5), SAR with its shape fixed at (0, 1, -1), and SAR
  without its third term (secondary).

It writes `SignalAdaptiveStudy_benchmarks/<database>/`; send the zip that the last message names. Its results also hold calibration, the
AUC's standard error, and the metrics on test rows with full follow-up.

## Comparing Cyclops's fitting algorithms

`Rscript extras/production/compareFitAlgorithms.R CCAE` fits one benchmark task (HFT2DM, both feature sets, one ordering, all budgets) with
SAR, SAR with a fixed shape, LASSO, ridge and elastic net, each by coordinate descent (`_ccd`) and by Cyclops's newer automatic algorithm
(`_auto`). It reuses the cohorts and data of `runBenchmarks.R` and writes `SignalAdaptiveStudy_fit_comparison/<database>/`. Send the zip that
the last message names: it includes how far apart the coefficients of the two algorithms are.

## Logs and timings

`outputFolder` holds a timestamped log of every `execute()` call (`log_<database>.txt`), the duration of each step
(`steps_<database>.csv`: cohort generation, extraction, the plpData cache, population, splits, model runs) and, per
model run, PLP's log (`plpLog.txt` in the run folder). Each worker process logs to
`workers/log_<database>_<pid>.txt` and keeps its temporary files in `outputFolder/tmp` (emptied at the next start).
`Rscript extras/production/timingReport.R` summarises where the time goes and writes `timing_<database>.csv`
(aggregate, fine to send).

## Restarting after a change

Every cached and finished output records the settings it was made with (`fingerprint.json`): a task's data cache
(cohort definitions, covariate settings, sample), its counts and each model run (also the population rule, split and
method settings). After a change of settings, a data cache is extracted again and counts are recounted
automatically; finished model runs made with other settings stop `execute()` with an error naming the run folders and
what changed (move them away, or restore the settings). Outputs from before fingerprints are accepted as they are.

## What leaves the site

`exportResults()` writes one zip: per-run metrics (AUC, log loss, Brier), counts, selected hyperparameters,
runtimes and package versions, with counts below `minCellCount` shown as `-minCellCount`. A run with `nonZero = 0`
predicts the same risk for everyone (a constant model), whether the intercept-only model was selected
(`interceptOnly = TRUE`) or a fitted model has no non-zero coefficient. It contains no
patient-level data, predictions, models or credentials. `outputFolder/data` and `outputFolder/models` stay on site.
