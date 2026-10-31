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

## Logs and timings

`outputFolder` holds a timestamped log of every `execute()` call (`log_<database>.txt`), the duration of each step
(`steps_<database>.csv`: cohort generation, extraction, the plpData cache, population, splits, model runs) and, per
model run, PLP's log (`plpLog.txt` in the run folder). Each worker process logs to
`workers/log_<database>_<pid>.txt` and keeps its temporary files in `outputFolder/tmp` (emptied at the next start).
`Rscript extras/production/timingReport.R` summarises where the time goes and writes `timing_<database>.csv`
(aggregate, fine to send).

## What leaves the site

`exportResults()` writes one zip: per-run metrics (AUC, log loss, Brier), counts, selected hyperparameters,
runtimes and package versions, with counts below `minCellCount` shown as `-minCellCount`. It contains no
patient-level data, predictions, models or credentials. `outputFolder/data` and `outputFolder/models` stay on site.
