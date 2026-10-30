# SignalAdaptiveStudy: production run (README_PRODUCTION.md). Edit the lines marked EDIT, then run step by step.
# Nothing patient-level leaves outputFolder/data and outputFolder/models; only the zip from exportResults() is shared.

library(SignalAdaptiveStudy)

# ---- EDIT: database ---------------------------------------------------------------------------------------------
# Credentials from environment variables (or keyring), never in this file.
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = Sys.getenv("SAR_DBMS"),                          # e.g. "postgresql", "redshift", "spark", "duckdb"
  server = Sys.getenv("SAR_SERVER"),
  user = Sys.getenv("SAR_USER"),
  password = Sys.getenv("SAR_PASSWORD"),
  port = Sys.getenv("SAR_PORT", unset = NA),
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)
cdmDatabaseSchema <- "cdm"                                # EDIT: CDM and vocabulary schema
cohortDatabaseSchema <- "scratch"                         # EDIT: a schema you can write to
cohortTable <- "sar_study_cohort"                         # EDIT if the name is taken
databaseName <- "CCAE"                                    # EDIT: short name used in the results
options(sqlRenderTempEmulationSchema = NULL)              # EDIT on Spark/Oracle/BigQuery: a writable schema
outputFolder <- file.path(getwd(), "SignalAdaptiveStudy_output", databaseName)   # EDIT if wanted
minCellCount <- 5                                         # EDIT: your site's minimum cell count
# -----------------------------------------------------------------------------------------------------------------

dataSource <- createCdmDataSource(connectionDetails, cdmDatabaseSchema, cohortDatabaseSchema, databaseName,
                                  cohortTable = cohortTable)
tasks <- studyTasks()          # the 26 original tasks and 4 new catalogue tasks (PLP benchmark definitions)

# ---- Step 3b: feasibility counts, no models ------------------------------------------------------------------------
# Generates the cohorts, extracts each task's data (cached in outputFolder/data) and counts the training pool and test
# set. Feasible: >= 500 training events and >= 80 test events.
counts <- execute(dataSource, tasks, outputFolder, mode = "counts")
print(counts[, c("taskId", "n", "poolEvents", "testEvents", "feasible")])
exportResults(outputFolder, minCellCount = minCellCount)  # send this zip after 3b

# ---- Step 3c: models on the feasible tasks -----------------------------------------------------------------------
# SAR, SAL, LASSO and ridge, budgets 100/250/500, 3 orderings, within-database train/test.
# Restartable: rerun the same lines after an interruption and finished runs are skipped.
# Original tasks only: the new catalogue tasks get counts but no models, so they stay unseen.
methods <- studyMethods()
originalTasks <- originalCohortTasks()
feasible <- originalTasks[originalTasks$taskId %in% counts$taskId[counts$feasible], ]
results <- execute(dataSource, feasible, outputFolder, mode = "models", budgets = c(100, 250, 500), orderings = 1:3,
                   methods = methods)
exportResults(outputFolder, minCellCount = minCellCount)  # send this zip after 3c
