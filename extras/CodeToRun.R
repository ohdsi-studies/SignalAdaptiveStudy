# SignalAdaptiveStudy: production run (README_PRODUCTION.md). Edit the lines marked EDIT, then run step by step.
# Run once per database: change databaseName and cdmDatabaseSchema and run again.
# Nothing patient-level leaves outputFolder/data and outputFolder/models; only the zip from exportResults() is shared.

library(SignalAdaptiveStudy)

# ---- EDIT: connection ---------------------------------------------------------------------------------------------
# Credentials from environment variables (.Renviron or keyring), never in this file. Databricks:
#   SAR_CONNECTION_STRING = jdbc:databricks://<host>:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=<path>;
#   SAR_USER = token, SAR_PASSWORD = <personal access token>,
#   DATABASECONNECTOR_JAR_FOLDER = a folder with the driver (DatabaseConnector::downloadJdbcDrivers("spark"))
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "spark",                                         # EDIT for other platforms
  connectionString = Sys.getenv("SAR_CONNECTION_STRING"),
  user = Sys.getenv("SAR_USER"),
  password = Sys.getenv("SAR_PASSWORD"),
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)
scratchSchema <- "scratch"                                # EDIT: a schema you can write to (e.g. catalog.schema)
options(sqlRenderTempEmulationSchema = scratchSchema)     # Spark has no temp tables: they go to the scratch schema
minCellCount <- 5                                         # EDIT: your site's minimum cell count

# ---- EDIT: database (one run per database) ---------------------------------------------------------------------
databaseName <- "CCAE"                                    # EDIT: CCAE, MDCD, MDCR, OPTUM_EHR, OPTUM_CLAIMS, ...
cdmDatabaseSchema <- "cdm_ccae"                           # EDIT: that database's CDM and vocabulary schema
# -----------------------------------------------------------------------------------------------------------------

cohortTable <- paste0("sar_study_cohort_", tolower(databaseName))   # one cohort table per database
outputFolder <- file.path(getwd(), "SignalAdaptiveStudy_output", databaseName)
dataSource <- createCdmDataSource(connectionDetails, cdmDatabaseSchema, scratchSchema, databaseName,
                                  cohortTable = cohortTable)
tasks <- studyTasks()          # the 26 original tasks and 4 new catalogue tasks (PLP benchmark definitions)

# ---- Step 3b: feasibility counts, no models ------------------------------------------------------------------------
# Generates the cohorts (once per database) and counts each task's training pool and test set; cohorts and outcomes
# only, no covariates. Feasible: >= 500 training events and >= 80 test events. Restartable; a failing task leaves
# outputFolder/counts/<database>/<task>_error.txt and the others continue.
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
