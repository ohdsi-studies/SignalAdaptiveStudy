# SignalAdaptiveStudy: production run (README_PRODUCTION.md). Run step by step, in a fresh R session.
# Site settings live in extras/localSettings.R (see below); no need to edit this file. Run once per database:
# change databaseName in extras/localSettings.R and run again.
# Nothing patient-level leaves outputFolder/data and outputFolder/models; only the zip from exportResults() is shared.

# The Databricks JDBC driver (3.x) reads results with Apache Arrow, which on Java >= 16 needs this permission (without
# it, querySql() fails with a Java NullPointerException). It must be set before Java starts, i.e. before the library()
# call below, in a fresh R session.
options(java.parameters = "--add-opens=java.base/java.nio=org.apache.arrow.memory.core,ALL-UNNAMED")

library(SignalAdaptiveStudy)

# ---- Connection ---------------------------------------------------------------------------------------------------
# Credentials from environment variables, never in this file:
#   DATABASE_SERVER    Databricks host (xxx.cloud.databricks.com), or a full jdbc:databricks://... connection string
#   DATABASE_PATH      the SQL warehouse's HTTP path (/sql/1.0/warehouses/...); not needed with a full string
#   DATABASE_PASSWORD  a personal access token (user "token")
#   DATABASE_DRIVER    the folder with the JDBC driver (DatabaseConnector::downloadJdbcDrivers("spark", pathToDriver = ...))
server <- Sys.getenv("DATABASE_SERVER")
connectionString <- if (startsWith(server, "jdbc:")) server else
  sprintf("jdbc:databricks://%s:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=%s;", server,
          Sys.getenv("DATABASE_PATH"))
# Arrow result transfer, with the Java permission above (the combination that works with driver 3.4.3)
if (!grepl("EnableArrow=", connectionString, ignore.case = TRUE)) {
  connectionString <- paste0(sub(";?$", ";", connectionString), "EnableArrow=1;")
}
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "spark",                                         # other platforms: change here
  connectionString = connectionString,
  user = "token",
  password = Sys.getenv("DATABASE_PASSWORD"),
  pathToDriver = Sys.getenv("DATABASE_DRIVER")
)
# ---- Site settings ----------------------------------------------------------------------------------------------
# In extras/localSettings.R (git-ignored, so pulls never touch it): copy extras/localSettings_template.R to
# extras/localSettings.R and fill it in. It sets scratchSchema, minCellCount, cdmSchemas and databaseName (the
# database of this run).
source("extras/localSettings.R")
options(sqlRenderTempEmulationSchema = scratchSchema)     # Spark has no temp tables: they go to the scratch schema
cdmDatabaseSchema <- cdmSchemas[[databaseName]]
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
# SAR, SAL, LASSO and ridge, budgets 100/250/500/750/1000 events, 3 orderings, within-database train/test.
# Tasks with more than 2 million persons are sampled at extraction to 2 million, or more where needed for 1.5 x 1000
# expected training-pool events (addSampleSizes()); a task runs the budgets its pool supports.
# Original tasks only: the new catalogue tasks get counts but no models, so they stay unseen.
# Restartable: rerun the same lines after an interruption and finished runs are skipped.
counts <- execute(dataSource, tasks, outputFolder, mode = "counts")   # the 3b counts (read back, not recounted)
modelTasks <- addSampleSizes(originalCohortTasks(), counts[counts$feasible, ])
if (exists("onlyTasks")) modelTasks <- modelTasks[modelTasks$taskId %in% onlyTasks, ]   # optional, localSettings.R
for (i in seq_len(nrow(modelTasks))) {
  task <- modelTasks[i, ]
  execute(dataSource, task, outputFolder, mode = "models", budgets = task$budgets[[1]], orderings = 1:3,
          methods = studyMethods())
}
exportResults(outputFolder, minCellCount = minCellCount)  # send this zip after 3c
