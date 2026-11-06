# The PLP benchmark tasks with the final study's settings: two feature sets, budgets 100 to 1000, 3 orderings, SAR, LASSO and ridge,
# with PLP's default LASSO, elastic net and two reduced SARs as secondary methods. Writes a zip of aggregates per database, to send
# back. Restartable: rerun after an interruption and finished runs are skipped.
#
#   Rscript extras/production/runBenchmarks.R               every database in cdmSchemas, one after the other
#   Rscript extras/production/runBenchmarks.R MDCD MDCR     only these
# (in the project folder; settings as for CodeToRun.R). A database that fails does not stop the others.

# Needed by the Databricks JDBC driver, before Java starts (see CodeToRun.R)
options(java.parameters = "--add-opens=java.base/java.nio=org.apache.arrow.memory.core,ALL-UNNAMED")
library(SignalAdaptiveStudy)

server <- Sys.getenv("DATABASE_SERVER")
connectionString <- if (startsWith(server, "jdbc:")) server else
  sprintf("jdbc:databricks://%s:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=%s;", server, Sys.getenv("DATABASE_PATH"))
if (!grepl("EnableArrow=", connectionString, ignore.case = TRUE)) {
  connectionString <- paste0(sub(";?$", ";", connectionString), "EnableArrow=1;")
}
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "spark", connectionString = connectionString, user = "token", password = Sys.getenv("DATABASE_PASSWORD"),
  pathToDriver = Sys.getenv("DATABASE_DRIVER"))

source("extras/localSettings.R")   # scratchSchema, minCellCount, cdmSchemas, workers
options(sqlRenderTempEmulationSchema = scratchSchema)
args <- commandArgs(trailingOnly = TRUE)
databases <- if (length(args)) args else names(cdmSchemas)
stopifnot(all(databases %in% names(cdmSchemas)))
message("Cyclops ", utils::packageVersion("Cyclops"), " commit ", utils::packageDescription("Cyclops")$RemoteSha)

runDatabase <- function(database) {
  cohortTable <- paste0("sar_benchmark_cohort_", tolower(database))   # its own cohort table, next to any other run
  outputFolder <- file.path(getwd(), "SignalAdaptiveStudy_benchmarks", database)
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
  dataSource <- createCdmDataSource(connectionDetails, cdmSchemas[[database]], scratchSchema, database, cohortTable = cohortTable,
                                    covariateSettings = finalCovariateSettings(scratchSchema, cohortTable),
                                    cohortDefinitionSet = benchmarkCohortDefinitionSet(),
                                    incremental = TRUE, combinedCohorts = combinedFeatureCohorts())
  counts <- execute(dataSource, benchmarkTasks(), outputFolder, mode = "counts")
  print(counts[, c("taskId", "n", "poolEvents", "testEvents")])
  feasible <- feasibleCounts(counts)
  message(sprintf("%s: %d of %d benchmark tasks have >= 1000 pool and >= 333 test events", database, nrow(feasible), nrow(benchmarkTasks())))
  tasks <- addSampleSizes(benchmarkTasks(), feasible, budgets = finalBudgets())
  tasks <- expandFeatureSets(tasks)
  if (exists("onlyTasks")) tasks <- tasks[tasks$taskId %in% onlyTasks, ]   # optional, localSettings.R
  execute(dataSource, tasks, outputFolder, mode = "models", orderings = 1:3, methods = finalStudyMethods(),
          workers = if (exists("workers")) workers else 1)
  exportResults(outputFolder, minCellCount = minCellCount)
}
results <- vapply(databases, function(database) {
  tryCatch(as.character(runDatabase(database)), error = function(e) paste("FAILED:", conditionMessage(e)))
}, "")
message(paste(sprintf("%s: %s", databases, results), collapse = "\n"))
