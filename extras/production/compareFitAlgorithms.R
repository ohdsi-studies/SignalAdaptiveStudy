# Compares Cyclops's two ways of fitting logistic regression on one benchmark task: coordinate descent with the KKT criterion against
# the newer automatic algorithm (Newton and prox-Newton steps, stopped by a certified gap), both through the study's folds and search.
# Writes a zip of aggregates to send back, including how far apart the coefficients are, and Newton-CG's work per fit (Hessian-vector
# products, preconditioner). One ordering; by default all four budgets and both feature sets.
#
#   Rscript extras/production/compareFitAlgorithms.R CCAE
#   Rscript extras/production/compareFitAlgorithms.R CCAE --budgets=100,250,500 --featureSets=defaults
# (in the project folder; settings as for runBenchmarks.R). It reuses the cohorts and, if present, the extracted data of runBenchmarks.R
# and of earlier comparisons. The output folder is named after the Cyclops commit, so a new Cyclops starts afresh.

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
database <- args[1]
stopifnot(!is.na(database), database %in% names(cdmSchemas))
option <- function(name, default) {
  given <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(given)) strsplit(given, ",")[[1]] else default
}
budgets <- as.numeric(option("budgets", finalBudgets()))
featureSets <- option("featureSets", c("defaults", "phenotypes"))
stopifnot(all(budgets %in% finalBudgets()), all(featureSets %in% c("defaults", "phenotypes")))
cyclopsCommit <- utils::packageDescription("Cyclops")$RemoteSha %||% ""
message("Cyclops ", utils::packageVersion("Cyclops"), " commit ", cyclopsCommit, "; budgets ", paste(budgets, collapse = ", "),
        "; feature sets ", paste(featureSets, collapse = ", "))

cohortTable <- paste0("sar_benchmark_cohort_", tolower(database))   # shared with runBenchmarks.R, so the cohorts are reused
outputFolder <- file.path(getwd(), paste0("SignalAdaptiveStudy_fit_comparison_", if (nzchar(cyclopsCommit)) substr(cyclopsCommit, 1, 8) else "local"),
                          database)
dir.create(file.path(outputFolder, "data", database), recursive = TRUE, showWarnings = FALSE)
message("Output: ", outputFolder)

# Start from data extracted before (runBenchmarks.R, earlier comparisons); its fingerprint tells whether it still applies
earlierData <- c(file.path(getwd(), "SignalAdaptiveStudy_benchmarks", database, "data", database),
                 Sys.glob(file.path(getwd(), "SignalAdaptiveStudy_fit_comparison*", database, "data", database)))
for (folder in Sys.glob(file.path(setdiff(earlierData, file.path(outputFolder, "data", database)), "benchmark_HFT2DM_*"))) {
  if (!dir.exists(file.path(outputFolder, "data", database, basename(folder)))) {
    file.copy(folder, file.path(outputFolder, "data", database), recursive = TRUE)
  }
}

dataSource <- createCdmDataSource(connectionDetails, cdmSchemas[[database]], scratchSchema, database, cohortTable = cohortTable,
                                  covariateSettings = finalCovariateSettings(scratchSchema, cohortTable),
                                  cohortDefinitionSet = benchmarkCohortDefinitionSet(), incremental = TRUE,
                                  combinedCohorts = combinedFeatureCohorts())
base <- benchmarkTasks()
base <- base[base$benchmark == "HFT2DM", ]
counts <- execute(dataSource, base, outputFolder, mode = "counts")
tasks <- expandFeatureSets(addSampleSizes(base, feasibleCounts(counts), budgets = budgets), featureSets)
execute(dataSource, tasks, outputFolder, mode = "models", orderings = 1, methods = fitAlgorithmComparisonMethods(),
        workers = if (exists("workers")) workers else 1)
compareFitAlgorithms(outputFolder, database)
message("Send: ", exportResults(outputFolder, minCellCount = minCellCount))
