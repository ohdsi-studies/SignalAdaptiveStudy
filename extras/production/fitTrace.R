# Aggregates of one training set and of every Cyclops fit of chosen models on it, for building simulated designs that behave like
# the real data (for the Cyclops work on faster cross-validation). Runs on the data of an earlier fit-algorithm comparison
# (compareFitAlgorithms.R), so it needs no new extraction; the counts step is skipped too when its counts are there.
# Writes a zip of aggregates to send back: design.csv (rows, events, non-zeros per row and column, duplicate and near-duplicate
# columns), column_types.csv (per FeatureExtraction analysis), trace_<model>.csv (one row per Cyclops fit: rows, columns,
# prior variance, seconds, iterations, Newton steps) and trace_summary.csv. Nothing per patient. By default SAR, ridge, LASSO
# and elastic net, both algorithms, at k1000 on the default covariates (about 45 minutes; SAR is most of it).
#
#   Rscript extras/production/fitTrace.R CCAE
#   Rscript extras/production/fitTrace.R CCAE --budget=1000 --featureSet=defaults --methods=SAR_ccd,SAR_auto,Ridge_ccd,Ridge_auto
#   Rscript extras/production/fitTrace.R CCAE --minFraction=0.001   # all of it after PLP's frequency filter (its own folder and zip)
#   Rscript extras/production/fitTrace.R CCAE --methods=none        # design only
# (in the project folder; settings as for runBenchmarks.R).

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

source("extras/localSettings.R")   # scratchSchema, minCellCount, cdmSchemas
options(sqlRenderTempEmulationSchema = scratchSchema)
args <- commandArgs(trailingOnly = TRUE)
database <- args[1]
stopifnot(!is.na(database), database %in% names(cdmSchemas))
option <- function(name, default) {
  given <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(given)) strsplit(given, ",")[[1]] else default
}
budget <- as.numeric(option("budget", "1000"))
featureSet <- option("featureSet", "defaults")
methods <- fitAlgorithmComparisonMethods()
chosen <- option("methods", c("SAR_ccd", "SAR_auto", "Ridge_ccd", "Ridge_auto", "LASSO_ccd", "LASSO_auto", "ElasticNet_ccd", "ElasticNet_auto"))
if (identical(chosen, "none")) chosen <- character() # design only
minFraction <- as.numeric(option("minFraction", "0"))
stopifnot(length(budget) == 1, budget %in% finalBudgets(), featureSet %in% c("defaults", "phenotypes"), all(chosen %in% names(methods)))
message("Cyclops ", utils::packageVersion("Cyclops"), " commit ", utils::packageDescription("Cyclops")$RemoteSha %||% "",
        "; k", budget, ", ", featureSet, "; ", paste(chosen, collapse = ", "))

# The data and counts of the most recent fit-algorithm comparison that has this task's data
task <- paste0("benchmark_HFT2DM_", featureSet)
folders <- Sys.glob(file.path(getwd(), "SignalAdaptiveStudy_fit_comparison*", database))
folders <- folders[dir.exists(file.path(folders, "data", database, task))]
if (!length(folders)) stop("no fit-algorithm comparison with data for ", task, "; run compareFitAlgorithms.R first")
outputFolder <- folders[which.max(file.mtime(file.path(folders, "data", database, task)))]
message("Data from: ", outputFolder)

cohortTable <- paste0("sar_benchmark_cohort_", tolower(database))
dataSource <- createCdmDataSource(connectionDetails, cdmSchemas[[database]], scratchSchema, database, cohortTable = cohortTable,
                                  covariateSettings = finalCovariateSettings(scratchSchema, cohortTable),
                                  cohortDefinitionSet = benchmarkCohortDefinitionSet(), incremental = TRUE,
                                  combinedCohorts = combinedFeatureCohorts())
base <- benchmarkTasks()
base <- base[base$benchmark == "HFT2DM", ]
counts <- execute(dataSource, base, outputFolder, mode = "counts")
tasks <- expandFeatureSets(addSampleSizes(base, feasibleCounts(counts), budgets = finalBudgets()), featureSet) # as the comparison: same sample
reportFolder <- file.path(getwd(), paste0("SignalAdaptiveStudy_fit_trace", if (minFraction > 0) paste0("_minFraction", minFraction)), database)
zipFile <- SignalAdaptiveStudy:::fitTraceReport(dataSource, tasks[1, ], budget, methods[chosen], outputFolder, reportFolder,
                                                minCellCount = minCellCount, minFraction = minFraction)
message("Send: ", zipFile)
