# Aggregates of one training set and of every Cyclops fit of chosen models on it, for building simulated designs that behave like
# the real data (for the Cyclops work on faster cross-validation). Runs on the data of an earlier fit-algorithm comparison
# (compareFitAlgorithms.R), so it needs no new extraction; the counts step is skipped too when its counts are there.
# Writes a zip of aggregates to send back: design.csv (rows, events, non-zeros per row and column, duplicate and near-duplicate
# columns), column_types.csv (per FeatureExtraction analysis), trace_<model>.csv (one row per Cyclops fit: rows, columns,
# prior variance, seconds, iterations, Newton steps) and trace_summary.csv. Nothing per patient. By default SAR, ridge, LASSO
# and elastic net with the study's algorithm (auto; ccd's fold fits are the same), at k1000 on the default covariates.
#
#   Rscript extras/production/fitTrace.R CCAE
#   Rscript extras/production/fitTrace.R CCAE --budget=1000 --featureSet=defaults --methods=SAR_ccd,SAR_auto,Ridge_ccd,Ridge_auto
#   Rscript extras/production/fitTrace.R CCAE --minFraction=0   # all of it with another frequency filter than the study's (own folder and zip)
#   Rscript extras/production/fitTrace.R CCAE --methods=none        # design only
#
# One training set of the benchmark run (runBenchmarks.R's data, split and preprocessing: the same fits), e.g. to reproduce a
# fold fit that did not converge; --retry fits each such fit again with more iterations, from zero and by coordinate
# descent, and writes Cyclops's convergence details (trace_<model>_failed.txt):
#   Rscript extras/production/fitTrace.R CCAE --data=benchmarks --task=ISAF --featureSet=phenotypes --budget=1000 --ordering=3 --methods=SAL --retry
# Methods: those of the fit-algorithm comparison (SAR_auto, ...) or of the final study (SAR, SAL, LASSO, ...).
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
benchmark <- option("task", "HFT2DM")
ordering <- as.integer(option("ordering", "1"))
data <- option("data", "comparison") # "benchmarks": runBenchmarks.R's output folder
retry <- "--retry" %in% args
methods <- c(fitAlgorithmComparisonMethods(), finalStudyMethods())
chosen <- option("methods", c("SAR_auto", "Ridge_auto", "LASSO_auto", "ElasticNet_auto"))
if (identical(chosen, "none")) chosen <- character() # design only
minFraction <- option("minFraction", NULL) # NULL: the study's preprocessing
if (!is.null(minFraction)) minFraction <- as.numeric(minFraction)
stopifnot(length(budget) == 1, budget %in% finalBudgets(), featureSet %in% c("defaults", "phenotypes"), all(chosen %in% names(methods)),
          length(benchmark) == 1, benchmark %in% benchmarkTasks()$benchmark, length(ordering) == 1, ordering %in% 1:3,
          data %in% c("comparison", "benchmarks"), !"PLP_LASSO" %in% chosen) # PLP's own LASSO is not traced
message("Cyclops ", utils::packageVersion("Cyclops"), " commit ", utils::packageDescription("Cyclops")$RemoteSha %||% "",
        "; ", benchmark, " k", budget, " ordering ", ordering, ", ", featureSet, "; ", paste(chosen, collapse = ", "),
        if (retry) "; failed fits retried")

task <- paste0("benchmark_", benchmark, "_", featureSet)
if (data == "benchmarks") {
  outputFolder <- file.path(getwd(), "SignalAdaptiveStudy_benchmarks", database)
  if (!dir.exists(file.path(outputFolder, "data", database, task))) stop("no benchmark data for ", task, " in ", outputFolder)
} else {
  # The data and counts of the most recent fit-algorithm comparison that has this task's data
  folders <- Sys.glob(file.path(getwd(), "SignalAdaptiveStudy_fit_comparison*", database))
  folders <- folders[dir.exists(file.path(folders, "data", database, task))]
  if (!length(folders)) stop("no fit-algorithm comparison with data for ", task, "; run compareFitAlgorithms.R first")
  outputFolder <- folders[which.max(file.mtime(file.path(folders, "data", database, task)))]
}
message("Data from: ", outputFolder)

cohortTable <- paste0("sar_benchmark_cohort_", tolower(database))
dataSource <- createCdmDataSource(connectionDetails, cdmSchemas[[database]], scratchSchema, database, cohortTable = cohortTable,
                                  covariateSettings = finalCovariateSettings(scratchSchema, cohortTable),
                                  cohortDefinitionSet = benchmarkCohortDefinitionSet(), incremental = TRUE,
                                  combinedCohorts = combinedFeatureCohorts())
# the tasks as the run that made the data builds them (the comparison: HFT2DM alone; the benchmarks: all), so the sample,
# split, orderings and folds are the same; the counts are read from the output folder
base <- benchmarkTasks()
if (data == "comparison") base <- base[base$benchmark == benchmark, ]
counts <- execute(dataSource, base, outputFolder, mode = "counts")
tasks <- expandFeatureSets(addSampleSizes(base, feasibleCounts(counts), budgets = finalBudgets()), featureSet)
tasks <- tasks[tasks$taskId == task, ]
if (nrow(tasks) != 1) stop(task, " is not feasible in ", database)
label <- if (data == "benchmarks" || benchmark != "HFT2DM" || ordering != 1) sprintf("%s_%s_k%g_o%d", benchmark, featureSet, budget, ordering)
reportFolder <- file.path(getwd(), paste0("SignalAdaptiveStudy_fit_trace", if (!is.null(label)) paste0("_", label),
                                          if (!is.null(minFraction)) paste0("_minFraction", minFraction)), database)
zipFile <- SignalAdaptiveStudy:::fitTraceReport(dataSource, tasks, budget, methods[chosen], outputFolder, reportFolder,
                                                ordering = ordering, minCellCount = minCellCount, minFraction = minFraction,
                                                label = label, retry = retry)
message("Send: ", zipFile)
