# Settings of the final study: the PLP benchmark tasks, the two feature sets, the methods. See README_PRODUCTION.md.

#' PLP benchmark tasks
#'
#' The 13 problems of PLPBenchmarks 2.1.1 (commit 81489d6) with their cohorts (ID 3000 + the benchmark's ID), risk windows and first-entry rule.
#'
#' @return A data frame of tasks for [execute()] with a [createCdmDataSource()].
#' @export
benchmarkTasks <- function() {
  utils::read.csv(system.file("settings", "benchmarkTasks.csv", package = "SignalAdaptiveStudy"), stringsAsFactors = FALSE)
}

#' Cohorts of the benchmark run
#'
#' The cohorts of [benchmarkTasks()] and the phenotype features.
#'
#' @return A cohort definition set for CohortGenerator.
#' @export
benchmarkCohortDefinitionSet <- function() {
  all <- studyCohortDefinitionSet(c("study", "benchmark", "feature"))
  tasks <- benchmarkTasks()
  all[all$cohortId %in% c(tasks$targetId, tasks$outcomeId) | all$use == "feature", ]
}

#' One task per feature set
#'
#' @param tasks Tasks, from [benchmarkTasks()].
#' @param featureSets `"defaults"`, `"phenotypes"` or both: see [finalCovariateSettings()].
#' @return The tasks repeated for each feature set, with `featureSet` and the feature set at the end of `taskId`.
#' @export
expandFeatureSets <- function(tasks, featureSets = c("defaults", "phenotypes")) {
  do.call(rbind, lapply(featureSets, function(set) {
    tasks$featureSet <- set
    tasks$taskId <- paste(tasks$taskId, set, sep = "_")
    tasks
  }))
}

# The last day of the covariate window. It ends before the risk window starts: PLPBenchmarks ends at day -1 on all its
# problems; elsewhere a risk window from day 1 keeps FeatureExtraction's default end at day 0.
covariateEndDay <- function(task) {
  benchmark <- identical(task$populationRule, "benchmark")
  if (benchmark || (task$riskWindowStart %||% 1) < 1) -1 else 0
}

# FeatureExtraction's defaults, without its continuous scores, so that every covariate is binary.
defaultFeatureSettings <- function(endDay) {
  defaults <- FeatureExtraction::createDefaultCovariateSettings()
  on <- names(defaults)[vapply(defaults, isTRUE, logical(1))]
  on <- setdiff(on, c("CharlsonIndex", "Dcsi", "Chads2", "Chads2Vasc"))
  do.call(FeatureExtraction::createCovariateSettings,
          c(stats::setNames(as.list(rep(TRUE, length(on))), paste0("use", on)),
            list(longTermStartDays = defaults$longTermStartDays, mediumTermStartDays = defaults$mediumTermStartDays,
                 shortTermStartDays = defaults$shortTermStartDays, endDays = endDay)))
}

#' Combined phenotype features
#'
#' The 14 antibiotic classes of the phenotype features (PhenotypeLibrary 1201-1214) as one indicator: an episode of any of
#' them. For [createCdmDataSource()], which makes it after generating the cohorts.
#'
#' @return A data frame: `cohortId`, `cohortName`, `sourceIds` (a list column).
#' @export
combinedFeatureCohorts <- function() {
  data.frame(cohortId = 199001, cohortName = "Antibiotics (any of 14 classes)", sourceIds = I(list(100000 + 1201:1214)))
}

# The phenotype features as covariates: the 64 phenotypes of the PLP vignette "Constrained Predictors" with the antibiotic
# classes combined (51 indicators).
phenotypeFeatureCohorts <- function() {
  features <- studyCohortDefinitionSet("feature")
  combined <- combinedFeatureCohorts()
  features <- features[!features$cohortId %in% unlist(combined$sourceIds), c("cohortId", "cohortName")]
  rbind(features, data.frame(cohortId = combined$cohortId, cohortName = combined$cohortName))
}

# Gender, 5-year age groups and the phenotype features, each an indicator of a cohort episode overlapping the window,
# covariate ID = cohort ID x 1000 + 456.
phenotypeFeatureSettings <- function(endDay, cohortDatabaseSchema, cohortTable) {
  features <- phenotypeFeatureCohorts()
  list(FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE, useDemographicsAgeGroup = TRUE),
       FeatureExtraction::createCohortBasedCovariateSettings(
         analysisId = 456, covariateCohortDatabaseSchema = cohortDatabaseSchema, covariateCohortTable = cohortTable,
         covariateCohorts = data.frame(cohortId = features$cohortId, cohortName = features$cohortName),
         valueType = "binary", startDay = -365, endDay = endDay))
}

#' Covariates of the final study
#'
#' Two feature sets, all binary and all penalised. `"defaults"`: FeatureExtraction's
#' `createDefaultCovariateSettings()` without Charlson, DCSI, CHADS2 and CHA2DS2-VASc. `"phenotypes"`: gender, 5-year age
#' groups and the 64 phenotypes of the PLP vignette "Constrained Predictors", the 14 antibiotic classes as one (51
#' indicators; the data source must make the combined cohort, [combinedFeatureCohorts()]). The window starts at day -365 (FeatureExtraction's
#' default windows in the first set) and ends before the risk window starts.
#'
#' @param cohortDatabaseSchema,cohortTable Where the phenotype cohorts are generated: as in [createCdmDataSource()].
#' @return A function of a task (with `featureSet`) that gives its covariate settings, for [createCdmDataSource()].
#' @export
finalCovariateSettings <- function(cohortDatabaseSchema, cohortTable) {
  function(task) {
    set <- task$featureSet
    if (is.null(set)) return(NULL) # counts use no covariates
    switch(set,
      defaults = defaultFeatureSettings(covariateEndDay(task)),
      phenotypes = phenotypeFeatureSettings(covariateEndDay(task), cohortDatabaseSchema, cohortTable),
      stop("unknown feature set: ", set))
  }
}

#' Methods of the final study
#'
#' Primary: SAR, LASSO and ridge. Secondary: PLP's default LASSO, elastic net (mixing 0.5), and two reduced versions of SAR, the
#' shape fixed at (0, 1, -1) and without the third term.
#'
#' @param fitAlgorithm As in [setSar()]; PLP's LASSO uses PLP's own settings.
#' @return A named list of model settings.
#' @export
finalStudyMethods <- function(fitAlgorithm = "auto") {
  shapes <- sarShapes()
  list(SAR = setSar(fitAlgorithm = fitAlgorithm),
       LASSO = setSal(shapes = c(0, 0, 0), name = "LASSO", fitAlgorithm = fitAlgorithm),
       Ridge = setSar(shapes = c(0, 0, 0), name = "Ridge", fitAlgorithm = fitAlgorithm),
       PLP_LASSO = PatientLevelPrediction::setLassoLogisticRegression(seed = STUDY_SEEDS$folds, threads = 1),
       ElasticNet = setElasticNet(fitAlgorithm = fitAlgorithm),
       SAR_fixed = setSar(shapes = c(0, 1, -1), name = "SAR_fixed", fitAlgorithm = fitAlgorithm),
       SAR_no_c3 = setSar(shapes = shapes[shapes[, "c3"] == 0, , drop = FALSE], name = "SAR_no_c3", fitAlgorithm = fitAlgorithm))
}

#' Methods for comparing Cyclops's fitting algorithms
#'
#' SAR, SAR with its shape fixed, LASSO, ridge and elastic net, each fitted by coordinate descent (`_ccd`) and by Cyclops's newer
#' automatic algorithm (`_auto`), through the same folds and search.
#'
#' @return A named list of model settings.
#' @export
fitAlgorithmComparisonMethods <- function() {
  one <- function(name, make) stats::setNames(lapply(c("ccd", "auto"), make), paste(name, c("ccd", "auto"), sep = "_"))
  c(one("SAR", function(a) setSar(name = paste0("SAR_", a), fitAlgorithm = a)),
    one("SAR_fixed", function(a) setSar(shapes = c(0, 1, -1), name = paste0("SAR_fixed_", a), fitAlgorithm = a)),
    one("LASSO", function(a) setSal(shapes = c(0, 0, 0), name = paste0("LASSO_", a), fitAlgorithm = a)),
    one("Ridge", function(a) setSar(shapes = c(0, 0, 0), name = paste0("Ridge_", a), fitAlgorithm = a)),
    one("ElasticNet", function(a) setElasticNet(name = paste0("ElasticNet_", a), fitAlgorithm = a)))
}

#' Compare the coefficients of the two fitting algorithms
#'
#' For each pair of runs of [fitAlgorithmComparisonMethods()] on the same task, budget and ordering, from the saved models: aggregates
#' only, no coefficients.
#'
#' @param outputFolder The `outputFolder` of [execute()].
#' @param databaseName The database.
#' @return A data frame, one row per pair; also written to `outputFolder/coefficient_comparison.csv`.
#' @export
compareFitAlgorithms <- function(outputFolder, databaseName) {
  root <- file.path(outputFolder, "models", databaseName)
  folders <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  ccd <- folders[grepl("_ccd$", folders) & dir.exists(file.path(folders, "model"))]
  read <- function(folder) {
    model <- PatientLevelPrediction::loadPlpModel(file.path(folder, "model"))
    co <- as.data.frame(model$model$coefficientsExact, stringsAsFactors = FALSE)
    stats::setNames(as.numeric(co$betas), co$covariateIds)
  }
  rows <- lapply(ccd, function(folder) {
    other <- sub("_ccd$", "_auto", folder)
    if (!dir.exists(file.path(other, "model"))) return(NULL)
    a <- read(folder); b <- read(other)
    ids <- union(names(a), names(b))
    x <- ifelse(ids %in% names(a), a[ids], 0); y <- ifelse(ids %in% names(b), b[ids], 0)
    data.frame(database = databaseName, taskId = basename(dirname(folder)), run = sub("_ccd$", "", basename(folder)),
               columns = length(ids), maxAbsDifference = max(abs(x - y)), rmsDifference = sqrt(mean((x - y)^2)),
               correlation = if (stats::sd(x) > 0 && stats::sd(y) > 0) stats::cor(x, y) else NA_real_,
               nonZeroCcd = sum(x != 0), nonZeroAuto = sum(y != 0), nonZeroInOnlyOne = sum((x != 0) != (y != 0)))
  })
  result <- do.call(rbind, rows)
  if (NROW(result)) utils::write.csv(result, file.path(outputFolder, "coefficient_comparison.csv"), row.names = FALSE)
  result
}

#' Budgets of the final study
#'
#' @return Events in the training set.
#' @export
finalBudgets <- function() c(100, 250, 500, 1000)

#' Tasks a database can run
#'
#' At least 1,000 training-pool events, which k = 1000 needs, and 333 test events.
#'
#' @param counts Counts from [execute()] with `mode = "counts"`.
#' @return The rows of `counts` that pass.
#' @export
feasibleCounts <- function(counts) counts[counts$poolEvents >= 1000 & counts$testEvents >= 333, ]

#' Elastic net model settings
#'
#' Logistic regression with Cyclops's elastic net prior, with mixing parameter 0.5, through the same folds, search and code as LASSO.
#'
#' @param name Model name shown by PLP.
#' @param fitAlgorithm As in [setSar()].
#' @return An object of class `modelSettings`.
#' @export
setElasticNet <- function(name = "ElasticNet", fitAlgorithm = "auto") {
  signalAdaptiveSettings("elasticnet", c(0, 0, 0), interceptOnly = TRUE, name = name, searchSettings = list(fitAlgorithm = fitAlgorithm))
}
