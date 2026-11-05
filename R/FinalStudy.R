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

# Gender, 5-year age groups and the 64 phenotypes of the PLP vignette "Constrained Predictors" (Reps et al. 2025), each an
# indicator of a cohort episode overlapping the window, covariate ID = cohort ID x 1000 + 456.
phenotypeFeatureSettings <- function(endDay, cohortDatabaseSchema, cohortTable) {
  features <- studyCohortDefinitionSet("feature")
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
#' groups and the 64 phenotypes of the PLP vignette "Constrained Predictors". The window starts at day -365 (FeatureExtraction's
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
#' @return A named list of model settings.
#' @export
finalStudyMethods <- function() {
  shapes <- sarShapes()
  list(SAR = setSar(), LASSO = setSal(shapes = c(0, 0, 0), name = "LASSO"), Ridge = setSar(shapes = c(0, 0, 0), name = "Ridge"),
       PLP_LASSO = PatientLevelPrediction::setLassoLogisticRegression(seed = STUDY_SEEDS$folds, threads = 1),
       ElasticNet = setElasticNet(),
       SAR_fixed = setSar(shapes = c(0, 1, -1), name = "SAR_fixed"),
       SAR_no_c3 = setSar(shapes = shapes[shapes[, "c3"] == 0, , drop = FALSE], name = "SAR_no_c3"))
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
#' @return An object of class `modelSettings`.
#' @export
setElasticNet <- function(name = "ElasticNet") {
  signalAdaptiveSettings("elasticnet", c(0, 0, 0), interceptOnly = TRUE, name = name)
}
