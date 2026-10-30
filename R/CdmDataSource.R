# The production data source: a CDM database through DatabaseConnector. Cohorts come from the Circe definitions in
# inst/cohorts (generated with CohortGenerator from the shipped SQL, so no Java or CirceR is needed), covariates from
# FeatureExtraction, and the population rules from PLP's population settings. The study split applies.

#' The original tasks
#'
#' The 26 prediction tasks of the original SAR run (8 first defined on IPCI, 18 on CCAE): target and outcome
#' phenotypes, time at risk and outcome family.
#'
#' @return A data frame, one row per task.
#' @noRd
originalTasks <- function() {
  utils::read.csv(system.file("settings", "originalTasks.csv", package = "SignalAdaptiveStudy"), stringsAsFactors = FALSE)
}

#' Cohort definitions of the study
#'
#' Circe definitions shipped with the package: for the original tasks, one target cohort per target phenotype and one
#' outcome cohort per outcome phenotype, from the original concept sets (IDs 1001-2010,
#' `extras/production/buildOriginalCohorts.R`); for the new catalogue tasks, the PLP benchmark definitions (ID 3000
#' + the benchmark's cohort ID, `extras/production/addBenchmarkCohorts.R`).
#'
#' @return A cohort definition set (`cohortId`, `cohortName`, `json`, `sql`) for CohortGenerator, with `phenotype`
#'   and `role`.
#' @export
studyCohortDefinitionSet <- function() {
  cohorts <- utils::read.csv(system.file("settings", "cohorts.csv", package = "SignalAdaptiveStudy"))
  read <- function(id, ext) {
    paste(readLines(system.file("cohorts", paste0(id, ext), package = "SignalAdaptiveStudy")), collapse = "\n")
  }
  cohorts$json <- vapply(cohorts$cohortId, read, "", ext = ".json")
  cohorts$sql <- vapply(cohorts$cohortId, read, "", ext = ".sql")
  cohorts
}

#' Original tasks with production cohort IDs
#'
#' The original tasks with their cohort IDs from [studyCohortDefinitionSet()], and the original population rules.
#'
#' @param database Keep the tasks of this original database ("IPCI", "CCAE"), or `NULL` for all.
#' @return A data frame of tasks for [execute()] with a [createCdmDataSource()].
#' @export
originalCohortTasks <- function(database = NULL) {
  tasks <- originalTasks()
  if (!is.null(database)) tasks <- tasks[tasks$database %in% database, ]
  cohorts <- studyCohortDefinitionSet()
  cohorts <- cohorts[cohorts$cohortId < 3000, ]
  id <- function(phenotype, role) cohorts$cohortId[match(paste(phenotype, role), paste(cohorts$phenotype, cohorts$role))]
  tasks$targetId <- id(tasks$target, "target")
  tasks$outcomeId <- id(tasks$outcome, "outcome")
  if (anyNA(tasks$targetId) || anyNA(tasks$outcomeId)) stop("a task has no cohort definition")
  tasks$populationRule <- "original"
  tasks
}

#' New catalogue tasks
#'
#' Four tasks of the SAR task catalogue (SAR-017, SAR-175, SAR-287, SAR-288), defined as the PLP benchmark problems
#' they adapt (ISAF, CCA, BMDD, HFT2DM): the benchmark's target and outcome cohorts and population settings. For
#' feasibility counts only, so that they stay unseen for a later confirmatory study.
#'
#' @return A data frame of tasks for [execute()] with a [createCdmDataSource()].
#' @export
catalogueCohortTasks <- function() {
  utils::read.csv(system.file("settings", "catalogueTasks.csv", package = "SignalAdaptiveStudy"), stringsAsFactors = FALSE)
}

#' All study tasks
#'
#' [originalCohortTasks()] and [catalogueCohortTasks()] with their common columns.
#'
#' @return A data frame of tasks for [execute()] with a [createCdmDataSource()].
#' @export
studyTasks <- function() {
  columns <- c("taskId", "target", "outcome", "timeAtRiskDays", "targetId", "outcomeId", "populationRule")
  rbind(originalCohortTasks()[, columns], catalogueCohortTasks()[, columns])
}

#' Study covariates
#'
#' Binary indicators of conditions, drugs, procedures, devices, measurements and observations in the 365 days up to
#' and including index (as the original features), sex, and age in years.
#'
#' @return FeatureExtraction covariate settings.
#' @export
studyCovariateSettings <- function() {
  FeatureExtraction::createCovariateSettings(
    useDemographicsGender = TRUE, useDemographicsAge = TRUE,
    useConditionOccurrenceLongTerm = TRUE, useDrugExposureLongTerm = TRUE, useProcedureOccurrenceLongTerm = TRUE,
    useDeviceExposureLongTerm = TRUE, useMeasurementLongTerm = TRUE, useObservationLongTerm = TRUE,
    longTermStartDays = -365, endDays = 0)
}

# The PLP benchmark's population settings (PLPBenchmarks, BenchmarkModelDesignList): >= 365 days of prior
# observation, no outcome any time before, the outcome in days 1 to the horizon, cases kept without full follow-up.
benchmarkPopulationSettings <- function(timeAtRiskDays) {
  PatientLevelPrediction::createStudyPopulationSettings(
    binary = TRUE, includeAllOutcomes = TRUE, firstExposureOnly = TRUE, washoutPeriod = 365,
    removeSubjectsWithPriorOutcome = TRUE, priorOutcomeLookback = 999999, requireTimeAtRisk = FALSE,
    minTimeAtRisk = 1, riskWindowStart = 1, startAnchor = "cohort start", riskWindowEnd = timeAtRiskDays,
    endAnchor = "cohort start")
}

# A task's population settings: the original rules, or the benchmark's for the catalogue tasks.
taskPopulationSettings <- function(task) {
  rule <- if (is.null(task$populationRule)) "original" else task$populationRule
  switch(rule,
    original = originalProductionPopulationSettings(task$timeAtRiskDays),
    benchmark = benchmarkPopulationSettings(task$timeAtRiskDays),
    stop("unknown population rule: ", rule))
}

# The original task rules as PLP population settings: no outcome in days -365 to 0, the outcome in days 1 to the
# horizon, and follow-up to at least the horizon for everyone (cases too). PLP removes prior outcomes with
# daysToEvent > -priorOutcomeLookback, and keeps a person when min(horizon, daysToObsEnd) >= 1 + minTimeAtRisk.
originalProductionPopulationSettings <- function(timeAtRiskDays) {
  PatientLevelPrediction::createStudyPopulationSettings(
    binary = TRUE, includeAllOutcomes = FALSE, firstExposureOnly = TRUE, washoutPeriod = 0,
    removeSubjectsWithPriorOutcome = TRUE, priorOutcomeLookback = 366, requireTimeAtRisk = TRUE,
    minTimeAtRisk = timeAtRiskDays - 1, riskWindowStart = 1, startAnchor = "cohort start",
    riskWindowEnd = timeAtRiskDays, endAnchor = "cohort start")
}

#' Data source: a CDM database
#'
#' Generates the cohorts of [studyCohortDefinitionSet()] once (in `cohortDatabaseSchema.cohortTable`), extracts
#' each task's `plpData` with PLP and caches it in `outputFolder/data`. Counts ([execute()] with `mode = "counts"`)
#' extract the cohorts and outcomes only (no covariates), so feasibility takes minutes. The study split applies: PLP's default 75/25
#' split fixes the test set, each ordering is a seeded random order of the training pool cut at the k-th event, and
#' 3 stratified folds.
#'
#' @param connectionDetails DatabaseConnector connection details (from environment variables or keyring, not code).
#' @param cdmDatabaseSchema Schema of the CDM tables (and vocabulary).
#' @param cohortDatabaseSchema A writable schema for the cohort tables.
#' @param databaseName Name used in the results.
#' @param cohortTable Name of the cohort table (CohortGenerator also creates its statistics tables).
#' @param tempEmulationSchema For platforms without temp tables (e.g. Spark, Oracle).
#' @param covariateSettings FeatureExtraction covariate settings.
#' @param cohortDefinitionSet Cohorts to generate.
#' @return A data source for [execute()], with tasks from [studyTasks()].
#' @export
createCdmDataSource <- function(connectionDetails, cdmDatabaseSchema, cohortDatabaseSchema, databaseName,
                                cohortTable = "sar_study_cohort", tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                covariateSettings = studyCovariateSettings(),
                                cohortDefinitionSet = studyCohortDefinitionSet()) {
  state <- new.env()
  state$generated <- FALSE
  generate <- function() {
    if (state$generated) return(invisible())
    tables <- CohortGenerator::getCohortTableNames(cohortTable)
    CohortGenerator::createCohortTables(connectionDetails = connectionDetails, cohortDatabaseSchema = cohortDatabaseSchema,
                                        cohortTableNames = tables, incremental = TRUE)
    CohortGenerator::generateCohortSet(connectionDetails = connectionDetails, cdmDatabaseSchema = cdmDatabaseSchema,
                                       tempEmulationSchema = tempEmulationSchema,
                                       cohortDatabaseSchema = cohortDatabaseSchema, cohortTableNames = tables,
                                       cohortDefinitionSet = cohortDefinitionSet, incremental = FALSE)
    state$generated <- TRUE
  }
  extract <- function(task, covariateSettings) {
    generate()
    details <- PatientLevelPrediction::createDatabaseDetails(
      connectionDetails = connectionDetails, cdmDatabaseSchema = cdmDatabaseSchema, cdmDatabaseName = databaseName,
      cdmDatabaseId = databaseName, tempEmulationSchema = tempEmulationSchema,
      cohortDatabaseSchema = cohortDatabaseSchema, cohortTable = cohortTable,
      outcomeDatabaseSchema = cohortDatabaseSchema, outcomeTable = cohortTable,
      targetId = task$targetId, outcomeIds = task$outcomeId)
    PatientLevelPrediction::getPlpData(details, covariateSettings, PatientLevelPrediction::createRestrictPlpDataSettings())
  }
  populationSettings <- taskPopulationSettings
  newDataSource(
    databaseName, "study",
    loadPlpData = function(task, outputFolder) {
      folder <- file.path(outputFolder, "data", databaseName, task$taskId)
      if (!file.exists(file.path(folder, "covariates"))) {
        plpData <- extract(task, covariateSettings)
        dir.create(dirname(folder), recursive = TRUE, showWarnings = FALSE)
        PatientLevelPrediction::savePlpData(plpData, folder)
      }
      PatientLevelPrediction::loadPlpData(folder)
    },
    populationSettings = populationSettings,
    split = studySplit,
    counts = function(task) {
      # PLP needs some covariate: sex, the cheapest
      plpData <- extract(task, FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE))
      population <- PatientLevelPrediction::createStudyPopulation(plpData, task$outcomeId, populationSettings(task))
      s <- studySplit(population, task, NULL, 1)
      list(n = nrow(population), poolN = attr(s, "poolN"), poolEvents = attr(s, "poolEvents"),
           testN = attr(s, "testN"), testEvents = attr(s, "testEvents"))
    }
  )
}
