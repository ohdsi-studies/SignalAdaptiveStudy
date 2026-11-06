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
#' @param uses Which cohorts: `"study"` (the original tasks and the 4 catalogue tasks), `"benchmark"` (the other PLP
#'   benchmark cohorts) and `"feature"` (the phenotype features, PhenotypeLibrary IDs + 100000), as in
#'   `inst/settings/cohorts.csv`.
#' @return A cohort definition set (`cohortId`, `cohortName`, `json`, `sql`) for CohortGenerator, with `phenotype`
#'   and `role`.
#' @export
studyCohortDefinitionSet <- function(uses = "study") {
  cohorts <- utils::read.csv(system.file("settings", "cohorts.csv", package = "SignalAdaptiveStudy"))
  cohorts <- cohorts[cohorts$use %in% uses, ]
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

#' Sample sizes for the model runs
#'
#' Large populations are sampled at extraction (a person-level sample of the target cohort, before covariates are
#' extracted), so that a task needs neither the full covariates of millions of persons nor predictions for a test set
#' of millions. A task whose population is larger than `cap` is sampled to the larger of `cap` and the size that gives
#' `margin` x the largest budget in expected training-pool events; smaller tasks are not sampled. Persons are sampled
#' with probability `sampleFraction` = `sampleSize` / population, so the sampled population has about `sampleSize`
#' persons. Counts ([execute()] with `mode = "counts"`) always use the full population.
#'
#' @param tasks Tasks, e.g. from [originalCohortTasks()].
#' @param counts The counts of these tasks in one database, from [execute()] with `mode = "counts"`.
#' @param budgets The event budgets of the model runs.
#' @param cap Largest population that is not sampled.
#' @param margin Expected training-pool events per event of the largest budget.
#' @return `tasks` (those with counts) with `sampleSize` (expected population after sampling; `NA`: not sampled),
#'   `sampleFraction`, `expectedPoolEvents`, `budgets` (a list column: the budgets up to the expected training-pool
#'   events) and `eventsPerPerson` (training-pool events per person in the full population).
#' @export
addSampleSizes <- function(tasks, counts, budgets = c(100, 250, 500, 750, 1000), cap = 2e6, margin = 1.5) {
  countColumns <- c("n", "poolN", "poolEvents", "testN", "testEvents")
  tasks <- merge(tasks[, setdiff(names(tasks), countColumns)], counts[, c("taskId", "n", "poolN", "poolEvents")],
                 by = "taskId", sort = FALSE)
  eventsPerPerson <- tasks$poolEvents / tasks$n
  needed <- ceiling(margin * max(budgets) / eventsPerPerson)
  tasks$sampleSize <- ifelse(tasks$n <= cap | needed >= tasks$n, NA, pmax(cap, needed))
  tasks$sampleFraction <- tasks$sampleSize / tasks$n
  tasks$expectedPoolEvents <- ifelse(is.na(tasks$sampleSize), tasks$poolEvents, floor(tasks$sampleSize * eventsPerPerson))
  tasks$budgets <- lapply(tasks$expectedPoolEvents, function(e) budgets[budgets <= e])
  tasks$eventsPerPerson <- eventsPerPerson
  tasks[, setdiff(names(tasks), countColumns)]
}

# Person-level sample of a target cohort, in SQL that runs on every OHDSI dialect. A person is in the sample when a
# multiplicative (golden-ratio) hash of the person ID falls below fraction x modulus: deterministic, nested across
# fractions, and spread evenly over the range of person IDs (so over enrolment time). PLP's own sampleSize is not used:
# it orders cohort entries by (person_id * month % 123) * (year * day % 123) and person ID, so a small sample is all
# ties at 0 (e.g. a 3% sample: 58% index dates in 2009, 45% person IDs divisible by 41).
SAMPLE_HASH <- list(modulus = 2147483647, multiplier = 1327217885, seed = 20261031)

isSampled <- function(task) !is.null(task$sampleFraction) && !is.na(task$sampleFraction)

sampleThreshold <- function(fraction) round(fraction * SAMPLE_HASH$modulus)

# The hash in R, exact in double precision (the product is split so that no term exceeds 2^48).
personHash <- function(personId) {
  m <- SAMPLE_HASH$modulus; a <- SAMPLE_HASH$multiplier
  x <- abs(personId) %% m
  ((x * (a %/% 65536)) %% m * 65536 %% m + x * (a %% 65536) + SAMPLE_HASH$seed) %% m
}

# The same hash in SQL Server SQL for SqlRender. The products stay below 2^62, within BIGINT. SqlRender translates
# `CAST(x AS BIGINT) % m` for BigQuery but has no rule for Oracle and IRIS, which get MOD().
personHashSql <- function(column, dbms) {
  mod <- function(x) if (dbms %in% c("oracle", "iris")) sprintf("MOD(%s, %.0f)", x, SAMPLE_HASH$modulus) else
    sprintf("CAST(%s AS BIGINT) %% %.0f", x, SAMPLE_HASH$modulus)
  mod(sprintf("(%s) * %.0f + %.0f", mod(sprintf("ABS(%s)", column)), SAMPLE_HASH$multiplier, SAMPLE_HASH$seed))
}

# Writes the sample of targetId (fraction of its persons) to cohortDatabaseSchema.sampleTable, replacing it.
createTargetSample <- function(connection, cohortDatabaseSchema, cohortTable, sampleTable, targetId, fraction,
                               tempEmulationSchema = NULL) {
  sql <- "DROP TABLE IF EXISTS @cohort_database_schema.@sample_table;
SELECT cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
INTO @cohort_database_schema.@sample_table
FROM @cohort_database_schema.@cohort_table
WHERE cohort_definition_id = @target_id AND @hash < @threshold;"
  DatabaseConnector::renderTranslateExecuteSql(
    connection, sql, cohort_database_schema = cohortDatabaseSchema, sample_table = sampleTable,
    cohort_table = cohortTable, target_id = targetId, threshold = sprintf("%.0f", sampleThreshold(fraction)),
    hash = personHashSql("subject_id", DatabaseConnector::dbms(connection)), tempEmulationSchema = tempEmulationSchema,
    progressBar = FALSE, reportOverallTime = FALSE)
}

# The PLP benchmark's population settings (PLPBenchmarks, BenchmarkModelDesignList): >= 365 days of prior
# observation, no outcome any time before, the outcome in the risk window, cases kept without full follow-up. The
# risk window's start and first-entry rule differ between benchmarks: day 0 for MACE; every entry for the COVID, lung
# cancer and dementia problems. A task without those columns: day 1, first entry.
benchmarkPopulationSettings <- function(timeAtRiskDays, riskWindowStart = 1, firstExposureOnly = TRUE) {
  PatientLevelPrediction::createStudyPopulationSettings(
    binary = TRUE, includeAllOutcomes = TRUE, firstExposureOnly = firstExposureOnly, washoutPeriod = 365,
    removeSubjectsWithPriorOutcome = TRUE, priorOutcomeLookback = 999999, requireTimeAtRisk = FALSE,
    minTimeAtRisk = 1, riskWindowStart = riskWindowStart, startAnchor = "cohort start", riskWindowEnd = timeAtRiskDays,
    endAnchor = "cohort start")
}

# A task's population settings: the original rules, or the benchmark's for the catalogue tasks.
taskPopulationSettings <- function(task) {
  rule <- if (is.null(task$populationRule)) "original" else task$populationRule
  switch(rule,
    original = originalProductionPopulationSettings(task$timeAtRiskDays),
    benchmark = benchmarkPopulationSettings(task$timeAtRiskDays, task$riskWindowStart %||% 1,
                                            as.logical(task$firstExposureOnly %||% TRUE)),
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
#' each task's `plpData` with PLP and caches it in `outputFolder/data`, from a person-level sample of the target
#' cohort (in `cohortTable_sample`) if the task has a `sampleFraction` ([addSampleSizes()]). Counts ([execute()] with `mode = "counts"`) extract the full population's cohorts and
#' outcomes only (no covariates), so feasibility takes minutes. The study split applies: PLP's default 75/25
#' split fixes the test set, each ordering is a seeded random order of the training pool cut at the k-th event, and
#' 3 stratified folds.
#'
#' @param connectionDetails DatabaseConnector connection details (from environment variables or keyring, not code).
#' @param cdmDatabaseSchema Schema of the CDM tables (and vocabulary).
#' @param cohortDatabaseSchema A writable schema for the cohort tables.
#' @param databaseName Name used in the results.
#' @param cohortTable Name of the cohort table (CohortGenerator also creates its statistics tables).
#' @param tempEmulationSchema For platforms without temp tables (e.g. Spark, Oracle).
#' @param covariateSettings FeatureExtraction covariate settings, or a function of a task that returns them.
#' @param cohortDefinitionSet Cohorts to generate.
#' @param combinedCohorts Cohorts made after generation as the union of generated ones: a data frame with `cohortId`,
#'   `cohortName` and `sourceIds` (a list column), as from [combinedFeatureCohorts()]. `NULL`: none.
#' @param incremental Skip cohorts already generated with the same definition (CohortGenerator keeps the record in the
#'   database). `FALSE`: every cohort is generated again in each run.
#' @return A data source for [execute()], with tasks from [studyTasks()].
#' @export
createCdmDataSource <- function(connectionDetails, cdmDatabaseSchema, cohortDatabaseSchema, databaseName,
                                cohortTable = "sar_study_cohort", tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                covariateSettings = studyCovariateSettings(),
                                cohortDefinitionSet = studyCohortDefinitionSet(), incremental = FALSE,
                                combinedCohorts = NULL) {
  taskCovariates <- function(task) if (is.function(covariateSettings)) covariateSettings(task) else covariateSettings
  combinedFingerprint <- if (NROW(combinedCohorts)) lapply(seq_len(nrow(combinedCohorts)), function(i) {
    list(cohortId = as.numeric(combinedCohorts$cohortId[i]), sourceIds = as.numeric(combinedCohorts$sourceIds[[i]]))
  })
  state <- new.env()
  state$generated <- FALSE
  generate <- function() {
    if (state$generated) return(invisible())
    tables <- CohortGenerator::getCohortTableNames(cohortTable)
    timed("all", "generate cohorts", {
      CohortGenerator::createCohortTables(connectionDetails = connectionDetails, cohortDatabaseSchema = cohortDatabaseSchema,
                                          cohortTableNames = tables, incremental = TRUE)
      CohortGenerator::generateCohortSet(connectionDetails = connectionDetails, cdmDatabaseSchema = cdmDatabaseSchema,
                                         tempEmulationSchema = tempEmulationSchema,
                                         cohortDatabaseSchema = cohortDatabaseSchema, cohortTableNames = tables,
                                         cohortDefinitionSet = cohortDefinitionSet, incremental = incremental)
      if (NROW(combinedCohorts)) {
        connection <- DatabaseConnector::connect(connectionDetails)
        on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
        createCombinedCohorts(connection, cohortDatabaseSchema, cohortTable, combinedCohorts, tempEmulationSchema)
      }
    })
    state$generated <- TRUE
  }
  sampleTable <- paste0(cohortTable, "_sample")
  extract <- function(task, covariateSettings, sample = FALSE) {
    generate()
    if (sample) {
      timed(task$taskId, "sample target", {
        connection <- DatabaseConnector::connect(connectionDetails)
        on.exit(DatabaseConnector::disconnect(connection))
        createTargetSample(connection, cohortDatabaseSchema, cohortTable, sampleTable, task$targetId,
                           task$sampleFraction, tempEmulationSchema)
      })
    }
    details <- PatientLevelPrediction::createDatabaseDetails(
      connectionDetails = connectionDetails, cdmDatabaseSchema = cdmDatabaseSchema, cdmDatabaseName = databaseName,
      cdmDatabaseId = databaseName, tempEmulationSchema = tempEmulationSchema,
      cohortDatabaseSchema = cohortDatabaseSchema, cohortTable = if (sample) sampleTable else cohortTable,
      outcomeDatabaseSchema = cohortDatabaseSchema, outcomeTable = cohortTable,
      targetId = task$targetId, outcomeIds = task$outcomeId)
    PatientLevelPrediction::getPlpData(details, covariateSettings, PatientLevelPrediction::createRestrictPlpDataSettings())
  }
  sampleLabel <- function(task) {
    if (isSampled(task)) sprintf("person hash %.0f, fraction %.12g", SAMPLE_HASH$seed, task$sampleFraction) else "none"
  }
  cohortSql <- function(cohortId) cohortDefinitionSet$sql[cohortDefinitionSet$cohortId == cohortId]
  # What a task's data depend on (R/Fingerprint.R).
  dataFingerprint <- function(task) {
    list(database = databaseName, cdm = cdmDatabaseSchema, targetId = as.numeric(task$targetId),
         outcomeId = as.numeric(task$outcomeId),
         targetCohort = cohortSql(task$targetId), outcomeCohort = cohortSql(task$outcomeId),
         covariates = taskCovariates(task), combined = combinedFingerprint, sample = sampleLabel(task))
  }
  # The folder of the task's cached plpData, extracted first if it is not there or was made with other settings (its
  # fingerprint). A cache from before fingerprints is checked by its sample only (sample.txt), as before.
  cache <- function(task, outputFolder) {
    folder <- file.path(outputFolder, "data", databaseName, task$taskId)
    fp <- fingerprint(dataFingerprint(task))
    changes <- fingerprintChanges(fp, file.path(folder, FINGERPRINT_FILE))
    current <- if (is.null(changes)) {
      sampleFile <- file.path(folder, "sample.txt")
      cached <- if (file.exists(sampleFile)) readLines(sampleFile) else if (file.exists(file.path(folder, "covariates"))) "none" else ""
      if (identical(cached, sampleLabel(task))) {
        ParallelLogger::logInfo(sprintf("%s: data cache from before fingerprints, same sample: used as it is", task$taskId))
      }
      identical(cached, sampleLabel(task))
    } else {
      if (length(changes)) ParallelLogger::logInfo(sprintf("%s: data cache made with other settings (%s), extracted again",
                                                           task$taskId, paste(changes, collapse = ", ")))
      !length(changes)
    }
    if (!current) {
      plpData <- timed(task$taskId, "extract", extract(task, taskCovariates(task), sample = isSampled(task)))
      # saved under a temporary name and renamed when complete, so a crash mid-save leaves no half-written cache
      partial <- paste0(folder, "_partial")
      unlink(c(partial, folder), recursive = TRUE)
      dir.create(dirname(folder), recursive = TRUE, showWarnings = FALSE)
      timed(task$taskId, "save plpData", PatientLevelPrediction::savePlpData(plpData, partial))
      writeLines(sampleLabel(task), file.path(partial, "sample.txt"))
      writeFingerprint(fp, file.path(partial, FINGERPRINT_FILE))
      if (!file.rename(partial, folder)) stop("could not rename ", partial, " to ", folder)
    }
    folder
  }
  populationSettings <- taskPopulationSettings
  newDataSource(
    databaseName, "study",
    loadPlpData = function(task, outputFolder) {
      timed(task$taskId, "load plpData", PatientLevelPrediction::loadPlpData(cache(task, outputFolder)))
    },
    plpDataFolder = function(task, outputFolder) cache(task, outputFolder),
    dataFingerprint = dataFingerprint,
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

# A cohort that is the union of others already in the cohort table: every episode of any of them, once. Made again on
# each run (cheap), so it always follows its sources.
createCombinedCohorts <- function(connection, cohortDatabaseSchema, cohortTable, combinedCohorts, tempEmulationSchema) {
  for (i in seq_len(nrow(combinedCohorts))) {
    DatabaseConnector::renderTranslateExecuteSql(connection, "
      DELETE FROM @schema.@cohort_table WHERE cohort_definition_id = @cohort_id;
      INSERT INTO @schema.@cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
      SELECT DISTINCT CAST(@cohort_id AS BIGINT), subject_id, cohort_start_date, cohort_end_date
      FROM @schema.@cohort_table WHERE cohort_definition_id IN (@source_ids);",
      schema = cohortDatabaseSchema, cohort_table = cohortTable, cohort_id = combinedCohorts$cohortId[i],
      source_ids = as.numeric(combinedCohorts$sourceIds[[i]]), tempEmulationSchema = tempEmulationSchema,
      progressBar = FALSE, reportOverallTime = FALSE)
  }
}
