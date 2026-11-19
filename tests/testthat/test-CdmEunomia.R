# The production path end to end on Eunomia (a small public CDM; downloaded on first use): cohort generation from the
# shipped SQL, light counts for all 26 original tasks (most have too few events there, which must not stop the run),
# and the export.

test_that("counts run end to end on a CDM (Eunomia)", {
  skip_if_not_installed("Eunomia")
  skip_on_cran()
  details <- tryCatch(Eunomia::getEunomiaConnectionDetails(databaseFile = tempfile(fileext = ".duckdb"), dbms = "duckdb"),
                      error = function(e) skip(paste("Eunomia not available:", conditionMessage(e))))
  outputFolder <- tempfile("sar_eunomia_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  dataSource <- createCdmDataSource(details, cdmDatabaseSchema = "main", cohortDatabaseSchema = "main",
                                    databaseName = "EUNOMIA")
  tasks <- studyTasks()
  counts <- suppressWarnings(execute(dataSource, tasks, outputFolder, mode = "counts"))
  errors <- list.files(file.path(outputFolder, "counts", "EUNOMIA"), "_error\\.txt$")
  expect_identical(length(errors) + NROW(counts), nrow(tasks)) # every task is either counted or has its error recorded
  # tasks with no target persons or no outcomes are kept with no events and no split (NA)
  split <- counts[!is.na(counts$poolN), , drop = FALSE]
  if (NROW(split) > 0) expect_true(all(split$n == split$poolN + split$testN))
  expect_true(all(counts$poolEvents[is.na(counts$poolN)] == 0))
  expect_false(any(file.exists(file.path(outputFolder, "data")))) # counts extract no covariates and cache nothing
  zipFile <- exportResults(outputFolder)
  expect_true("counts.csv" %in% utils::unzip(zipFile, list = TRUE)$Name || NROW(counts) == 0)
})

test_that("restricting the extraction to first entries with washout gives the same population and split", {
  skip_if_not_installed("Eunomia")
  skip_on_cran()
  details <- tryCatch(Eunomia::getEunomiaConnectionDetails(databaseFile = tempfile(fileext = ".duckdb"), dbms = "duckdb"),
                      error = function(e) skip(paste("Eunomia not available:", conditionMessage(e))))
  connection <- DatabaseConnector::connect(details)
  # target: every viral sinusitis or pharyngitis start (several per person, some without 365 days before); outcome:
  # acute bronchitis
  DatabaseConnector::executeSql(connection, "
    CREATE TABLE main.restrict_cohort AS
    SELECT DISTINCT CAST(1 AS INTEGER) AS cohort_definition_id, person_id AS subject_id,
           condition_start_date AS cohort_start_date, condition_start_date AS cohort_end_date
    FROM main.condition_occurrence WHERE condition_concept_id IN (40481087, 4112343)
    UNION ALL
    SELECT DISTINCT CAST(2 AS INTEGER), person_id, condition_start_date, condition_start_date
    FROM main.condition_occurrence WHERE condition_concept_id = 260139;", progressBar = FALSE, reportOverallTime = FALSE)
  DatabaseConnector::disconnect(connection)
  population <- benchmarkPopulationSettings(365, 1, TRUE)
  extract <- function(restrict) {
    databaseDetails <- PatientLevelPrediction::createDatabaseDetails(
      connectionDetails = details, cdmDatabaseSchema = "main", cohortDatabaseSchema = "main", cohortTable = "restrict_cohort",
      outcomeDatabaseSchema = "main", outcomeTable = "restrict_cohort", targetId = 1, outcomeIds = 2)
    PatientLevelPrediction::getPlpData(databaseDetails, FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE),
                                       restrict)
  }
  full <- suppressMessages(extract(PatientLevelPrediction::createRestrictPlpDataSettings()))
  restricted <- suppressMessages(extract(extractionRestriction(population)))
  expect_lt(nrow(restricted$cohorts), nrow(full$cohorts) / 2) # many entries were not first, or had no washout
  study <- function(plpData) suppressMessages(PatientLevelPrediction::createStudyPopulation(plpData, 2, population))
  a <- study(full); b <- study(restricted)
  expect_gt(sum(a$outcomeCount > 0), 30)
  columns <- c("subjectId", "cohortStartDate", "outcomeCount", "timeAtRisk", "daysFromObsStart", "ageYear", "gender")
  expect_identical(as.data.frame(a[, intersect(columns, names(a))]), as.data.frame(b[, intersect(columns, names(b))]))
  task <- data.frame(taskId = "t", outcomeId = 2, timeAtRiskDays = 365)
  for (ordering in 1:2) {
    sa <- studySplit(a, task, budget = 20, ordering = ordering); sb <- studySplit(b, task, budget = 20, ordering = ordering)
    expect_identical(a$subjectId[match(sa$rowId, a$rowId)], b$subjectId[match(sb$rowId, b$rowId)])
    expect_identical(sa$index, sb$index)
  }
})

test_that("mode data extracts each task's data once, for the models to use; a database held by another run waits", {
  skip_if_not_installed("Eunomia")
  skip_on_cran()
  details <- tryCatch(Eunomia::getEunomiaConnectionDetails(databaseFile = tempfile(fileext = ".duckdb"), dbms = "duckdb"),
                      error = function(e) skip(paste("Eunomia not available:", conditionMessage(e))))
  outputFolder <- tempfile("sar_eunomia_data_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  # two cohorts of plain SQL: viral sinusitis or pharyngitis starts (several per person) and acute bronchitis
  cohortSql <- function(concepts) sprintf(paste(
    "DELETE FROM @target_database_schema.@target_cohort_table WHERE cohort_definition_id = @target_cohort_id;",
    "INSERT INTO @target_database_schema.@target_cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)",
    "SELECT DISTINCT @target_cohort_id, person_id, condition_start_date, condition_start_date",
    "FROM @cdm_database_schema.condition_occurrence WHERE condition_concept_id IN (%s);"), concepts)
  cohorts <- data.frame(cohortId = c(1, 2), cohortName = c("sinusitis", "bronchitis"), json = "",
                        sql = c(cohortSql("40481087, 4112343"), cohortSql("260139")))
  dataSource <- createCdmDataSource(details, cdmDatabaseSchema = "main", cohortDatabaseSchema = "main", databaseName = "EUNOMIA",
                                    cohortTable = "data_mode_cohort", cohortDefinitionSet = cohorts,
                                    covariateSettings = FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE))
  task <- data.frame(taskId = "sinusitis_bronchitis", targetId = 1, outcomeId = 2, timeAtRiskDays = 365, riskWindowStart = 1,
                     firstExposureOnly = TRUE, populationRule = "benchmark")
  folders <- suppressWarnings(execute(dataSource, task, outputFolder, mode = "data"))
  expect_true(dir.exists(folders[[1]]) && file.exists(file.path(folders[[1]], "fingerprint.json")))
  again <- suppressWarnings(execute(dataSource, task, outputFolder, mode = "data"))
  expect_identical(again, folders)
  steps <- utils::read.csv(file.path(outputFolder, "steps_EUNOMIA.csv"))
  expect_equal(sum(steps$step == "extract"), 1) # the second call found the data
  expect_false(dir.exists(file.path(outputFolder, "busy_EUNOMIA.lock")))
  # extracted with the task's first-entry and washout rules: one entry per person
  cohortsExtracted <- PatientLevelPrediction::loadPlpData(folders[[1]])$cohorts
  expect_gt(nrow(cohortsExtracted), 100)
  expect_false(anyDuplicated(cohortsExtracted$subjectId) > 0)
  # another run holds the database: wait (here: stop after the timeout)
  dir.create(file.path(outputFolder, "busy_EUNOMIA.lock"))
  old <- options(sarLockTimeoutSeconds = 0.5)
  on.exit(options(old), add = TRUE)
  expect_error(execute(dataSource, task, outputFolder, mode = "data"), "held by another run")
})
