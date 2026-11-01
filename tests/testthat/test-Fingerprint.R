# Settings fingerprints (R/Fingerprint.R): finished runs made with other settings stop execute(); counts made with other
# settings are counted again; a data cache made with other settings is extracted again; outputs from before
# fingerprints are accepted.

simulatedSource <- function(riskWindowEnd = 365) {
  profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
  source <- PatientLevelPrediction::simulatePlpData(profile, n = 2000, seed = 11)
  prevalence <- source$covariateData$covariates %>% dplyr::count(.data$covariateId) %>% dplyr::collect()
  signal <- utils::head(prevalence$covariateId[prevalence$n > 0.05 * 2000 & prevalence$n < 0.4 * 2000 &
                                                 prevalence$covariateId != 1002], 6)
  coefficients <- c("(Intercept)" = -3, stats::setNames(rep(c(1.5, -1.2, 1.0), length.out = length(signal)), signal))
  PatientLevelPrediction::simulatePlpBenchmarkData(source, coefficients, n = 8000, riskWindowStart = 1,
                                                   riskWindowEnd = 365, outcomeId = 3, targetOutcomeRate = 0.06,
                                                   seed = 12, returnTruth = FALSE)
}

test_that("finished runs made with other settings stop execute(); counts are recounted; old outputs are accepted", {
  skip_on_cran()
  plpData <- simulatedSource()
  settings <- function(end) PatientLevelPrediction::createStudyPopulationSettings(riskWindowStart = 1,
                                                                                    riskWindowEnd = end, requireTimeAtRisk = FALSE)
  dataSource <- createPlpDataSource(plpData, "SIM", settings(365))
  tasks <- data.frame(taskId = "sim_task", outcomeId = 3, timeAtRiskDays = 365)
  outputFolder <- tempfile("sar_fingerprint_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  lasso <- list(LASSO = studyMethods()$LASSO)
  run <- file.path(outputFolder, "models", "SIM", "sim_task", "k100_o1_LASSO")

  execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1, methods = lasso)
  expect_true(file.exists(file.path(run, "fingerprint.json")))
  stamp <- file.mtime(file.path(run, "result.csv"))
  execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1, methods = lasso) # same: skipped
  expect_identical(file.mtime(file.path(run, "result.csv")), stamp)

  # other method settings under the same name, and another population rule: execute() stops and says what changed
  changed <- lasso; changed$LASSO$param$interceptOnly <- FALSE
  expect_error(execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1, methods = changed),
               "k100_o1_LASSO: method")
  expect_error(execute(createPlpDataSource(plpData, "SIM", settings(730)), tasks, outputFolder, mode = "models",
                       budgets = 100, orderings = 1, methods = lasso), "k100_o1_LASSO: population")
  expect_identical(file.mtime(file.path(run, "result.csv")), stamp)
  # a run from before fingerprints is accepted as it is
  unlink(file.path(run, "fingerprint.json"))
  execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1, methods = changed)
  expect_identical(file.mtime(file.path(run, "result.csv")), stamp)
  expect_match(readLines(file.path(outputFolder, "log_SIM.txt")), "1 finished run\\(s\\) from before fingerprints", all = FALSE)

  # counts: the same settings are read back, other settings are counted again
  first <- execute(dataSource, tasks, outputFolder, mode = "counts")
  countsFile <- file.path(outputFolder, "counts", "SIM", "sim_task.csv")
  expect_true(file.exists(file.path(outputFolder, "counts", "SIM", "sim_task_fingerprint.json")))
  stamp <- file.mtime(countsFile)
  Sys.sleep(1.1)
  execute(dataSource, tasks, outputFolder, mode = "counts")
  expect_identical(file.mtime(countsFile), stamp)
  again <- execute(createPlpDataSource(plpData, "SIM", settings(180)), tasks, outputFolder, mode = "counts")
  expect_false(identical(file.mtime(countsFile), stamp))
  expect_lt(again$poolEvents, first$poolEvents) # a shorter time at risk: fewer events
})

test_that("a CDM data cache is reused only with the same settings; one from before fingerprints by its sample", {
  skip_on_cran()
  plpData <- simulatedSource()
  outputFolder <- tempfile("sar_cache_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  task <- originalCohortTasks()[1, ]; task$sampleFraction <- NA
  folder <- file.path(outputFolder, "data", "CDM", task$taskId)
  # no connection: any extraction fails, so a returned folder means the cache was reused
  source <- function(covariateSettings = studyCovariateSettings()) {
    createCdmDataSource(connectionDetails = NULL, cdmDatabaseSchema = "cdm", cohortDatabaseSchema = "scratch",
                        databaseName = "CDM", cohortTable = "x", covariateSettings = covariateSettings)
  }
  PatientLevelPrediction::savePlpData(plpData, folder)
  writeLines("none", file.path(folder, "sample.txt"))
  expect_identical(source()$plpDataFolder(task, outputFolder), folder)      # from before fingerprints, same sample
  sampled <- task; sampled$sampleFraction <- 0.5
  expect_error(source()$plpDataFolder(sampled, outputFolder))              # another sample: extracted again
  writeFingerprint(fingerprint(source()$dataFingerprint(task)), file.path(folder, "fingerprint.json"))
  expect_identical(source()$plpDataFolder(task, outputFolder), folder)      # same settings
  fewer <- FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE, useDemographicsAge = TRUE)
  expect_error(source(fewer)$plpDataFolder(task, outputFolder))            # other covariates: extracted again
  expect_identical(fingerprintChanges(fingerprint(source(fewer)$dataFingerprint(task)), file.path(folder, "fingerprint.json")),
                   "covariates")
})
