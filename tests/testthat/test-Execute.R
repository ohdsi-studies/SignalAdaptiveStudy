# Test 5 (PLAN.md): execute() end to end on simulated data, all four methods, one budget; restart and export.
# Only the study API: createPlpDataSource(), execute(), studyMethods(), collectResults(), exportResults().

test_that("execute() runs end to end on simulated data", {
  profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
  source <- PatientLevelPrediction::simulatePlpData(profile, n = 2000, seed = 11)
  prevalence <- source$covariateData$covariates %>% dplyr::count(.data$covariateId) %>% dplyr::collect()
  binary <- prevalence$covariateId[prevalence$n > 0.05 * 2000 & prevalence$n < 0.4 * 2000 &
                                     prevalence$covariateId != 1002] # not age
  signal <- utils::head(binary, 9)
  expect_gte(length(signal), 3)
  coefficients <- c("(Intercept)" = -3, stats::setNames(rep(c(1.5, -1.2, 1.0), length.out = length(signal)), signal))
  plpData <- PatientLevelPrediction::simulatePlpBenchmarkData(source, coefficients, n = 12000, riskWindowStart = 1,
                                                              riskWindowEnd = 365, outcomeId = 3,
                                                              targetOutcomeRate = 0.06, seed = 12, returnTruth = FALSE)
  dataSource <- createPlpDataSource(plpData, "SIM", PatientLevelPrediction::createStudyPopulationSettings(
    riskWindowStart = 1, riskWindowEnd = 365, requireTimeAtRisk = FALSE))
  tasks <- data.frame(taskId = "sim_task", outcomeId = 3, timeAtRiskDays = 365)
  outputFolder <- tempfile("sar_execute_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)

  counts <- execute(dataSource, tasks, outputFolder, mode = "counts")
  expect_identical(nrow(counts), 1L)
  expect_true(counts$feasible)
  expect_identical(counts$n, counts$poolN + counts$testN)
  # a failing task (no such outcome) leaves error.txt and does not stop the others; done tasks are not recounted
  stamp <- file.mtime(file.path(outputFolder, "counts", "SIM", "sim_task.csv"))
  both <- execute(dataSource, rbind(data.frame(taskId = "no_outcome", outcomeId = 99, timeAtRiskDays = 365), tasks),
                  outputFolder, mode = "counts")
  expect_identical(both$taskId, "sim_task")
  expect_true(file.exists(file.path(outputFolder, "counts", "SIM", "no_outcome_error.txt")))
  expect_identical(file.mtime(file.path(outputFolder, "counts", "SIM", "sim_task.csv")), stamp)

  results <- execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1)
  expect_setequal(results$method, c("SAR", "SAL", "Ridge", "LASSO"))
  expect_true(all(results$eventsTrain == 100))
  expect_identical(length(unique(results$nTrain)), 1L)
  expect_true(all(results$nTest == counts$testN) && all(results$eventsTest == counts$testEvents))
  expect_true(all(results$auc >= 0.5 & results$auc <= 1)) # 0.5: a penalty that zeroes every coefficient is a valid choice
  expect_gt(results$auc[results$method == "SAR"], 0.6)
  expect_true(all(is.finite(results$logLoss) & is.finite(results$nullLogLoss) & is.finite(results$brier)))
  comparators <- results[results$method %in% c("Ridge", "LASSO") & !results$interceptOnly, ]
  expect_true(all(comparators$c1 == 0 & comparators$c2 == 0 & comparators$c3 == 0))
  expect_true(all(results$search == "auto" & results$split == "study"))
  expect_true(all(grepl("Cyclops 4.0.0", results$versions)))
  expect_false(any(file.exists(list.files(file.path(outputFolder, "models"), "error.txt", recursive = TRUE, full.names = TRUE))))

  # restart: nothing is refitted
  stamps <- file.mtime(list.files(file.path(outputFolder, "models"), "result.csv", recursive = TRUE, full.names = TRUE))
  again <- execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1)
  expect_identical(nrow(again), 4L)
  expect_identical(file.mtime(list.files(file.path(outputFolder, "models"), "result.csv", recursive = TRUE, full.names = TRUE)), stamps)

  # export: aggregates only, small counts censored
  zipFile <- exportResults(outputFolder, minCellCount = 200)
  expect_setequal(utils::unzip(zipFile, list = TRUE)$Name, c("results.csv", "counts.csv", "export.txt", "steps_SIM.csv"))
  exported <- utils::read.csv(unz(zipFile, "results.csv"))
  expect_true(all(exported$eventsTrain == -200))
  expect_identical(exported$auc, collectResults(outputFolder)$auc)

  # a failed fit records its error
  broken <- setSar(); broken$param$search <- "none"
  execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1, methods = list(Broken = broken))
  error <- readLines(file.path(outputFolder, "models", "SIM", "sim_task", "k100_o1_Broken", "error.txt"))
  expect_match(error, "unknown search: none", all = FALSE)

  # another PLP model runs through the same harness: the SAR fields are NA, the non-zeros come from its coefficients
  other <- execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1,
                   methods = list(PlpLasso = PatientLevelPrediction::setLassoLogisticRegression(seed = 1)))
  plpLasso <- other[other$method == "PlpLasso", ]
  expect_identical(nrow(plpLasso), 1L)
  expect_true(is.na(plpLasso$c1) && is.na(plpLasso$search) && is.na(plpLasso$lambda))
  expect_true(plpLasso$nonZero >= 0 && is.finite(plpLasso$fitSeconds) && plpLasso$auc > 0.5)
  if (requireNamespace("BrokenAdaptiveRidge", quietly = TRUE)) { # its variance is one per coefficient: one row still
    bar <- suppressWarnings(execute( # BAR notes that it excludes the intercept from regularization
      dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1,
      methods = list(BAR = PatientLevelPrediction::setBrokenAdaptiveRidge(seed = 1, threads = 1))))
    expect_identical(sum(bar$method == "BAR"), 1L)
    expect_true(is.na(bar$variance[bar$method == "BAR"]))
  }
})

test_that("runs on the shared split and preprocessing give exactly what runPlp() gives, for every method in turn", {
  profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
  source <- PatientLevelPrediction::simulatePlpData(profile, n = 2000, seed = 11)
  prevalence <- source$covariateData$covariates %>% dplyr::count(.data$covariateId) %>% dplyr::collect()
  signal <- utils::head(prevalence$covariateId[prevalence$n > 100 & prevalence$n < 800 & prevalence$covariateId != 1002], 9)
  coefficients <- c("(Intercept)" = -3, stats::setNames(rep(c(1.5, -1.2, 1.0), length.out = length(signal)), signal))
  plpData <- PatientLevelPrediction::simulatePlpBenchmarkData(source, coefficients, n = 12000, riskWindowStart = 1,
                                                              riskWindowEnd = 365, outcomeId = 3,
                                                              targetOutcomeRate = 0.06, seed = 12, returnTruth = FALSE)
  populationSettings <- PatientLevelPrediction::createStudyPopulationSettings(riskWindowStart = 1, riskWindowEnd = 365,
                                                                              requireTimeAtRisk = FALSE)
  dataSource <- createPlpDataSource(plpData, "SIM", populationSettings)
  task <- data.frame(taskId = "sim_task", outcomeId = 3, timeAtRiskDays = 365)
  population <- PatientLevelPrediction::createStudyPopulation(plpData, 3, populationSettings)
  split <- studySplit(population, task, 100, 1)
  data <- prepareRunData(plpData, population, split)
  folder <- tempfile("sar_identity_")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  for (method in names(studyMethods())) {
    runFolder <- file.path(folder, method)
    dir.create(runFolder, recursive = TRUE)
    row <- runOne(data, population, task, split, studyMethods()[[method]], method, runFolder, dataSource, 100, 1)
    shared <- PatientLevelPrediction::loadPlpModel(file.path(runFolder, "model"))
    reference <- PatientLevelPrediction::runPlp(
      plpData = plpData, outcomeId = 3, analysisId = method, populationSettings = populationSettings,
      splitSettings = PatientLevelPrediction::createExistingSplitSettings(data.frame(rowId = split$rowId, index = split$index)),
      preprocessSettings = studyPreprocessSettings(), modelSettings = studyMethods()[[method]],
      executeSettings = PatientLevelPrediction::createExecuteSettings(
        runSplitData = TRUE, runSampleData = FALSE, runFeatureEngineering = FALSE, runPreprocessData = TRUE,
        runModelDevelopment = TRUE, runCovariateSummary = FALSE),
      saveDirectory = file.path(folder, "reference"))
    test <- reference$prediction[reference$prediction$evaluationType == "Test", ]
    m <- testMetrics(as.numeric(test$outcomeCount > 0), test$value,
                     as.numeric(population$outcomeCount[match(split$rowId[split$index > 0], population$rowId)] > 0))
    expect_identical(row$auc, m$auc, info = method)
    expect_identical(row$logLoss, m$logLoss, info = method)
    expect_identical(row$nTest, nrow(test), info = method)
    PatientLevelPrediction::savePlpModel(reference$model, file.path(folder, "reference", method, "model"))
    saved <- PatientLevelPrediction::loadPlpModel(file.path(folder, "reference", method, "model")) # both through JSON
    expect_identical(shared$covariateImportance$covariateValue, saved$covariateImportance$covariateValue, info = method)
  }
})

test_that("the test AUC holds for large test sets (no integer overflow)", {
  # 8,000 events and 492,000 non-events (8,000 x 492,000 > 2^31); events score 1, 3/4 of non-events 0 and 1/4 2
  y <- c(rep(1, 8000), rep(0, 492000))
  p <- c(rep(0.5, 8000), rep(0.1, 369000), rep(0.9, 123000))
  expect_identical(aucMannWhitney(y, p), 0.75)
  expect_identical(testMetrics(y, p, y)$auc, 0.75)
})
