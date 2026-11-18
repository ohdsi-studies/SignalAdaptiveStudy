test_that("designSummary counts duplicates, near-duplicates and column types", {
  X <- Matrix::sparseMatrix(i = c(1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 6), j = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5), x = 1,
                            dims = c(10, 5))
  covariateData <- Andromeda::andromeda(covariateRef = data.frame(covariateId = 11:15, analysisId = c(1, 1, 2, 2, 2)),
                                        analysisRef = data.frame(analysisId = c(1, 2), analysisName = c("A", "B")))
  on.exit(Andromeda::close(covariateData))
  d <- designSummary(X, c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1), rep(1:2, 5), 11:15, covariateData, minCellCount = 2)
  expect_equal(d$design$exactDuplicateColumns, 1)      # columns 1 and 2
  expect_equal(d$design$pairsCorrelationExact1, 1)
  expect_equal(d$design$nonZeros, 13)
  expect_equal(d$design$rowNonZeroTopMean, 3)         # rows 1-3 have 3 each
  expect_equal(d$columnTypes$columns[d$columnTypes$analysis == "B"], 3)
  expect_equal(sum(d$columnTypes$nonZeroShare), 1)
  d5 <- designSummary(X, c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1), rep(1:2, 5), 11:15, covariateData, minCellCount = 5)
  expect_equal(d5$design$columnRowsP10, -5)             # 1-4 persons per column: censored
  expect_true(all(d5$columnTypes$rowsMedian == -5))
})

test_that("fitTraceReport traces every fit of a model on the study's training set", {
  skip_on_cran()
  profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
  simulated <- PatientLevelPrediction::simulatePlpData(profile, n = 2000, seed = 11)
  prevalence <- simulated$covariateData$covariates %>% dplyr::count(.data$covariateId) %>% dplyr::collect()
  binary <- prevalence$covariateId[prevalence$n > 0.05 * 2000 & prevalence$n < 0.4 * 2000 & prevalence$covariateId != 1002]
  signal <- utils::head(binary, 9)
  coefficients <- c("(Intercept)" = -3, stats::setNames(rep(c(1.5, -1.2, 1.0), length.out = length(signal)), signal))
  plpData <- PatientLevelPrediction::simulatePlpBenchmarkData(simulated, coefficients, n = 8000, riskWindowStart = 1, riskWindowEnd = 365,
                                                              outcomeId = 3, targetOutcomeRate = 0.06, seed = 12, returnTruth = FALSE)
  settings <- PatientLevelPrediction::createStudyPopulationSettings(riskWindowStart = 1, riskWindowEnd = 365, requireTimeAtRisk = FALSE)
  source <- createPlpDataSource(plpData, "SIM", settings)
  task <- data.frame(taskId = "sim_task_defaults", outcomeId = 3, timeAtRiskDays = 365, featureSet = "defaults")
  folder <- tempfile("sar_trace_")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  methods <- fitAlgorithmComparisonMethods()[c("Ridge_ccd", "SAR_fixed_auto")]
  zipFile <- fitTraceReport(source, task, 100, methods, folder, folder)
  expect_true(file.exists(zipFile))
  design <- utils::read.csv(file.path(folder, "design.csv"))
  expect_equal(design$events, 100)
  summary <- utils::read.csv(file.path(folder, "trace_summary.csv"))
  expect_setequal(summary$method, names(methods))
  trace <- utils::read.csv(file.path(folder, "trace_Ridge_ccd.csv"))
  expect_setequal(unique(trace$stage), c("fold", "final"))
  expect_equal(sum(trace$stage == "final"), 1)
  expect_true(all(trace$iterations > 0 & trace$seconds >= 0 & trace$status == "SUCCESS"))
  expect_equal(summary$foldFits[summary$method == "Ridge_ccd"], sum(trace$stage == "fold"))
  expect_null(getOption("SignalAdaptiveStudy.fitTrace"))
})
