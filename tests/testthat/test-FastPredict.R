test_that("predictions from the test design equal predictPlp()'s to the last bit", {
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
  task <- data.frame(taskId = "sim", outcomeId = 3, timeAtRiskDays = 365)
  population <- PatientLevelPrediction::createStudyPopulation(plpData, 3, settings)
  data <- prepareRunData(plpData, population, studySplit(population, task, 100, 1))
  design <- testDesign(data$Test)
  # the simulated data hold a non-binary covariate (age), so the normalisation is exercised too
  expect_false(is.null(design$value))
  expect_true(any(data$Test$covariateData$covariates %>% dplyr::filter(.data$covariateValue != 1) %>% dplyr::count() %>% dplyr::pull() > 0))
  cache <- new.env()
  expect_identical(cachedTestDesign(cache, data$Test), design)
  for (settings in list(setSar(shapes = c(0, 1, -1), name = "SAR_fixed"), setSal(shapes = c(0, 0, 0), name = "LASSO"))) {
    model <- PatientLevelPrediction::fitPlp(data$Train, settings, analysisId = "a", analysisPath = tempfile())
    expect_equal(sum(model$prediction$evaluationType == "Train"), nrow(data$Train$labels))
    viaPlp <- PatientLevelPrediction::predictPlp(model, data$Test, data$Test$labels)
    fast <- predictFromDesign(model, design, data$Test$labels)
    expect_identical(fast$rowId, viaPlp$rowId)
    expect_identical(fast$value, viaPlp$value)
    # the training predictions from the matrix agree with the DuckDB route (summed in another order: equal, not identical)
    viaDuckdb <- predictSignalAdaptive(model, data$Train, data$Train$labels)
    train <- model$prediction[model$prediction$evaluationType == "Train", ]
    expect_equal(train$value[match(viaDuckdb$rowId, train$rowId)], viaDuckdb$value, tolerance = 1e-12)
  }
  expect_true(isStudyModel(setSar()))
  expect_false(isStudyModel(PatientLevelPrediction::setLassoLogisticRegression()))
})
