# Smoke test for a production machine: the whole pipeline on simulated data (no database, no download), about a
# minute. Run after renv::restore(); send back the printed table and the zip (simulated data only).
#
#   Rscript extras/production/smokeTest.R [outputFolder]

library(SignalAdaptiveStudy)
args <- commandArgs(trailingOnly = TRUE)
outputFolder <- if (length(args)) args[1] else file.path(tempdir(), "sar_smoke_test")

profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
source <- PatientLevelPrediction::simulatePlpData(profile, n = 2000, seed = 11)
prevalence <- source$covariateData$covariates |> dplyr::count(covariateId) |> dplyr::collect()
signal <- utils::head(prevalence$covariateId[prevalence$n > 100 & prevalence$n < 800 & prevalence$covariateId != 1002], 9)
coefficients <- c("(Intercept)" = -3, stats::setNames(rep(c(1.5, -1.2, 1.0), length.out = length(signal)), signal))
plpData <- PatientLevelPrediction::simulatePlpBenchmarkData(source, coefficients, n = 12000, riskWindowStart = 1,
                                                            riskWindowEnd = 365, outcomeId = 3, targetOutcomeRate = 0.06,
                                                            seed = 12, returnTruth = FALSE)
dataSource <- createPlpDataSource(plpData, "SIMULATED", PatientLevelPrediction::createStudyPopulationSettings(
  riskWindowStart = 1, riskWindowEnd = 365, requireTimeAtRisk = FALSE))
tasks <- data.frame(taskId = "smoke_test", outcomeId = 3, timeAtRiskDays = 365)

start <- Sys.time()
counts <- execute(dataSource, tasks, outputFolder, mode = "counts")
results <- execute(dataSource, tasks, outputFolder, mode = "models", budgets = 100, orderings = 1,
                   methods = studyMethods())
zipFile <- exportResults(outputFolder)
print(results[, c("method", "nTrain", "eventsTrain", "auc", "nonZero", "fitSeconds")], row.names = FALSE)
message(sprintf("smoke test %s in %.1f min; zip: %s", if (all(results$auc > 0.6)) "PASSED" else "FAILED (AUC <= 0.6)",
                as.numeric(difftime(Sys.time(), start, units = "mins")), zipFile))
message(unique(results$versions))
