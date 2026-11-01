# repairTestAuc(): runs recorded without a test AUC get it from their saved model, identical to the AUC of the run, in
# one process and in workers; a run whose saved model does not reproduce the stored metrics is left alone.

test_that("missing test AUCs are repaired from the saved models, only when the stored metrics are reproduced", {
  skip_on_cran()
  profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
  source <- PatientLevelPrediction::simulatePlpData(profile, n = 2000, seed = 11)
  prevalence <- source$covariateData$covariates %>% dplyr::count(.data$covariateId) %>% dplyr::collect()
  signal <- utils::head(prevalence$covariateId[prevalence$n > 0.05 * 2000 & prevalence$n < 0.4 * 2000 &
                                                 prevalence$covariateId != 1002], 9)
  coefficients <- c("(Intercept)" = -3, stats::setNames(rep(c(1.5, -1.2, 1.0), length.out = length(signal)), signal))
  plpData <- PatientLevelPrediction::simulatePlpBenchmarkData(source, coefficients, n = 12000, riskWindowStart = 1,
                                                              riskWindowEnd = 365, outcomeId = 3,
                                                              targetOutcomeRate = 0.06, seed = 12, returnTruth = FALSE)
  dataSource <- createPlpDataSource(plpData, "SIM", PatientLevelPrediction::createStudyPopulationSettings(
    riskWindowStart = 1, riskWindowEnd = 365, requireTimeAtRisk = FALSE))
  tasks <- data.frame(taskId = c("task_a", "task_b"), outcomeId = 3, timeAtRiskDays = 365)
  outputFolder <- tempfile("sar_repair_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  results <- execute(dataSource, tasks, outputFolder, mode = "models", budgets = c(100, 250), orderings = 1:2,
                     methods = studyMethods()[c("SAR", "LASSO")])
  files <- list.files(file.path(outputFolder, "models"), "^result\\.csv$", recursive = TRUE, full.names = TRUE)
  stored <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  # as before the fix: AUC NA (here: every run of task_a and one of task_b)
  broken <- grepl("task_a", files) | grepl("task_b/k250_o2_SAR", files, fixed = TRUE)
  breakAuc <- function() for (i in which(broken)) { r <- stored[[i]]; r$auc <- NA; utils::write.csv(r, files[i], row.names = FALSE) }
  aucs <- function() vapply(files, function(f) utils::read.csv(f)$auc, numeric(1))

  breakAuc()
  repaired <- suppressMessages(repairTestAuc(dataSource, tasks, outputFolder))
  expect_identical(unname(aucs()), vapply(stored, `[[`, numeric(1), "auc"))
  expect_true(all(file.exists(file.path(dirname(files[broken]), "result_na_auc.csv"))))
  expect_false(any(file.exists(file.path(dirname(files[!broken]), "result_na_auc.csv"))))
  expect_false(anyNA(repaired$auc))

  # in workers; a run whose stored log loss the model does not give is not repaired
  breakAuc()
  tampered <- which(broken)[1]
  r <- utils::read.csv(files[tampered]); r$logLoss <- r$logLoss * 1.01; utils::write.csv(r, files[tampered], row.names = FALSE)
  suppressMessages(repairTestAuc(dataSource, tasks, outputFolder, workers = 2))
  expected <- vapply(stored, `[[`, numeric(1), "auc"); expected[tampered] <- NA
  expect_identical(unname(aucs()), expected)
  expect_match(readLines(file.path(outputFolder, "log_SIM.txt")), "AUC not repaired: the saved model does not reproduce",
               all = FALSE)
})
