# The production data source: cohort definitions shipped with the package, task mapping, and the original task rules
# as PLP population settings (checked at the day boundaries). Against the original SQL on a real CDM:
# extras/production/validateCdmSource.R.

test_that("the cohort definitions cover every task", {
  cohorts <- studyCohortDefinitionSet()
  expect_identical(nrow(cohorts), 26L)
  expect_false(anyDuplicated(cohorts$cohortId) > 0)
  expect_true(all(nzchar(cohorts$json) & grepl("@target_cohort_id", cohorts$sql)))
  tasks <- originalCohortTasks()
  expect_identical(nrow(tasks), 26L)
  expect_true(all(tasks$targetId %in% cohorts$cohortId[cohorts$role == "target"]))
  expect_true(all(tasks$outcomeId %in% cohorts$cohortId[cohorts$role == "outcome"]))
  expect_identical(nrow(originalCohortTasks("IPCI")), 8L)
  catalogue <- catalogueCohortTasks()
  expect_setequal(catalogue$catalogueTask, c("SAR-017", "SAR-175", "SAR-287", "SAR-288"))
  expect_true(all(catalogue$targetId %in% cohorts$cohortId[cohorts$role == "target"]))
  expect_true(all(catalogue$outcomeId %in% cohorts$cohortId[cohorts$role == "outcome"]))
  all <- studyTasks()
  expect_identical(nrow(all), 30L)
  expect_false(anyDuplicated(all$taskId) > 0)
  expect_identical(taskPopulationSettings(all[all$taskId == "SAR-287_MDD_to_bipolar_1y", ])$priorOutcomeLookback, 999999)
  expect_identical(taskPopulationSettings(all[1, ])$priorOutcomeLookback, 366)
})

test_that("population settings follow the original rules at the day boundaries", {
  horizon <- 365
  # one person per case: (days of the outcome relative to index, days of observation after index)
  cases <- data.frame(
    outcomeDay = c(NA, -366, -365, 0, 1, horizon, horizon + 1, NA, 10),
    obsEnd = c(horizon, horizon, horizon, horizon, horizon, horizon, horizon, horizon - 1, horizon - 1),
    kept = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    outcome = c(0, 0, NA, NA, 1, 1, 0, NA, NA))
  n <- nrow(cases)
  cohorts <- data.frame(rowId = seq_len(n), subjectId = seq_len(n), targetId = 1,
                        cohortStartDate = as.Date("2015-01-01"), daysFromObsStart = 400, daysToCohortEnd = 0,
                        daysToObsEnd = cases$obsEnd, ageYear = 50, gender = 8507)
  attr(cohorts, "metaData") <- list(targetId = 1, attrition = data.frame(outcomeId = 2, description = "test",
                                                                          targetCount = n, uniquePeople = n, outcomes = 1))
  has <- !is.na(cases$outcomeDay)
  outcomes <- data.frame(rowId = which(has), outcomeId = 2, outcomeCount = 1, daysToEvent = cases$outcomeDay[has])
  attr(outcomes, "metaData") <- data.frame(outcomeIds = 2)
  covariateData <- Andromeda::andromeda(
    covariates = data.frame(rowId = 1, covariateId = 1, covariateValue = 1),
    covariateRef = data.frame(covariateId = 1, covariateName = "x", analysisId = 1, conceptId = 0),
    analysisRef = data.frame(analysisId = 1, analysisName = "x", domainId = "x", startDay = 0, endDay = 0,
                             isBinary = "Y", missingMeansZero = "Y"))
  class(covariateData) <- "CovariateData"
  plpData <- structure(list(cohorts = cohorts, outcomes = outcomes, covariateData = covariateData,
                            metaData = list(databaseDetails = list(), restrictPlpDataSettings = list())),
                       class = "plpData")
  population <- suppressMessages(PatientLevelPrediction::createStudyPopulation(
    plpData, 2, originalProductionPopulationSettings(horizon)))
  expect_setequal(population$rowId, which(cases$kept))
  expect_identical(as.numeric(population$outcomeCount[order(population$rowId)] > 0), cases$outcome[cases$kept])
})
