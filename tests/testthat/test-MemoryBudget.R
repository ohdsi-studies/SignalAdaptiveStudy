test_that("unit memory estimates grow with the training set and the test set", {
  shape <- list(persons = 545505, entriesPerPerson = 311, eventRate = 1000 / 85467)
  estimates <- vapply(c(100, 250, 500, 1000), unitMemoryGB, numeric(1), shape = shape, duckdbGB = 1)
  expect_true(all(diff(estimates) > 0))
  expect_equal(estimates[4] - estimates[1], 311 * 900 / shape$eventRate * UNIT_MEMORY$trainBytes / 1e9, tolerance = 1e-6)
  # the training set cannot exceed the pool (75% of the persons)
  expect_equal(unitMemoryGB(1e9, shape, 0), unitMemoryGB(2e9, shape, 0))
  old <- options(sarUnitMemory = list(personBytes = 0, trainBytes = 0, testBytes = 0))
  expect_equal(unitMemoryGB(1000, shape, 1), UNIT_MEMORY$baseGB + 1)
  options(old)
  expect_gt(availableMemoryGB(), 0)
})

test_that("with a memory budget smaller than any unit, units run one at a time and all finish", {
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
  tasks <- data.frame(taskId = "task_a", outcomeId = 3, timeAtRiskDays = 365)
  tasks$budgets <- list(c(100, 250))
  folder <- tempfile("sar_budget_")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  old <- options(sarMemoryBudgetGB = 0.001)
  on.exit(options(old), add = TRUE)
  results <- suppressMessages(execute(dataSource, tasks, folder, mode = "models", orderings = 1:2,
                                      methods = studyMethods()[c("LASSO", "Ridge")], workers = 2))
  expect_equal(nrow(results), 8)
  log <- readLines(file.path(folder, "log_SIM.txt"))
  expect_true(any(grepl("memory rule live: estimates up to 0.0 GB", log)))
  expect_true(any(grepl("task_a: [0-9]+ non-zeros per person; estimated memory per unit", log)))
  expect_true(any(grepl("least memory available during the run: [0-9.]+ GB at", log)))
  # each unit's estimate against its worker's peak memory
  memory <- utils::read.csv(file.path(folder, "memory_SIM.csv"))
  expect_equal(nrow(memory), 4)
  expect_true(all(memory$peakAfterGB >= memory$peakBeforeGB & memory$peakAfterGB > 0.05))
  expect_true(all(memory$trainingRows > 0 & memory$entriesPerPerson > 0))
  stepMemory <- utils::read.csv(file.path(folder, "memory_steps_SIM.csv"))
  expect_true(all(c("prepare k100 o1", "k100_o1_LASSO") %in% stepMemory$step) && all(stepMemory$peakGB > 0))
  expect_true(any(grepl("(, peak [0-9.]+ GB|peak at most [0-9.]+ GB, the worker's peak from an earlier unit) \\(estimate [0-9.]+ GB\\) done", log)))
  # one unit at a time: each unit's steps (split, prepare, its runs) end before the next unit's begin
  steps <- utils::read.csv(file.path(folder, "steps_SIM.csv"))
  steps$unit <- ifelse(grepl("^(split|prepare) k", steps$step), sub("^(split|prepare) ", "", steps$step),
                       ifelse(grepl("^k[0-9]+_o[0-9]+_", steps$step), sub("^(k[0-9]+)_o([0-9]+)_.*", "\\1 o\\2", steps$step), NA))
  steps <- steps[!is.na(steps$unit), ]
  steps$start <- as.numeric(as.POSIXct(steps$started)); steps$end <- steps$start + steps$seconds
  spans <- do.call(rbind, lapply(split(steps, steps$unit), function(d) data.frame(start = min(d$start), end = max(d$end))))
  spans <- spans[order(spans$start), ]
  expect_equal(nrow(spans), 4)
  expect_true(all(spans$start[-1] >= spans$end[-nrow(spans)] - 1)) # started-times are rounded to the second
})

test_that("a unit starts when the memory free now holds it and what the running units may still add", {
  running <- data.frame(estimate = c(8, 4), started = c(0, 1000))
  # at t = 1100 the unit started at 1000 is recent (whole 4 GB to come), the other may still peak (0.6 x 8)
  expect_true(admitsUnit(5, running, available = 5 + 4 + 4.8 + 6, cap = 100, reserve = 6, now = 1100))
  expect_false(admitsUnit(5, running, available = 5 + 4 + 4.8 + 6 - 0.1, cap = 100, reserve = 6, now = 1100))
  # later both are past their start: only the largest may still peak
  expect_true(admitsUnit(5, running, available = 5 + 4.8 + 6, cap = 100, reserve = 6, now = 1300))
  # the cap on the estimates holds whatever is free
  expect_false(admitsUnit(5, running, available = 1000, cap = 16.9, reserve = 6, now = 1300))
  expect_true(admitsUnit(5, running, available = Inf, cap = 17, reserve = 6, now = 1300))
  expect_true(admitsUnit(5, running[0, ], available = 11, cap = 17, reserve = 6, now = 0))
  expect_gt(totalMemoryGB(), 0)
})

test_that("the largest waiting unit that fits starts first, the first when none fits", {
  expect_equal(unitToStart(c(2, 9, 5, 7), fits = c(1, 3, 4)), 4)
  expect_equal(unitToStart(c(2, 9, 5, 7), fits = integer()), 1)
  expect_equal(unitToStart(c(5, 5), fits = 1:2), 1)
})
