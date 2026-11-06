# execute() with worker processes (mirai): the same results as in one process, the step timings and worker logs, and a
# worker that dies is replaced and its unit run again.

test_that("workers give the same results as one process, and a worker that dies is replaced", {
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
  tasks$budgets <- list(c(100, 250), 100) # each task its own budgets
  methods <- studyMethods()[c("SAR", "LASSO")]
  folders <- c(serial = tempfile("sar_serial_"), workers = tempfile("sar_workers_"), crash = tempfile("sar_crash_"))
  on.exit(unlink(folders, recursive = TRUE), add = TRUE)

  serial <- execute(dataSource, tasks, folders[["serial"]], mode = "models", orderings = 1:2, methods = methods)
  parallel <- suppressMessages(execute(dataSource, tasks, folders[["workers"]], mode = "models", orderings = 1:2,
                                       methods = methods, workers = 2))
  expect_identical(nrow(serial), 12L)
  key <- function(r) r[order(r$taskId, r$budget, r$ordering, r$method), setdiff(names(r), grep("Seconds$", names(r), value = TRUE))]
  expect_identical(`rownames<-`(key(parallel), NULL), `rownames<-`(key(serial), NULL))
  # the workers' step timings reach the steps file; each worker has its log
  steps <- utils::read.csv(file.path(folders[["workers"]], "steps_SIM.csv"))
  expect_setequal(steps$step[grepl("^k", steps$step)], unique(sprintf("k%s_o%d_%s", serial$budget, serial$ordering, serial$method)))
  expect_true(all(steps$ok) && all(c("load plpData", "population", "split k250 o1", "prepare k250 o1") %in% steps$step))
  expect_gte(length(list.files(file.path(folders[["workers"]], "workers"), "^log_SIM_[0-9]+\\.txt$")), 1)

  # a worker that dies at its first fit (as on a stack overflow or out of memory) is replaced; the unit runs again
  flag <- normalizePath(tempfile("sar_crashed_"), winslash = "/", mustWork = FALSE)
  crashing <- methods$LASSO
  crashing$fitFunction <- sprintf("function(...) { if (!file.exists('%s')) { file.create('%s'); tools::pskill(Sys.getpid(), tools::SIGKILL) }; SignalAdaptiveStudy::fitSignalAdaptive(...) }", flag, flag)
  crash <- suppressMessages(execute(dataSource, tasks[1, ], folders[["crash"]], mode = "models", orderings = 1,
                                    methods = list(LASSO = crashing), workers = 2))
  expect_true(file.exists(flag))
  expect_identical(nrow(crash), 2L)
  expect_identical(crash$auc[order(crash$budget)],
                   serial$auc[serial$taskId == "task_a" & serial$ordering == 1 & serial$method == "LASSO"][
                     order(serial$budget[serial$taskId == "task_a" & serial$ordering == 1 & serial$method == "LASSO"])])
  expect_match(readLines(file.path(folders[["crash"]], "log_SIM.txt")), "the worker stopped", all = FALSE)
})
