# Repair of test AUCs that execute() recorded as NA before aucMannWhitney() counted in doubles: as integers, n1 x n0
# overflowed beyond 2^31 (test sets of ~500,000 rows with thousands of events). The saved model of each such run
# predicts the test set again (the test set is the same for every budget and ordering of a task); the predictions must
# give the stored log loss, Brier score and null log loss, and then the AUC is written into result.csv (the original
# kept as result_na_auc.csv). No refitting.

# The runs of a task without a test AUC, one unit per budget and ordering (as execute()'s units, so the worker pool
# runs them).
repairUnits <- function(task, dataSource, outputFolder) {
  files <- list.files(file.path(outputFolder, "models", dataSource$databaseName, task$taskId), "^result\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  rows <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  missing <- vapply(rows, function(r) is.na(r$auc), logical(1))
  if (!any(missing)) return(list())
  runs <- do.call(rbind, lapply(which(missing), function(i) data.frame(
    folder = dirname(files[i]), budget = rows[[i]]$budget, ordering = rows[[i]]$ordering, method = rows[[i]]$method,
    stringsAsFactors = FALSE)))
  keys <- unique(runs[, c("budget", "ordering")])
  source <- list(databaseName = dataSource$databaseName, splitLabel = dataSource$splitLabel, split = dataSource$split)
  lapply(seq_len(nrow(keys)), function(j) {
    list(work = "repair", task = task, budget = keys$budget[j], ordering = keys$ordering[j],
         runs = runs[runs$budget == keys$budget[j] & runs$ordering == keys$ordering[j], ],
         populationSettings = dataSource$populationSettings(task), source = source)
  })
}

# One unit: the test set (built once per task and process, kept with the task's data), then per run the saved model's
# predictions, checked against the stored metrics. Returns the failures as messages.
repairUnit <- function(unit, loadPlpData) {
  task <- unit$task; budget <- unit$budget; ordering <- unit$ordering
  current <- taskData(unit, loadPlpData)
  population <- current$population
  split <- unit$source$split(population, task, budget, ordering)
  testRows <- sort(split$rowId[split$index < 0])
  if (!identical(current$testRows, testRows)) {
    splitSettings <- PatientLevelPrediction::createExistingSplitSettings(data.frame(rowId = split$rowId, index = split$index))
    current$test <- timed(task$taskId, "test data", PatientLevelPrediction::splitData(
      plpData = current$plpData, population = population, splitSettings = splitSettings)$Test)
    current$testRows <- testRows
  }
  trainingRows <- split$rowId[split$index > 0]
  trainingY <- as.numeric(population$outcomeCount[match(trainingRows, population$rowId)] > 0)
  failures <- character()
  for (r in seq_len(nrow(unit$runs))) {
    folder <- unit$runs$folder[r]
    outcome <- tryCatch(timed(task$taskId, paste("repair", basename(folder)), {
      model <- PatientLevelPrediction::loadPlpModel(file.path(folder, "model"))
      test <- PatientLevelPrediction::predictPlp(plpModel = model, plpData = current$test,
                                                 population = current$test$labels)
      m <- testMetrics(as.numeric(test$outcomeCount > 0), test$value, trainingY)
      row <- utils::read.csv(file.path(folder, "result.csv"), stringsAsFactors = FALSE)
      same <- function(a, b) isTRUE(abs(a - b) <= 1e-9 * abs(b))
      if (nrow(test) != row$nTest || !same(m$logLoss, row$logLoss) || !same(m$brier, row$brier) ||
          !same(m$nullLogLoss, row$nullLogLoss)) {
        stop(sprintf("the saved model does not reproduce the stored metrics (log loss %.12g, stored %.12g)",
                     m$logLoss, row$logLoss))
      }
      if (!file.exists(file.path(folder, "result_na_auc.csv"))) {
        file.copy(file.path(folder, "result.csv"), file.path(folder, "result_na_auc.csv"))
      }
      row$auc <- m$auc
      utils::write.csv(row, file.path(folder, "result.csv"), row.names = FALSE)
    }), error = function(e) e)
    if (inherits(outcome, "error")) {
      failures <- c(failures, sprintf("%s %s %s: AUC not repaired: %s", unit$source$databaseName, task$taskId,
                                      basename(folder), conditionMessage(outcome)))
    }
  }
  failures
}

# Repairs every run of these tasks whose result.csv has no test AUC, in `workers` processes as execute() runs models;
# logged to log_<database>.txt and steps_<database>.csv. Returns the result rows.
repairTestAuc <- function(dataSource, tasks, outputFolder, workers = 1) {
  database <- dataSource$databaseName
  logger <- ParallelLogger::createLogger(name = "SignalAdaptiveStudy", threshold = "TRACE", appenders = list(
    ParallelLogger::createFileAppender(layout = ParallelLogger::layoutTimestamp,
                                       fileName = file.path(outputFolder, sprintf("log_%s.txt", database)))))
  ParallelLogger::unregisterLogger("SignalAdaptiveStudy", silent = TRUE)
  ParallelLogger::registerLogger(logger)
  .timing$file <- file.path(outputFolder, sprintf("steps_%s.csv", database))
  on.exit({
    ParallelLogger::unregisterLogger("SignalAdaptiveStudy", silent = TRUE)
    .timing$file <- NULL
    clearTaskData()
  }, add = TRUE)
  if (is.null(getOption("duckdb.home")) && Sys.getenv("DUCKDB_R_HOME") == "") { # see execute()
    Sys.setenv(DUCKDB_R_HOME = file.path(tempdir(), "duckdb"))
    on.exit(Sys.unsetenv("DUCKDB_R_HOME"), add = TRUE)
  }
  clearTaskData()
  unitsOf <- function(task) repairUnits(task, dataSource, outputFolder)
  report(sprintf("repairTestAuc(%s): %d run(s) without a test AUC", database,
                 sum(vapply(seq_len(nrow(tasks)), function(i) sum(vapply(unitsOf(tasks[i, ]), function(u) nrow(u$runs),
                                                                         numeric(1))), numeric(1)))))
  if (workers > 1) {
    runInWorkers(tasks, unitsOf, dataSource, outputFolder, workers)
  } else {
    for (i in seq_len(nrow(tasks))) {
      task <- tasks[i, ]
      for (unit in unitsOf(task)) {
        for (failure in repairUnit(unit, function() dataSource$loadPlpData(task, outputFolder))) report(failure, warn = TRUE)
      }
    }
  }
  results <- collectResults(outputFolder)
  report(sprintf("repairTestAuc(%s): %d run(s) still without a test AUC", database, sum(is.na(results$auc))))
  results
}
