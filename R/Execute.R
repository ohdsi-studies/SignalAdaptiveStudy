# execute(): counts and models for task x database x budget x ordering x method, restartable.
#
# outputFolder/
#   data/<database>/<taskId>/            cached plpData (patient-level; never exported)
#   counts/<database>/<taskId>.csv       population and split counts
#   models/<database>/<taskId>/<run>/    result.csv (one row) and model/ (coefficients; private)
# A run whose result.csv exists is skipped; a failed run leaves error.txt and is tried again next time.
#
# A data source supplies, per task, the plpData, the population settings and the split of a budget and ordering
# (and optionally a faster counts function); execute() itself has no data- or split-specific choices.

newDataSource <- function(databaseName, splitLabel, loadPlpData, populationSettings, split, counts = NULL) {
  structure(list(databaseName = databaseName, splitLabel = splitLabel, loadPlpData = loadPlpData,
                 populationSettings = populationSettings, split = split, counts = counts), class = "sarDataSource")
}

#' Data source: an existing plpData object
#'
#' For simulated data and tests. The study split applies: PLP's default 75/25 split fixes the test set, each
#' ordering is a seeded random order of the training pool cut at the k-th event, and folds are stratified.
#'
#' @param plpData A `plpData` object.
#' @param databaseName Name used in the results.
#' @param populationSettings PLP population settings used for every task.
#' @return A data source for [execute()].
#' @export
createPlpDataSource <- function(plpData, databaseName, populationSettings) {
  newDataSource(databaseName, "study",
                loadPlpData = function(task, outputFolder) plpData,
                populationSettings = function(task) populationSettings,
                split = studySplit)
}

#' The four study methods
#'
#' SAR and SAL over the 12 shapes; ridge and LASSO as SAR and SAL with c = (0, 0, 0), through the same code.
#'
#' @return A named list of model settings.
#' @export
studyMethods <- function() {
  list(SAR = setSar(), SAL = setSal(),
       Ridge = setSar(shapes = c(0, 0, 0), name = "Ridge"), LASSO = setSal(shapes = c(0, 0, 0), name = "LASSO"))
}

# No frequency filter and no redundancy removal (all non-constant columns kept). normalize = TRUE divides by each
# covariate's maximum, which is 1 for binary covariates; FALSE would trigger
# https://github.com/OHDSI/PatientLevelPrediction/issues/688.
studyPreprocessSettings <- function() {
  PatientLevelPrediction::createPreprocessSettings(minFraction = 0, normalize = TRUE, removeRedundancy = FALSE)
}

aucMannWhitney <- function(y, p) {
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

testMetrics <- function(y, p, trainingY) {
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  b0 <- stats::plogis(smoothedIntercept(trainingY)) # the null model: the training (support) rate
  list(auc = aucMannWhitney(y, p), logLoss = -mean(y * log(p) + (1 - y) * log(1 - p)),
       nullLogLoss = -mean(y * log(b0) + (1 - y) * log(1 - b0)), brier = mean((p - y)^2))
}

packageVersions <- function() {
  cyclops <- utils::packageDescription("Cyclops")
  paste0("SignalAdaptiveStudy ", utils::packageVersion("SignalAdaptiveStudy"),
         "; PatientLevelPrediction ", utils::packageVersion("PatientLevelPrediction"),
         "; Cyclops ", cyclops$Version, if (!is.null(cyclops$RemoteSha)) paste0("@", substr(cyclops$RemoteSha, 1, 8)) else "",
         "; ", R.version.string)
}

# runPlp logs errors and carries on; the ERROR lines of its log in saveDirectory.
plpLogErrors <- function(saveDirectory) {
  lines <- unlist(lapply(list.files(saveDirectory, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE), readLines))
  errors <- grep("ERROR", lines, value = TRUE)
  if (length(errors) == 0) "no error in the runPlp log" else paste(errors, collapse = " | ")
}

# PLP's trainingTime text ("12.3 secs") in seconds.
trainingSeconds <- function(text) {
  if (is.null(text)) return(NULL)
  parts <- strsplit(text, " ")[[1]]
  as.numeric(parts[1]) * c(secs = 1, mins = 60, hours = 3600, days = 86400)[[parts[2]]]
}

runOne <- function(plpData, population, task, split, modelSettings, method, runFolder, dataSource, budget, ordering) {
  start <- Sys.time()
  scratch <- tempfile("sar_run_")
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  result <- PatientLevelPrediction::runPlp(
    plpData = plpData, outcomeId = task$outcomeId, analysisId = basename(runFolder), analysisName = basename(runFolder),
    populationSettings = dataSource$populationSettings(task),
    splitSettings = PatientLevelPrediction::createExistingSplitSettings(data.frame(rowId = split$rowId, index = split$index)),
    preprocessSettings = studyPreprocessSettings(),
    modelSettings = modelSettings,
    logSettings = PatientLevelPrediction::createLogSettings(verbosity = "WARN", logName = "runPlp"),
    executeSettings = PatientLevelPrediction::createExecuteSettings(
      runSplitData = TRUE, runSampleData = FALSE, runFeatureEngineering = FALSE, runPreprocessData = TRUE,
      runModelDevelopment = TRUE, runCovariateSummary = FALSE),
    saveDirectory = scratch
  )
  if (is.null(result$model) || isTRUE(result$model$noModel)) stop("model fit failed: ", plpLogErrors(scratch))
  test <- result$prediction[result$prediction$evaluationType == "Test", ]
  if (nrow(test) == 0) stop("no test predictions: ", plpLogErrors(scratch))
  trainingRows <- split$rowId[split$index > 0]
  trainingY <- as.numeric(population$outcomeCount[match(trainingRows, population$rowId)] > 0)
  m <- testMetrics(as.numeric(test$outcomeCount > 0), test$value, trainingY)
  fp <- result$model$trainDetails$finalModelParameters
  # other PLP models (e.g. Broken Adaptive Ridge) lack the SAR fields: NA, non-zeros from their coefficients
  value <- function(x) if (length(x) == 1) x else NA # BAR's variance, e.g., is one per coefficient
  shape <- unlist(fp$shape) %||% rep(NA, 3)
  if (is.null(fp$nonZero)) {
    co <- result$model$model$coefficients
    fp$nonZero <- if (is.data.frame(co)) sum(co$betas != 0 & co$covariateIds != "(Intercept)") else NA
  }
  row <- data.frame(
    database = dataSource$databaseName, taskId = task$taskId, budget = if (is.null(budget)) NA else budget,
    ordering = ordering, method = method, split = dataSource$splitLabel, search = value(fp$search),
    cSelection = value(modelSettings$param$cSelection),
    nTrain = attr(split, "trainingN"), eventsTrain = attr(split, "trainingEvents"),
    nTest = nrow(test), eventsTest = sum(test$outcomeCount > 0),
    auc = m$auc, logLoss = m$logLoss, nullLogLoss = m$nullLogLoss, brier = m$brier,
    interceptOnly = isTRUE(fp$interceptOnly), c1 = shape[1], c2 = shape[2], c3 = shape[3],
    lambda = value(fp$lambda), variance = value(fp$variance), cvLoss = value(fp$cvLoss),
    ordinateLoss = value(fp$ordinateLoss), interceptCvLoss = value(fp$interceptCvLoss), nonZero = fp$nonZero,
    kkt = value(fp$kkt), fitSeconds = value(fp$cvSeconds %||% trainingSeconds(result$model$trainDetails$trainingTime)),
    runSeconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    versions = packageVersions(), stringsAsFactors = FALSE
  )
  PatientLevelPrediction::savePlpModel(result$model, file.path(runFolder, "model"))
  utils::write.csv(row, file.path(runFolder, "result.csv"), row.names = FALSE)
  row
}

taskCounts <- function(dataSource, task, outputFolder) {
  if (!is.null(dataSource$counts)) return(dataSource$counts(task))
  plpData <- dataSource$loadPlpData(task, outputFolder)
  population <- PatientLevelPrediction::createStudyPopulation(plpData, task$outcomeId, dataSource$populationSettings(task))
  s <- dataSource$split(population, task, NULL, 1)
  list(n = nrow(population), poolN = attr(s, "poolN"), poolEvents = attr(s, "poolEvents"),
       testN = attr(s, "testN"), testEvents = attr(s, "testEvents"))
}

#' Run the study
#'
#' `mode = "counts"`: population, training-pool and test counts per task and the feasibility rule (>= 500 training
#' events and >= 80 test events); no models. `mode = "models"`: every task x budget x ordering x method, skipping
#' runs that have already finished.
#'
#' @param dataSource A data source, e.g. from [createPlpDataSource()].
#' @param tasks Data frame with one row per task: `taskId`, `outcomeId`, `timeAtRiskDays`.
#' @param outputFolder Folder for the data cache, counts and results.
#' @param mode `"counts"` or `"models"`.
#' @param budgets Event budgets k.
#' @param orderings Orderings (1, 2, 3).
#' @param methods Named list of model settings (default [studyMethods()]).
#' @return For counts, a data frame of counts; for models, the result rows of this and earlier runs.
#' @export
execute <- function(dataSource, tasks, outputFolder, mode = c("counts", "models"), budgets = c(100, 250, 500),
                    orderings = 1:3, methods = studyMethods()) {
  mode <- match.arg(mode)
  if (!inherits(dataSource, "sarDataSource")) stop("dataSource must come from a create...DataSource() function")
  if (is.null(names(methods)) || any(names(methods) == "")) stop("methods must be a named list")
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
  database <- dataSource$databaseName
  if (mode == "counts") {
    # restartable like the models: a task with counts is skipped; a failed task leaves error.txt and is tried again
    folder <- file.path(outputFolder, "counts", database)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    rows <- lapply(seq_len(nrow(tasks)), function(i) {
      task <- tasks[i, ]
      file <- file.path(folder, paste0(task$taskId, ".csv"))
      if (file.exists(file)) return(utils::read.csv(file))
      errorFile <- file.path(folder, paste0(task$taskId, "_error.txt"))
      counts <- tryCatch(taskCounts(dataSource, task, outputFolder), error = function(e) {
        writeLines(conditionMessage(e), errorFile)
        message(database, " ", task$taskId, " counts failed: ", conditionMessage(e))
        NULL
      })
      if (is.null(counts)) return(NULL)
      unlink(errorFile)
      row <- data.frame(database = database, taskId = task$taskId, n = counts$n, poolN = counts$poolN,
                        poolEvents = counts$poolEvents, testN = counts$testN, testEvents = counts$testEvents,
                        feasible = counts$poolEvents >= 500 && counts$testEvents >= 80)
      utils::write.csv(row, file, row.names = FALSE)
      row
    })
    return(do.call(rbind, rows))
  }
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    todo <- expand.grid(budget = budgets, ordering = orderings, method = names(methods), stringsAsFactors = FALSE)
    todo$folder <- file.path(outputFolder, "models", database, task$taskId,
                             sprintf("k%s_o%d_%s", todo$budget, todo$ordering, todo$method))
    todo <- todo[!file.exists(file.path(todo$folder, "result.csv")), , drop = FALSE]
    if (nrow(todo) == 0) next
    plpData <- dataSource$loadPlpData(task, outputFolder)
    population <- PatientLevelPrediction::createStudyPopulation(plpData, task$outcomeId, dataSource$populationSettings(task))
    keys <- unique(todo[, c("budget", "ordering")])
    for (j in seq_len(nrow(keys))) {
      budget <- keys$budget[j]; ordering <- keys$ordering[j]
      splitIds <- tryCatch(dataSource$split(population, task, budget, ordering), error = function(e) e)
      runs <- todo[todo$budget == budget & todo$ordering == ordering, ]
      for (r in seq_len(nrow(runs))) {
        dir.create(runs$folder[r], recursive = TRUE, showWarnings = FALSE)
        unlink(file.path(runs$folder[r], "error.txt"))
        outcome <- if (inherits(splitIds, "error")) splitIds else tryCatch(
          runOne(plpData, population, task, splitIds, methods[[runs$method[r]]], runs$method[r], runs$folder[r],
                 dataSource, budget, ordering),
          error = function(e) e)
        if (inherits(outcome, "error")) {
          writeLines(conditionMessage(outcome), file.path(runs$folder[r], "error.txt"))
          message(sprintf("%s %s k=%s o=%d %s failed: %s", database, task$taskId, budget, ordering,
                          runs$method[r], conditionMessage(outcome)))
        }
      }
    }
  }
  collectResults(outputFolder)
}

#' Collect the result rows
#'
#' @param outputFolder The `outputFolder` of [execute()].
#' @return All result rows found, one per finished run.
#' @export
collectResults <- function(outputFolder) {
  files <- list.files(file.path(outputFolder, "models"), pattern = "^result\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(data.frame())
  do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
}
