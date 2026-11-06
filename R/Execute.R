# execute(): counts and models for task x database x budget x ordering x method, restartable.
#
# outputFolder/
#   data/<database>/<taskId>/            cached plpData (patient-level; never exported)
#   counts/<database>/<taskId>.csv       population and split counts
#   models/<database>/<taskId>/<run>/    result.csv (one row), plpLog.txt, model/ (coefficients; private) and
#                                        fingerprint.json (the settings of the run, R/Fingerprint.R)
#   log_<database>.txt                   timestamped log of every execute() call (PLP, FeatureExtraction, step times)
#   steps_<database>.csv                 duration of each step: cohorts, extraction, cache, population, split, the
#                                        data of a budget and ordering (prepare), run
#   workers/log_<database>_<pid>.txt     with workers > 1: the log of each worker process
#   tmp/                                 with workers > 1: the workers' temporary files (emptied at the next start)
# A run whose result.csv exists is skipped; a failed run leaves error.txt and is tried again next time. A finished run
# whose recorded settings differ from the current ones stops execute() (R/Fingerprint.R).
#
# A data source supplies, per task, the plpData, the population settings and the split of a budget and ordering
# (and optionally a faster counts function and the folder of the saved plpData); execute() itself has no data- or
# split-specific choices.

# Step timings of the current execute() call: appended to .timing$file (steps_<database>.csv) when execute() set it,
# or collected in .timing$rows in a worker, which returns them to execute(), the only writer of the file; also written
# to the log. A step that fails is recorded too (ok = FALSE).
.timing <- new.env()

timed <- function(taskId, step, expr) {
  start <- Sys.time()
  ok <- FALSE
  on.exit({
    seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    ParallelLogger::logInfo(sprintf("%s: %s %s %.1f secs", taskId, step, if (ok) "took" else "failed after", seconds))
    row <- data.frame(taskId = taskId, step = step, started = format(start, "%Y-%m-%d %H:%M:%S"),
                      seconds = round(seconds, 2), ok = ok)
    if (!is.null(.timing$rows)) {
      .timing$rows[[length(.timing$rows) + 1]] <- row
    } else if (!is.null(.timing$file)) appendSteps(row, .timing$file)
  })
  value <- expr
  ok <- TRUE
  value
}

appendSteps <- function(rows, file) {
  if (is.null(rows) || nrow(rows) == 0) return(invisible())
  if (file.exists(file)) {
    utils::write.table(rows, file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  } else utils::write.csv(rows, file, row.names = FALSE)
}

# plpDataFolder(task, outputFolder), optional: the folder of the task's saved plpData (extracted first if need be), which
# workers load; without it, execute() with workers > 1 saves the plpData under outputFolder/tmp.
# dataFingerprint(task), optional: a named list of what the task's data depend on (cohorts, covariates, sample), for
# the fingerprints of counts and runs (R/Fingerprint.R).
newDataSource <- function(databaseName, splitLabel, loadPlpData, populationSettings, split, counts = NULL,
                          plpDataFolder = NULL, dataFingerprint = NULL) {
  structure(list(databaseName = databaseName, splitLabel = splitLabel, loadPlpData = loadPlpData,
                 populationSettings = populationSettings, split = split, counts = counts,
                 plpDataFolder = plpDataFolder, dataFingerprint = dataFingerprint), class = "sarDataSource")
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

# Counts as doubles: as integers, n1 * n0 overflows (NA) beyond 2^31, e.g. 8,000 events in a test set of 500,000.
aucMannWhitney <- function(y, p) {
  r <- rank(p)
  n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# DeLong's standard error of an AUC (Sun & Xu's fast form), from midranks.
aucDeLongSe <- function(y, p) {
  n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  if (n1 < 2 || n0 < 2) return(NA_real_)
  all <- rank(p)
  v10 <- (all[y == 1] - rank(p[y == 1])) / n0
  v01 <- 1 - (all[y == 0] - rank(p[y == 0])) / n1
  sqrt(stats::var(v10) / n1 + stats::var(v01) / n0)
}

# Calibration on the test rows: logistic slope and intercept on the model's linear predictor, the intercept with the
# slope fixed at 1 (calibration in the large), and expected over observed. NA where the predictions are constant.
calibrationMetrics <- function(y, p) {
  eta <- stats::qlogis(p)
  out <- list(calibrationSlope = NA_real_, calibrationIntercept = NA_real_, calibrationInTheLarge = NA_real_,
              expectedOverObserved = if (sum(y) > 0) sum(p) / sum(y) else NA_real_)
  if (length(unique(y)) < 2) return(out)
  fit <- function(...) suppressWarnings(tryCatch(stats::glm.fit(..., family = stats::binomial()), error = function(e) NULL))
  large <- fit(rep(1, length(y)), y, offset = eta)
  if (!is.null(large) && large$converged) out$calibrationInTheLarge <- unname(large$coefficients[1])
  if (stats::sd(eta) > 0) {
    free <- fit(cbind(1, eta), y)
    if (!is.null(free) && free$converged) {
      out$calibrationIntercept <- unname(free$coefficients[1])
      out$calibrationSlope <- unname(free$coefficients[2])
    }
  }
  out
}

testMetrics <- function(y, p, trainingY) {
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  b0 <- stats::plogis(smoothedIntercept(trainingY)) # the null model: the training (support) rate
  c(list(auc = aucMannWhitney(y, p), aucSe = aucDeLongSe(y, p), logLoss = -mean(y * log(p) + (1 - y) * log(1 - p)),
         nullLogLoss = -mean(y * log(b0) + (1 - y) * log(1 - b0)), brier = mean((p - y)^2)),
    calibrationMetrics(y, p))
}

# The same metrics on the test rows with the full time at risk: the secondary test set without patients lost to
# follow-up (cases and non-cases), for the rule that keeps them in the primary one.
followUpMetrics <- function(test, population, settings, trainingY) {
  horizon <- settings$riskWindowEnd - settings$riskWindowStart + 1
  atRisk <- population$timeAtRisk[match(test$rowId, population$rowId)]
  keep <- !is.na(atRisk) & atRisk >= horizon
  y <- as.numeric(test$outcomeCount[keep] > 0)
  out <- list(nTestFullFollowUp = sum(keep), eventsTestFullFollowUp = sum(y), aucFullFollowUp = NA_real_,
              logLossFullFollowUp = NA_real_, calibrationInTheLargeFullFollowUp = NA_real_)
  if (length(unique(y)) == 2) {
    m <- testMetrics(y, test$value[keep], trainingY)
    out$aucFullFollowUp <- m$auc; out$logLossFullFollowUp <- m$logLoss
    out$calibrationInTheLargeFullFollowUp <- m$calibrationInTheLarge
  }
  out
}

packageVersions <- function() {
  cyclops <- utils::packageDescription("Cyclops")
  paste0("SignalAdaptiveStudy ", utils::packageVersion("SignalAdaptiveStudy"),
         "; PatientLevelPrediction ", utils::packageVersion("PatientLevelPrediction"),
         "; Cyclops ", cyclops$Version, if (!is.null(cyclops$RemoteSha)) paste0("@", substr(cyclops$RemoteSha, 1, 8)) else "",
         "; ", R.version.string)
}

# PLP's trainingTime text ("12.3 secs") in seconds.
trainingSeconds <- function(text) {
  if (is.null(text)) return(NULL)
  parts <- strsplit(text, " ")[[1]]
  as.numeric(parts[1]) * c(secs = 1, mins = 60, hours = 3600, days = 86400)[[parts[2]]]
}

# The data of one budget and ordering, shared by its runs (one per method): PLP's split (the study split as an existing
# split) and PLP's preprocessing of the training set, as runPlp() does them before the fit.
prepareRunData <- function(plpData, population, split) {
  splitSettings <- PatientLevelPrediction::createExistingSplitSettings(data.frame(rowId = split$rowId, index = split$index))
  data <- PatientLevelPrediction::splitData(plpData = plpData, population = population, splitSettings = splitSettings)
  data$Train$covariateData <- PatientLevelPrediction::preprocessData(data$Train$covariateData, studyPreprocessSettings())
  data
}

# One run: what runPlp() does after the preprocessing, the fit and the test predictions (fitPlp(), predictPlp()), but
# not its evaluation and saved result, which the study does not use.
runOne <- function(data, population, task, split, modelSettings, method, runFolder, dataSource, budget, ordering,
                   fingerprint = NULL) {
  start <- Sys.time()
  scratch <- tempfile("sar_run_")
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  # PLP's timestamped log of the run, in the run folder
  logger <- ParallelLogger::createLogger(name = "SignalAdaptiveStudyRun", threshold = "INFO", appenders = list(
    ParallelLogger::createFileAppender(layout = ParallelLogger::layoutTimestamp,
                                       fileName = file.path(runFolder, "plpLog.txt"))))
  ParallelLogger::registerLogger(logger)
  on.exit(ParallelLogger::unregisterLogger("SignalAdaptiveStudyRun", silent = TRUE), add = TRUE)
  # Warnings of the fit (e.g. Cyclops's failed cross-validation folds) go to the result row and the run's log.
  fitWarnings <- character()
  model <- withCallingHandlers(
    PatientLevelPrediction::fitPlp(trainData = data$Train, modelSettings = modelSettings,
                                   analysisId = basename(runFolder), analysisPath = scratch),
    warning = function(w) {
      fitWarnings <<- c(fitWarnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  for (w in unique(fitWarnings)) ParallelLogger::logWarn(paste("fit warning:", w))
  model$prediction <- NULL # the CV predictions, which runPlp() also takes out of the model
  test <- PatientLevelPrediction::predictPlp(plpModel = model, plpData = data$Test, population = data$Test$labels)
  if (nrow(test) == 0) stop("no test predictions")
  trainingRows <- split$rowId[split$index > 0]
  trainingY <- as.numeric(population$outcomeCount[match(trainingRows, population$rowId)] > 0)
  m <- testMetrics(as.numeric(test$outcomeCount > 0), test$value, trainingY)
  full <- followUpMetrics(test, population, dataSource$populationSettings(task), trainingY)
  fp <- model$trainDetails$finalModelParameters
  # other PLP models (e.g. Broken Adaptive Ridge) lack the SAR fields: NA, non-zeros from their coefficients
  value <- function(x) if (length(x) == 1) x else NA # BAR's variance, e.g., is one per coefficient
  shape <- unlist(fp$shape) %||% rep(NA, 3)
  if (is.null(fp$nonZero)) {
    co <- model$model$coefficients
    fp$nonZero <- if (is.data.frame(co)) sum(co$betas != 0 & co$covariateIds != "(Intercept)") else NA
  }
  # Convergence: the study's search reports its fold and final fits; PLP's Cyclops models their status, where
  # "REACHED MAXIMUM NUMBER OF ITERATIONS" or "ILL CONDITIONED" means PLP returned an all-zero model.
  status <- fp$finalStatus %||% model$model$modelStatus %||% NA_character_
  if (!is.na(status) && !status %in% c("SUCCESS", "OK", "INTERCEPT_ONLY")) {
    ParallelLogger::logWarn(sprintf("%s %s: model status %s", task$taskId, basename(runFolder), status))
  }
  row <- data.frame(
    database = dataSource$databaseName, taskId = task$taskId, featureSet = task$featureSet %||% NA,
    budget = if (is.null(budget)) NA else budget,
    ordering = ordering, method = method, split = dataSource$splitLabel, search = value(fp$search),
    cSelection = value(modelSettings$param$cSelection),
    nTrain = attr(split, "trainingN"), eventsTrain = attr(split, "trainingEvents"),
    nTest = nrow(test), eventsTest = sum(test$outcomeCount > 0),
    auc = m$auc, aucSe = m$aucSe, logLoss = m$logLoss, nullLogLoss = m$nullLogLoss, brier = m$brier,
    calibrationSlope = m$calibrationSlope, calibrationIntercept = m$calibrationIntercept,
    calibrationInTheLarge = m$calibrationInTheLarge, expectedOverObserved = m$expectedOverObserved,
    nTestFullFollowUp = full$nTestFullFollowUp, eventsTestFullFollowUp = full$eventsTestFullFollowUp,
    aucFullFollowUp = full$aucFullFollowUp, logLossFullFollowUp = full$logLossFullFollowUp,
    calibrationInTheLargeFullFollowUp = full$calibrationInTheLargeFullFollowUp,
    interceptOnly = isTRUE(fp$interceptOnly), c1 = shape[1], c2 = shape[2], c3 = shape[3],
    lambda = value(fp$lambda), variance = value(fp$variance), cvLoss = value(fp$cvLoss),
    ordinateLoss = value(fp$ordinateLoss), interceptCvLoss = value(fp$interceptCvLoss), nonZero = fp$nonZero,
    kkt = value(fp$kkt), modelStatus = status, finalIterations = value(fp$finalIterations),
    foldFitsFailed = value(fp$foldFitsFailed), foldFitStatus = value(fp$foldFitStatus),
    fitWarnings = paste(unique(fitWarnings), collapse = " | "), fitSeconds = value(fp$cvSeconds %||% trainingSeconds(model$trainDetails$trainingTime)),
    runSeconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    versions = packageVersions(), stringsAsFactors = FALSE
  )
  PatientLevelPrediction::savePlpModel(model, file.path(runFolder, "model"))
  if (!is.null(fingerprint)) writeFingerprint(fingerprint, file.path(runFolder, FINGERPRINT_FILE))
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

# A task's runs (budget x ordering x method) and their folders.
runFolders <- function(task, budgets, orderings, methods, dataSource, outputFolder) {
  runs <- expand.grid(budget = budgets, ordering = orderings, method = names(methods), stringsAsFactors = FALSE)
  runs$folder <- file.path(outputFolder, "models", dataSource$databaseName, task$taskId,
                           sprintf("k%s_o%d_%s", runs$budget, runs$ordering, runs$method))
  runs
}

# Finished runs whose recorded settings differ from the current ones: one line each (folder: settings that changed).
# Finished runs without a fingerprint (from before fingerprints) are counted in attr "legacy".
changedRuns <- function(task, budgets, orderings, methods, dataSource, outputFolder) {
  runs <- runFolders(task, budgets, orderings, methods, dataSource, outputFolder)
  runs <- runs[file.exists(file.path(runs$folder, "result.csv")), , drop = FALSE]
  legacy <- 0
  changed <- character()
  for (r in seq_len(nrow(runs))) {
    fp <- runFingerprint(dataSource, task, methods[[runs$method[r]]], runs$budget[r], runs$ordering[r])
    changes <- fingerprintChanges(fp, file.path(runs$folder[r], FINGERPRINT_FILE))
    if (is.null(changes)) legacy <- legacy + 1
    else if (length(changes)) changed <- c(changed, sprintf("%s: %s", runs$folder[r], paste(changes, collapse = ", ")))
  }
  structure(changed, legacy = legacy)
}

# The runs still to do of a task, grouped by budget and ordering: one unit per budget and ordering, with what its runs
# need (and nothing of the data source beyond that, so a unit can go to a worker process).
modelUnits <- function(task, budgets, orderings, methods, dataSource, outputFolder) {
  todo <- runFolders(task, budgets, orderings, methods, dataSource, outputFolder)
  todo <- todo[!file.exists(file.path(todo$folder, "result.csv")), , drop = FALSE]
  keys <- unique(todo[, c("budget", "ordering")])
  source <- list(databaseName = dataSource$databaseName, splitLabel = dataSource$splitLabel, split = dataSource$split,
                 populationSettings = dataSource$populationSettings)
  lapply(seq_len(nrow(keys)), function(j) {
    runs <- todo[todo$budget == keys$budget[j] & todo$ordering == keys$ordering[j], c("method", "folder")]
    fingerprints <- lapply(seq_len(nrow(runs)), function(r) {
      runFingerprint(dataSource, task, methods[[runs$method[r]]], keys$budget[j], keys$ordering[j])
    })
    list(task = task, budget = keys$budget[j], ordering = keys$ordering[j], runs = runs,
         methods = methods[unique(runs$method)], populationSettings = dataSource$populationSettings(task),
         source = source, fingerprints = stats::setNames(fingerprints, runs$folder))
  })
}

# The plpData and population of the task in hand, kept between its budgets and orderings (in execute()'s process, or
# in each worker); emptied at the start and end of execute().
.taskData <- new.env()

clearTaskData <- function() {
  rm(list = ls(.taskData, all.names = TRUE), envir = .taskData)
  invisible(gc(verbose = FALSE)) # closes the Andromeda of a plpData loaded here
}

taskData <- function(unit, loadPlpData) {
  key <- paste(unit$source$databaseName, unit$task$taskId)
  if (!identical(.taskData$key, key)) {
    clearTaskData()
    plpData <- loadPlpData()
    .taskData$population <- timed(unit$task$taskId, "population", PatientLevelPrediction::createStudyPopulation(
      plpData, unit$task$outcomeId, unit$populationSettings))
    .taskData$plpData <- plpData
    .taskData$key <- key
  }
  .taskData
}

# One budget and ordering of a task: its split and preprocessed data, shared by its runs (one per method still to do).
# loadPlpData() gives the task's plpData when it is not already in hand. Returns the failures as messages (a failed
# run leaves error.txt).
runUnit <- function(unit, loadPlpData) {
  task <- unit$task; budget <- unit$budget; ordering <- unit$ordering
  current <- taskData(unit, loadPlpData)
  population <- current$population
  splitIds <- tryCatch(timed(task$taskId, sprintf("split k%s o%d", budget, ordering),
                             unit$source$split(population, task, budget, ordering)), error = function(e) e)
  data <- if (inherits(splitIds, "error")) splitIds else tryCatch(
    timed(task$taskId, sprintf("prepare k%s o%d", budget, ordering), prepareRunData(current$plpData, population, splitIds)),
    error = function(e) e)
  failures <- character()
  for (r in seq_len(nrow(unit$runs))) {
    folder <- unit$runs$folder[r]; method <- unit$runs$method[r]
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    unlink(file.path(folder, "error.txt"))
    outcome <- if (inherits(data, "error")) data else tryCatch(
      timed(task$taskId, basename(folder), runOne(data, population, task, splitIds, unit$methods[[method]], method,
                                                  folder, unit$source, budget, ordering, unit$fingerprints[[folder]])),
      error = function(e) e)
    if (inherits(outcome, "error")) {
      writeLines(conditionMessage(outcome), file.path(folder, "error.txt"))
      failures <- c(failures, sprintf("%s %s k=%s o=%d %s failed: %s", unit$source$databaseName, task$taskId, budget,
                                      ordering, method, conditionMessage(outcome)))
    }
  }
  failures
}

# Marks the unfinished runs of a unit as failed (error.txt), for a unit whose worker failed outside its runs.
failUnit <- function(unit, reason) {
  for (folder in unit$runs$folder[!file.exists(file.path(unit$runs$folder, "result.csv"))]) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    writeLines(reason, file.path(folder, "error.txt"))
  }
  sprintf("%s %s k=%s o=%d failed: %s", unit$source$databaseName, unit$task$taskId, unit$budget, unit$ordering, reason)
}

# Message on the console and in the log.
report <- function(text, warn = FALSE) {
  message(text)
  if (warn) ParallelLogger::logWarn(text) else ParallelLogger::logInfo(text)
}

# The source folder of the package when it was loaded with pkgload::load_all() (tests), for the workers to load it
# the same way; NULL when installed.
devPackagePath <- function() {
  ns <- asNamespace("SignalAdaptiveStudy")
  if (exists(".__DEVTOOLS__", envir = ns, inherits = FALSE)) getNamespaceInfo(ns, "path") else NULL
}

# Memory limit (GB) of each Andromeda (duckdb) in a worker: Andromeda's default is 20% of the RAM for every Andromeda,
# and a worker has up to three open (the plpData, its split, the preprocessed training data). Here the Andromedas of all
# workers together are held to half the RAM, beyond which duckdb writes to disk; the option andromedaMemoryLimit
# overrides it.
workerMemoryLimit <- function(workers) {
  if (!is.null(getOption("andromedaMemoryLimit"))) return(getOption("andromedaMemoryLimit"))
  ram <- tryCatch(as.numeric(memuse::Sys.meminfo()$totalram) / 1024^3, error = function(e) 16)
  round(0.5 * ram / workers / 3, 2)
}

# A unit in a worker process (runInWorkers()): one thread and a memory limit for each Andromeda, the worker's own log,
# and the step timings returned to execute(), which writes them. The random number generator is execute()'s: mirai
# daemons use L'Ecuyer-CMRG, under which the study's seeds would give other orderings and folds.
workerUnit <- function(payload, setup) {
  unit <- unserialize(payload)
  do.call(RNGkind, as.list(setup$rngKind))
  options(andromedaThreads = 1, andromedaMemoryLimit = setup$memoryLimit)
  dir.create(file.path(setup$outputFolder, "workers"), showWarnings = FALSE)
  ParallelLogger::unregisterLogger("SignalAdaptiveStudyWorker", silent = TRUE)
  ParallelLogger::registerLogger(ParallelLogger::createLogger(name = "SignalAdaptiveStudyWorker", threshold = "TRACE",
    appenders = list(ParallelLogger::createFileAppender(layout = ParallelLogger::layoutTimestamp, fileName = file.path(
      setup$outputFolder, "workers", sprintf("log_%s_%d.txt", setup$database, Sys.getpid()))))))
  .timing$rows <- list()
  on.exit(.timing$rows <- NULL)
  work <- if (identical(unit$work, "repair")) repairUnit else runUnit # repairUnit(): R/Repair.R
  failures <- work(unit, function() timed(unit$task$taskId, "load plpData",
                                          PatientLevelPrediction::loadPlpData(unit$plpDataFolder)))
  list(steps = do.call(rbind, .timing$rows), failures = failures)
}

# execute() with workers > 1: each unit (a budget and ordering of a task) goes to the next free worker process (mirai
# daemons), the largest first. The data stay with this process: each task's plpData is extracted (or found in the
# cache) here, one task at a time, and its units are queued as soon as it is ready, so the workers start on the first
# task while the next is extracted. A worker runs one thread for duckdb, BLAS and OpenMP, has a memory limit for each
# Andromeda and writes its temporary files under outputFolder/tmp. A worker that dies (e.g. out of memory) is replaced
# and its unit tried once more.
runInWorkers <- function(tasks, unitsOf, dataSource, outputFolder, workers) {
  profile <- "SignalAdaptiveStudy"
  database <- dataSource$databaseName
  tmp <- file.path(outputFolder, "tmp")
  unlink(tmp, recursive = TRUE) # what workers that died left
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  # environment variables the workers inherit when they start (TMPDIR: R's temporary folder, and so Andromeda's)
  env <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", TMPDIR = normalizePath(tmp))
  old <- Sys.getenv(names(env), unset = NA)
  do.call(Sys.setenv, as.list(env))
  on.exit({
    mirai::daemons(0, .compute = profile)
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(old[[name]]), name))
    }
  }, add = TRUE)
  mirai::daemons(workers, .compute = profile)
  setup <- list(libPaths = .libPaths(), dev = devPackagePath(), outputFolder = outputFolder, database = database,
                memoryLimit = workerMemoryLimit(workers), rngKind = RNGkind())
  report(sprintf("%d workers, Andromeda memory limit %.2f GB each", workers, setup$memoryLimit))
  stepsFile <- .timing$file
  submit <- function(unit, tries = 1) {
    # the unit goes serialized: a worker loads the package (and its library paths) before it reads the unit
    job <- mirai::mirai({
      .libPaths(setup$libPaths)
      if (!is.null(setup$dev) && !isNamespaceLoaded("SignalAdaptiveStudy")) {
        getExportedValue("pkgload", "load_all")(setup$dev, quiet = TRUE)
      }
      get("workerUnit", envir = asNamespace("SignalAdaptiveStudy"))(payload, setup)
    }, payload = serialize(unit, NULL), setup = setup, .compute = profile)
    list(unit = unit, tries = tries, job = job)
  }
  pending <- list()
  finished <- 0
  queued <- 0
  handle <- function(p) {
    value <- p$job$data
    label <- sprintf("%s k%s o%d", p$unit$task$taskId, p$unit$budget, p$unit$ordering)
    if (mirai::is_mirai_error(value)) { # an error outside the runs, e.g. loading the data
      report(failUnit(p$unit, trimws(as.character(value))[1]), warn = TRUE)
    } else if (mirai::is_error_value(value)) { # the worker died (19: connection reset)
      report(sprintf("%s: the worker stopped (%s)", label, format(value)), warn = TRUE)
      if (identical(as.integer(value), 19L)) mirai::launch_local(1, .compute = profile)
      if (p$tries < 2) {
        pending[[length(pending) + 1]] <<- submit(p$unit, p$tries + 1)
        return(invisible())
      }
      report(failUnit(p$unit, "the worker stopped twice (out of memory?)"), warn = TRUE)
    } else {
      appendSteps(value$steps, stepsFile)
      for (failure in value$failures) report(failure, warn = TRUE)
    }
    finished <<- finished + 1
    report(sprintf("%s done (%d of %d queued units)", label, finished, queued))
  }
  collect <- function(wait) {
    repeat {
      ready <- vapply(pending, function(p) !mirai::unresolved(p$job), logical(1))
      done <- pending[ready]
      pending <<- pending[!ready]
      for (p in done) handle(p)
      if (!wait || length(pending) == 0) return(invisible())
      if (length(done) == 0) Sys.sleep(1)
    }
  }
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    units <- unitsOf(task)
    if (length(units) == 0) next
    folder <- tryCatch(
      if (!is.null(dataSource$plpDataFolder)) dataSource$plpDataFolder(task, outputFolder) else {
        folder <- file.path(tmp, "plpData", task$taskId)
        PatientLevelPrediction::savePlpData(dataSource$loadPlpData(task, outputFolder), folder)
        folder
      }, error = function(e) e)
    if (inherits(folder, "error")) {
      for (unit in units) report(failUnit(unit, conditionMessage(folder)), warn = TRUE)
      next
    }
    for (unit in units) {
      unit$plpDataFolder <- folder
      pending[[length(pending) + 1]] <- submit(unit)
    }
    queued <- queued + length(units)
    collect(wait = FALSE)
  }
  collect(wait = TRUE)
}

#' Run the study
#'
#' `mode = "counts"`: population, training-pool and test counts per task and the feasibility rule (>= 500 training
#' events and >= 80 test events); no models. `mode = "models"`: every task x budget x ordering x method, skipping
#' runs that have already finished.
#'
#' With `workers > 1`, the models run in that many R processes ([mirai](https://mirai.r-lib.org) daemons), one budget
#' and ordering of a task each, the largest first; the results are the same as with one. Each worker uses one thread
#' (duckdb in Andromeda, BLAS, OpenMP) and a memory limit for each Andromeda (half the RAM shared by all workers, or the
#' option `andromedaMemoryLimit`), with temporary files under `outputFolder/tmp`. The data are extracted in this
#' process.
#'
#' @param dataSource A data source, e.g. from [createPlpDataSource()].
#' @param tasks Data frame with one row per task: `taskId`, `outcomeId`, `timeAtRiskDays`; optionally `budgets` (a
#'   list of each task's budgets, as from [addSampleSizes()], instead of `budgets`) and `eventsPerPerson` (to start the
#'   largest units first).
#' @param outputFolder Folder for the data cache, counts and results.
#' @param mode `"counts"` or `"models"`.
#' @param budgets Event budgets k.
#' @param orderings Orderings (1, 2, 3).
#' @param methods Named list of model settings (default [studyMethods()]).
#' @param workers Number of worker processes for the models (1: in this process).
#' @return For counts, a data frame of counts; for models, the result rows of this and earlier runs.
#' @export
execute <- function(dataSource, tasks, outputFolder, mode = c("counts", "models"), budgets = c(100, 250, 500),
                    orderings = 1:3, methods = studyMethods(), workers = 1) {
  mode <- match.arg(mode)
  if (!inherits(dataSource, "sarDataSource")) stop("dataSource must come from a create...DataSource() function")
  if (is.null(names(methods)) || any(names(methods) == "")) stop("methods must be a named list")
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
  database <- dataSource$databaseName
  # timestamped log of this call (appended), and the step timings
  logger <- ParallelLogger::createLogger(name = "SignalAdaptiveStudy", threshold = "TRACE", appenders = list(
    ParallelLogger::createFileAppender(layout = ParallelLogger::layoutTimestamp,
                                       fileName = file.path(outputFolder, sprintf("log_%s.txt", database)))))
  ParallelLogger::unregisterLogger("SignalAdaptiveStudy", silent = TRUE)
  ParallelLogger::registerLogger(logger)
  .timing$file <- file.path(outputFolder, sprintf("steps_%s.csv", database))
  on.exit({
    ParallelLogger::unregisterLogger("SignalAdaptiveStudy", silent = TRUE)
    .timing$file <- NULL
  }, add = TRUE)
  # duckdb >= 2.0 announces where it keeps downloaded extensions and secrets at every new Andromeda (up to 60 times per
  # process, in each worker too) unless the place is chosen; Andromeda uses neither, so a temporary one (inherited by
  # the workers)
  if (is.null(getOption("duckdb.home")) && Sys.getenv("DUCKDB_R_HOME") == "") {
    Sys.setenv(DUCKDB_R_HOME = file.path(tempdir(), "duckdb"))
    on.exit(Sys.unsetenv("DUCKDB_R_HOME"), add = TRUE)
  }
  ParallelLogger::logInfo(sprintf("execute(%s, mode = %s): %d task(s)", database, mode, nrow(tasks)))
  if (mode == "counts") {
    # restartable like the models: a task with counts is skipped; a failed task leaves error.txt and is tried again
    folder <- file.path(outputFolder, "counts", database)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    rows <- lapply(seq_len(nrow(tasks)), function(i) {
      task <- tasks[i, ]
      file <- file.path(folder, paste0(task$taskId, ".csv"))
      fp <- countsFingerprint(dataSource, task)
      fpFile <- file.path(folder, paste0(task$taskId, "_", FINGERPRINT_FILE))
      if (file.exists(file)) {
        changes <- fingerprintChanges(fp, fpFile)
        if (!length(changes)) return(utils::read.csv(file)) # same settings, or counts from before fingerprints
        ParallelLogger::logInfo(sprintf("%s: counts made with other settings (%s), counted again", task$taskId,
                                        paste(changes, collapse = ", ")))
      }
      errorFile <- file.path(folder, paste0(task$taskId, "_error.txt"))
      counts <- tryCatch(timed(task$taskId, "counts", taskCounts(dataSource, task, outputFolder)), error = function(e) {
        writeLines(conditionMessage(e), errorFile)
        message(database, " ", task$taskId, " counts failed: ", conditionMessage(e))
        NULL
      })
      if (is.null(counts)) return(NULL)
      unlink(errorFile)
      row <- data.frame(database = database, taskId = task$taskId, n = counts$n, poolN = counts$poolN,
                        poolEvents = counts$poolEvents, testN = counts$testN, testEvents = counts$testEvents,
                        feasible = counts$poolEvents >= 500 && counts$testEvents >= 80)
      writeFingerprint(fp, fpFile)
      utils::write.csv(row, file, row.names = FALSE)
      row
    })
    return(do.call(rbind, rows))
  }
  clearTaskData()
  on.exit(clearTaskData(), add = TRUE)
  taskBudgets <- function(task) if ("budgets" %in% names(task)) task$budgets[[1]] else budgets
  # finished runs made with other settings: stop before anything runs, rather than mix two designs
  changed <- lapply(seq_len(nrow(tasks)), function(i) changedRuns(tasks[i, ], taskBudgets(tasks[i, ]), orderings, methods,
                                                                  dataSource, outputFolder))
  legacy <- sum(vapply(changed, attr, numeric(1), "legacy"))
  if (legacy > 0) ParallelLogger::logInfo(sprintf("%d finished run(s) from before fingerprints: accepted as they are", legacy))
  changed <- unlist(changed)
  if (length(changed)) {
    stop(sprintf(paste0("%d finished run(s) were made with other settings than now (folder: settings that changed):\n%s%s\n",
                        "Move these run folders away (or restore the settings) and run again."),
                 length(changed), paste(utils::head(changed, 10), collapse = "\n"),
                 if (length(changed) > 10) sprintf("\n... and %d more", length(changed) - 10) else ""), call. = FALSE)
  }
  unitsOf <- function(task) modelUnits(task, taskBudgets(task), orderings, methods, dataSource, outputFolder)
  if (workers > 1) {
    # largest first: the tasks by the training rows of their largest budget, and within a task the largest budget first
    size <- vapply(seq_len(nrow(tasks)), function(i) max(c(0, taskBudgets(tasks[i, ]))) /
                     (if ("eventsPerPerson" %in% names(tasks)) tasks$eventsPerPerson[i] else 1), numeric(1))
    largestFirst <- function(task) {
      units <- unitsOf(task)
      units[order(-vapply(units, function(unit) unit$budget, numeric(1)))]
    }
    runInWorkers(tasks[order(-size), , drop = FALSE], largestFirst, dataSource, outputFolder, workers)
  } else {
    for (i in seq_len(nrow(tasks))) {
      task <- tasks[i, ]
      for (unit in unitsOf(task)) {
        for (failure in runUnit(unit, function() dataSource$loadPlpData(task, outputFolder))) message(failure)
      }
    }
  }
  results <- collectResults(outputFolder)
  convergenceReport(results)
  results
}

# One message on runs whose fits had trouble: a model status other than success, failed fold fits, or fit warnings.
convergenceReport <- function(results) {
  if (!NROW(results) || !"modelStatus" %in% names(results)) return(invisible())
  bad <- !is.na(results$modelStatus) & !results$modelStatus %in% c("SUCCESS", "OK", "INTERCEPT_ONLY")
  folds <- !is.na(results$foldFitsFailed) & results$foldFitsFailed > 0
  warned <- !is.na(results$fitWarnings) & nzchar(results$fitWarnings)
  if (any(bad | folds | warned)) {
    message(sprintf("Convergence: %d run(s) with a model status other than success, %d with failed fold fits, %d with fit warnings. See the modelStatus, foldFitsFailed, foldFitStatus and fitWarnings columns.",
                    sum(bad), sum(folds), sum(warned)))
  } else message("Convergence: every fit converged, no fold fit failed, no fit warnings.")
  invisible()
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
