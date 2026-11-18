# Aggregates of one training set and of every Cyclops fit of chosen models on it, for building simulated designs that
# behave like the real data (extras/production/fitTrace.R). Nothing per patient: counts of columns, quantiles over rows
# and columns, and one row per fit (its size, prior variance, iterations, seconds). A single person's value (the most
# covariates in a row) is replaced by the mean of the top minCellCount rows, and counts of persons below minCellCount are
# shown as -minCellCount.

# The design: rows, events, folds, columns, non-zeros per row and per column, duplicate columns, column types.
designSummary <- function(X, y, folds, covariateIds, covariateData, minCellCount = 5, nearDuplicateColumns = 3000) {
  X <- methods::as(methods::as(X, "CsparseMatrix"), "generalMatrix")
  n <- nrow(X); p <- ncol(X)
  perRow <- Matrix::rowSums(X != 0)
  perColumn <- diff(X@p)
  prevalence <- perColumn / n
  q <- function(v, probs) stats::setNames(as.numeric(stats::quantile(v, probs, names = FALSE)), paste0("p", probs * 100))
  topMean <- function(v) mean(sort(v, decreasing = TRUE)[seq_len(min(length(v), minCellCount))])
  foldSizes <- as.integer(table(folds))
  design <- data.frame(
    rows = n, events = sum(y), folds = length(foldSizes), foldRowsMin = min(foldSizes), foldRowsMax = max(foldSizes),
    fitRowsPerFold = n - stats::median(foldSizes), columns = p, nonZeros = length(X@i),
    allBinary = all(X@x == 1),
    rowNonZeroMean = mean(perRow), rowNonZeroMedian = stats::median(perRow),
    rowNonZeroP90 = q(perRow, 0.9), rowNonZeroP99 = q(perRow, 0.99), rowNonZeroTopMean = topMean(perRow),
    columnRowsP10 = q(perColumn, 0.1), columnRowsP50 = q(perColumn, 0.5), columnRowsP90 = q(perColumn, 0.9),
    columnRowsP99 = q(perColumn, 0.99), columnRowsMax = max(perColumn),
    columnsBelowMinCell = sum(perColumn < minCellCount),
    columnsAbove1pct = sum(prevalence > 0.01), columnsAbove10pct = sum(prevalence > 0.1), columnsAbove50pct = sum(prevalence > 0.5)
  )

  # exact duplicate columns (same rows and values): a cheap key, then checked exactly within each key
  j <- rep(seq_len(p), perColumn)
  sums <- rowsum(cbind(X@i + 1, (X@i + 1)^2 %% 1000003, X@x), j, reorder = TRUE)
  key <- rep(NA_character_, p)
  key[as.integer(rownames(sums))] <- paste(perColumn[as.integer(rownames(sums))], sums[, 1], sums[, 2], sums[, 3])
  groups <- split(seq_len(p), key)
  groups <- groups[lengths(groups) > 1]
  rowsOf <- function(k) X@i[(X@p[k] + 1):X@p[k + 1]]
  duplicates <- 0L; duplicateGroups <- 0L
  for (g in groups) {
    signature <- vapply(g, function(k) paste(rowsOf(k), collapse = ","), character(1))
    tab <- table(signature)
    duplicates <- duplicates + sum(tab[tab > 1] - 1L); duplicateGroups <- duplicateGroups + sum(tab > 1)
  }
  design$exactDuplicateColumns <- duplicates # columns that repeat an earlier one
  design$exactDuplicateGroups <- duplicateGroups

  # near-duplicates among the most prevalent columns: correlation of the 0/1 columns from their overlaps
  top <- order(perColumn, decreasing = TRUE)[seq_len(min(p, nearDuplicateColumns))]
  B <- X[, top, drop = FALSE]; B@x[] <- 1
  overlap <- as.matrix(Matrix::crossprod(B))
  a <- diag(overlap)
  correlation <- (n * overlap - outer(a, a)) / sqrt(outer(a * (n - a), a * (n - a)))
  upper <- correlation[upper.tri(correlation)]
  design$nearDuplicateColumnsChecked <- length(top)
  design$nearDuplicateColumnsCheckedMinRows <- min(perColumn[top])
  design$pairsCorrelationAbove0.95 <- sum(upper > 0.95 & upper < 1 - 1e-12, na.rm = TRUE)
  design$pairsCorrelationAbove0.8 <- sum(upper > 0.8 & upper < 1 - 1e-12, na.rm = TRUE)
  design$pairsCorrelationExact1 <- sum(upper >= 1 - 1e-12, na.rm = TRUE)

  # column types: FeatureExtraction's analyses (names of analyses, not of covariates)
  ref <- as.data.frame(dplyr::collect(covariateData$covariateRef))
  analyses <- tryCatch(as.data.frame(dplyr::collect(covariateData$analysisRef)), error = function(e) NULL)
  analysis <- ref$analysisId[match(covariateIds, ref$covariateId)]
  name <- as.character(analysis)
  if (!is.null(analyses) && all(c("analysisId", "analysisName") %in% names(analyses))) {
    m <- match(analysis, analyses$analysisId)
    name[!is.na(m)] <- as.character(analyses$analysisName[m[!is.na(m)]])
  }
  name[is.na(name)] <- "unknown"
  types <- do.call(rbind, lapply(split(seq_len(p), name), function(k) data.frame(
    analysis = name[k[1]], columns = length(k), nonZeroShare = sum(perColumn[k]) / length(X@i),
    rowsMedian = stats::median(perColumn[k]), rowsP90 = q(perColumn[k], 0.9), rowsMax = max(perColumn[k]),
    columnsAbove1pct = sum(prevalence[k] > 0.01), rowsWithAny = sum(Matrix::rowSums(X[, k, drop = FALSE] != 0) > 0))))
  types <- types[order(-types$nonZeroShare), ]
  rownames(types) <- NULL
  # counts of persons below minCellCount shown as -minCellCount, as exportResults() does
  censor <- function(v) ifelse(v > 0 & v < minCellCount, -minCellCount, v)
  for (k in c("columnRowsP10", "columnRowsP50", "columnRowsP90", "columnRowsP99", "columnRowsMax", "nearDuplicateColumnsCheckedMinRows")) {
    design[[k]] <- censor(design[[k]])
  }
  for (k in c("rowsMedian", "rowsP90", "rowsMax", "rowsWithAny")) types[[k]] <- censor(types[[k]])
  list(design = design, columnTypes = types)
}

# The training set of one task, budget and ordering, as the study's runs see it (execute()'s split and preprocessing).
# minFraction: PLP's frequency filter (covariates in fewer than minFraction x rows removed, as FeatureExtraction's
# tidyCovariateData()) instead of the study's; NULL: the study's preprocessing.
studyTrainingSet <- function(dataSource, task, budget, ordering, outputFolder, minFraction = NULL) {
  plpData <- dataSource$loadPlpData(task, outputFolder)
  population <- PatientLevelPrediction::createStudyPopulation(plpData, task$outcomeId, dataSource$populationSettings(task))
  split <- dataSource$split(population, task, budget, ordering)
  preprocess <- studyPreprocessSettings()
  if (!is.null(minFraction)) preprocess$minFraction <- minFraction
  train <- prepareRunData(plpData, population, split, preprocess)$Train
  sparse <- PatientLevelPrediction::toSparseM(train, train$labels)
  labels <- sparse$labels
  list(train = train, X = sparse$dataMatrix, y = as.numeric(labels$outcomeCount > 0),
       folds = train$folds$index[match(labels$originalRowId, train$folds$rowId)],
       covariateIds = sparse$covariateMap$covariateId[match(seq_len(ncol(sparse$dataMatrix)), sparse$covariateMap$columnId)])
}

# Fits the given models (as fitPlp() would) with every Cyclops fit traced, and summarises the trace per model.
traceModels <- function(train, methods, folder) {
  summaries <- list()
  for (name in names(methods)) {
    file <- file.path(folder, sprintf("trace_%s.csv", name))
    unlink(file)
    old <- options(SignalAdaptiveStudy.fitTrace = file)
    started <- Sys.time()
    model <- tryCatch(fitSignalAdaptive(train, methods[[name]]), error = function(e) e)
    options(old)
    if (inherits(model, "error")) { message(name, " failed: ", conditionMessage(model)); next }
    trace <- utils::read.csv(file, stringsAsFactors = FALSE)
    trace$method <- name
    utils::write.csv(trace, file, row.names = FALSE)
    summaries[[name]] <- traceSummary(trace, name, as.numeric(difftime(Sys.time(), started, units = "secs")))
  }
  do.call(rbind, summaries)
}

traceSummary <- function(trace, name, wallSeconds) {
  fold <- trace[trace$stage == "fold", ]; final <- trace[trace$stage == "final", ]
  q <- function(v, pr) as.numeric(stats::quantile(v, pr, names = FALSE, na.rm = TRUE))
  perShape <- table(fold$shape) / max(1, length(unique(fold$rows))) # variances tried per shape (one fit per fold each)
  smallest <- fold[fold$varianceMedian == min(fold$varianceMedian), ]; largest <- fold[fold$varianceMedian == max(fold$varianceMedian), ]
  data.frame(
    method = name, wallSeconds = wallSeconds, shapes = length(unique(fold$shape)),
    variancesTriedPerShapeMedian = stats::median(as.numeric(perShape)), foldFits = nrow(fold),
    foldSeconds = sum(fold$seconds), finalSeconds = sum(final$seconds),
    foldIterationsP10 = q(fold$iterations, 0.1), foldIterationsP50 = q(fold$iterations, 0.5),
    foldIterationsP90 = q(fold$iterations, 0.9), foldIterationsMax = max(fold$iterations),
    foldSecondsP10 = q(fold$seconds, 0.1), foldSecondsP50 = q(fold$seconds, 0.5), foldSecondsP90 = q(fold$seconds, 0.9),
    foldSecondsMax = max(fold$seconds), foldNewtonSteps = sum(fold$newton, na.rm = TRUE),
    foldNotSuccess = sum(fold$status != "SUCCESS"),
    varianceMin = min(fold$varianceMedian), varianceMax = max(fold$varianceMedian),
    nonZeroAtSmallestVariance = stats::median(smallest$nonZero), nonZeroAtLargestVariance = stats::median(largest$nonZero),
    selectedShape = if (nrow(final)) final$shape[1] else NA_character_,
    selectedVarianceMedian = if (nrow(final)) final$varianceMedian[1] else NA_real_,
    selectedVarianceSpread = if (nrow(final)) final$varianceSpread[1] else NA_real_,
    finalNonZero = if (nrow(final)) final$nonZero[1] else NA_integer_,
    finalIterations = if (nrow(final)) final$iterations[1] else NA_integer_,
    finalNewton = if (nrow(final)) final$newton[1] else NA_integer_,
    finalSolver = if (nrow(final)) final$solver[1] else NA_character_,
    finalCertificate = if (nrow(final)) final$certificate[1] else NA_character_,
    finalStatus = if (nrow(final)) final$status[1] else NA_character_)
}

# The whole report for one task: design aggregates, traced fits, a zip to send back. minFraction: all of it with that
# frequency filter instead of the study's (see studyTrainingSet()); the zip's name says so.
fitTraceReport <- function(dataSource, task, budget, methods, outputFolder, reportFolder, ordering = 1, minCellCount = 5,
                           minFraction = NULL) {
  dir.create(reportFolder, recursive = TRUE, showWarnings = FALSE)
  message("Preparing the training set: ", task$taskId, ", k", budget, ", ordering ", ordering,
          sprintf(", covariates in at least %g of the rows", minFraction %||% studyPreprocessSettings()$minFraction))
  set <- studyTrainingSet(dataSource, task, budget, ordering, outputFolder, minFraction)
  message("Design aggregates")
  d <- designSummary(set$X, set$y, set$folds, set$covariateIds, set$train$covariateData, minCellCount = minCellCount)
  d$design <- cbind(database = dataSource$databaseName, taskId = task$taskId, budget = budget, ordering = ordering,
                    minFraction = minFraction %||% studyPreprocessSettings()$minFraction, d$design)
  utils::write.csv(d$design, file.path(reportFolder, "design.csv"), row.names = FALSE)
  utils::write.csv(d$columnTypes, file.path(reportFolder, "column_types.csv"), row.names = FALSE)
  print(t(d$design)); print(d$columnTypes, row.names = FALSE)
  summary <- if (length(methods)) traceModels(set$train, methods, reportFolder)
  if (!is.null(summary)) {
    utils::write.csv(summary, file.path(reportFolder, "trace_summary.csv"), row.names = FALSE)
    print(t(summary))
  }
  writeLines(c(paste("Exported", format(Sys.time(), usetz = TRUE)), packageVersions(), paste("minCellCount", minCellCount)),
             file.path(reportFolder, "export.txt"))
  zipFile <- normalizePath(file.path(reportFolder, sprintf("SignalAdaptiveStudy_fit_trace_%s%s_%s.zip", dataSource$databaseName,
                                                          if (!is.null(minFraction)) sprintf("_minFraction%g", minFraction) else "",
                                                          format(Sys.Date(), "%Y%m%d"))), mustWork = FALSE)
  unlink(zipFile)
  zip::zipr(zipFile, list.files(reportFolder, "\\.(csv|txt)$", full.names = TRUE))
  zipFile
}
