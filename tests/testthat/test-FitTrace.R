test_that("designSummary counts duplicates, near-duplicates and column types", {
  X <- Matrix::sparseMatrix(i = c(1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 6), j = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5), x = 1,
                            dims = c(10, 5))
  covariateData <- Andromeda::andromeda(covariateRef = data.frame(covariateId = 11:15, analysisId = c(1, 1, 2, 2, 2)),
                                        analysisRef = data.frame(analysisId = c(1, 2), analysisName = c("A", "B")))
  on.exit(Andromeda::close(covariateData))
  d <- designSummary(X, c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1), rep(1:2, 5), 11:15, covariateData, minCellCount = 2)
  expect_equal(d$design$exactDuplicateColumns, 1)      # columns 1 and 2
  expect_equal(d$design$pairsCorrelationExact1, 1)
  expect_equal(d$design$nonZeros, 13)
  expect_equal(d$design$rowNonZeroTopMean, 3)         # rows 1-3 have 3 each
  expect_equal(d$columnTypes$columns[d$columnTypes$analysis == "B"], 3)
  expect_equal(sum(d$columnTypes$nonZeroShare), 1)
  d5 <- designSummary(X, c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1), rep(1:2, 5), 11:15, covariateData, minCellCount = 5)
  expect_equal(d5$design$columnRowsP10, -5)             # 1-4 persons per column: censored
  expect_true(all(d5$columnTypes$rowsMedian == -5))
})

test_that("fitTraceReport traces every fit of a model on the study's training set", {
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
  source <- createPlpDataSource(plpData, "SIM", settings)
  task <- data.frame(taskId = "sim_task_defaults", outcomeId = 3, timeAtRiskDays = 365, featureSet = "defaults")
  folder <- tempfile("sar_trace_")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  methods <- fitAlgorithmComparisonMethods()[c("Ridge_ccd", "SAR_fixed_auto")]
  zipFile <- fitTraceReport(source, task, 100, methods, folder, folder)
  expect_true(file.exists(zipFile))
  design <- utils::read.csv(file.path(folder, "design.csv"))
  expect_equal(design$events, 100)
  summary <- utils::read.csv(file.path(folder, "trace_summary.csv"))
  expect_setequal(summary$method, names(methods))
  trace <- utils::read.csv(file.path(folder, "trace_Ridge_ccd.csv"))
  expect_setequal(unique(trace$stage), c("fold", "final"))
  expect_equal(sum(trace$stage == "final"), 1)
  expect_true(all(trace$iterations > 0 & trace$seconds >= 0 & trace$status == "SUCCESS"))
  expect_equal(summary$foldFits[summary$method == "Ridge_ccd"], sum(trace$stage == "fold"))
  expect_setequal(unique(trace$fold[trace$stage == "fold"]), c("1", "2", "3"))
  expect_equal(trace$fold[trace$stage == "final"], "all")
  expect_true(all(trace$start[trace$stage == "fold"] %in% c("cold", "search")))
  sar <- utils::read.csv(file.path(folder, "trace_SAR_fixed_auto.csv"))
  expect_true(all(c("gapEstimate", "varianceMin", "varianceMax", "maxIterations") %in% names(sar)))
  expect_false(file.exists(file.path(folder, "trace_Ridge_ccd_failed.txt"))) # every fit converged
  expect_null(getOption("SignalAdaptiveStudy.fitTrace"))
  # after the frequency filter: fewer columns, each in at least 5% of the rows, and the models fitted on them
  filteredFolder <- tempfile("sar_trace_filtered_")
  on.exit(unlink(filteredFolder, recursive = TRUE), add = TRUE)
  zipFiltered <- fitTraceReport(source, task, 100, methods["Ridge_ccd"], filteredFolder, filteredFolder, minFraction = 0.05)
  expect_match(basename(zipFiltered), "minFraction0.05")
  filtered <- utils::read.csv(file.path(filteredFolder, "design.csv"))
  expect_equal(filtered$rows, design$rows)
  expect_lt(filtered$columns, design$columns)
  expect_gte(filtered$columnRowsP10, 0.05 * filtered$rows)
  filteredTrace <- utils::read.csv(file.path(filteredFolder, "trace_Ridge_ccd.csv"))
  expect_equal(max(filteredTrace$columns), filtered$columns)
  only <- tempfile("sar_trace_design_")
  on.exit(unlink(only, recursive = TRUE), add = TRUE)
  fitTraceReport(source, task, 100, list(), only, only)
  expect_false(any(grepl("^trace", list.files(only))))
})

test_that("a traced fit that does not converge is described and retried with variants", {
  set.seed(5)
  X <- Matrix::rsparsematrix(400, 30, density = 0.2, rand.x = function(n) rep(1, n))
  y <- stats::rbinom(400, 1, stats::plogis(-1 + as.numeric(X[, 1:5] %*% rep(1, 5))))
  problem <- createFitProblem(X, y, rep(TRUE, 400))
  keep <- rep(TRUE, length(problem$ids))
  control <- Cyclops::createControl(fitter = "super", algorithm = "auto", convergenceType = "kkt", tolerance = 1e-8,
                                    maxIterations = 1, threads = 1)
  file <- tempfile("trace_", fileext = ".csv")
  on.exit(unlink(c(file, sub("\\.csv$", "_failed.txt", file))), add = TRUE)
  old <- options(SignalAdaptiveStudy.fitTrace = file, SignalAdaptiveStudy.fitTraceRetry = TRUE)
  on.exit(options(old), add = TRUE)
  .cyclopsTiming$stage <- "fold"; .cyclopsTiming$shape <- "0 0 -2"; .cyclopsTiming$fold <- 2; .cyclopsTiming$start <- "shape"
  fit <- fitProblem(problem, "lasso", rep(0.05, length(problem$ids)), keep, NULL, control)
  resetCyclopsTiming()
  expect_false(fit$converged)
  trace <- utils::read.csv(file, stringsAsFactors = FALSE)
  expect_equal(trace$status[1], "MAX_ITERATIONS")
  expect_equal(as.character(trace$fold[1]), "2")
  expect_equal(trace$stage[-1], paste("fold", c("retry: maxIterations 100000", "retry: cold start", "retry: ccd",
                                               "retry: ccd, maxIterations 100000")))
  expect_equal(trace$maxIterations, c(1, 100000, 1, 1, 100000))
  expect_equal(trace$algorithm[4:5], c("ccd", "ccd"))
  expect_equal(trace$status[2], "SUCCESS")              # with enough iterations it converges
  expect_equal(trace$shape, rep("0 0 -2", 5))
  failed <- readLines(sub("\\.csv$", "_failed.txt", file))
  expect_match(failed[1], "status=MAX_ITERATIONS", fixed = TRUE)
  expect_match(failed[1], "fold=2", fixed = TRUE)
})
