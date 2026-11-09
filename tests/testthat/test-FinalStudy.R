test_that("the benchmark tasks follow PLPBenchmarks' risk windows and first-entry rule", {
  tasks <- benchmarkTasks()
  expect_equal(nrow(tasks), 13)
  get <- function(name) tasks[tasks$benchmark == name, ]
  expect_equal(c(get("MACEMNCS")$riskWindowStart, get("MACEMNCS_Cr")$riskWindowStart), c(0, 0))
  expect_true(all(tasks$riskWindowStart[!grepl("MACE", tasks$benchmark)] == 1))
  expect_false(any(tasks$firstExposureOnly[tasks$benchmark %in% c("coverF", "coverH", "coverI", "LCOV", "DOV")]))
  expect_true(all(tasks$firstExposureOnly[!tasks$benchmark %in% c("coverF", "coverH", "coverI", "LCOV", "DOV")]))
  s <- taskPopulationSettings(get("MACEMNCS"))
  expect_equal(s$riskWindowStart, 0); expect_equal(s$riskWindowEnd, 30); expect_true(s$firstExposureOnly)
  expect_false(s$requireTimeAtRisk); expect_true(s$includeAllOutcomes); expect_equal(s$washoutPeriod, 365)
  s <- taskPopulationSettings(get("coverF"))
  expect_equal(s$riskWindowStart, 1); expect_false(s$firstExposureOnly)
})

test_that("catalogue tasks keep their earlier population settings", {
  task <- catalogueCohortTasks()[1, ]
  expect_identical(unclass(taskPopulationSettings(task)), unclass(benchmarkPopulationSettings(task$timeAtRiskDays, 1, TRUE)))
})

test_that("the benchmark run has all its cohorts, and the old study set is unchanged", {
  expect_equal(nrow(studyCohortDefinitionSet()), 26)
  set <- benchmarkCohortDefinitionSet()
  tasks <- benchmarkTasks()
  expect_true(all(c(tasks$targetId, tasks$outcomeId) %in% set$cohortId))
  expect_equal(sum(set$use == "feature"), 64)
  expect_equal(sort(set$cohortId[set$use == "feature"]), 100000 + 1152:1215)
  expect_true(all(nchar(set$sql) > 100) && all(nchar(set$json) > 100))
})

test_that("the covariate window ends before the risk window starts, and at day -1 on benchmarks", {
  expect_equal(covariateEndDay(benchmarkTasks()[1, ]), -1)
  expect_equal(covariateEndDay(benchmarkTasks()[6, ]), -1)
  expect_equal(covariateEndDay(list(populationRule = "other", riskWindowStart = 1)), 0)
  expect_equal(covariateEndDay(list(populationRule = "other", riskWindowStart = 90)), 0)
  expect_equal(covariateEndDay(list(populationRule = "other", riskWindowStart = 0)), -1)
  expect_equal(covariateEndDay(list(populationRule = "other")), 0)
})

test_that("the default feature set is FeatureExtraction's defaults without the scores", {
  f <- finalCovariateSettings("s", "t")
  s <- f(expandFeatureSets(benchmarkTasks()[1, ], "defaults"))
  d <- FeatureExtraction::createDefaultCovariateSettings()
  on <- function(x) sort(names(x)[vapply(x, isTRUE, logical(1))])
  expect_equal(on(s), setdiff(on(d), c("CharlsonIndex", "Dcsi", "Chads2", "Chads2Vasc")))
  for (keep in c("DemographicsAgeGroup", "DemographicsGender", "DemographicsRace", "DemographicsEthnicity", "DemographicsIndexYear", "DemographicsIndexMonth")) {
    expect_true(isTRUE(s[[keep]]))
  }
  expect_false(isTRUE(s[["DemographicsAge"]]))
  expect_equal(s$endDays, -1)
  expect_equal(c(s$longTermStartDays, s$shortTermStartDays), c(d$longTermStartDays, d$shortTermStartDays))
  expect_equal(f(list(featureSet = "defaults", populationRule = "other", riskWindowStart = 1))$endDays, 0)
})

test_that("the phenotype feature set has gender, age groups and 51 binary cohort indicators, the antibiotics combined", {
  s <- finalCovariateSettings("schema", "table")(expandFeatureSets(benchmarkTasks()[1, ], "phenotypes"))
  expect_length(s, 2)
  expect_true(isTRUE(s[[1]]$DemographicsGender)); expect_true(isTRUE(s[[1]]$DemographicsAgeGroup)); expect_false(isTRUE(s[[1]][["DemographicsAge"]]))
  expect_equal(nrow(s[[2]]$covariateCohorts), 51)
  expect_false(any(s[[2]]$covariateCohorts$cohortId %in% (100000 + 1201:1214)))
  expect_true(199001 %in% s[[2]]$covariateCohorts$cohortId)
  expect_equal(s[[2]]$valueType, "binary"); expect_equal(s[[2]]$startDay, -365); expect_equal(s[[2]]$endDay, -1)
  expect_equal(s[[2]]$covariateCohortTable, "table")
})

test_that("a task runs once per feature set with its own settings", {
  tasks <- expandFeatureSets(benchmarkTasks())
  expect_equal(nrow(tasks), 26)
  expect_false(anyDuplicated(tasks$taskId) > 0)
  expect_equal(sum(tasks$featureSet == "phenotypes"), 13)
})

test_that("the gate keeps tasks with 1000 pool and 333 test events", {
  counts <- data.frame(taskId = c("a", "b", "c"), poolEvents = c(1000, 999, 5000), testEvents = c(333, 900, 332))
  expect_equal(feasibleCounts(counts)$taskId, "a")
})

test_that("calibration and the AUC's standard error match their definitions", {
  set.seed(3)
  n <- 4000
  eta <- rnorm(n)
  y <- rbinom(n, 1, plogis(-1 + 2 * eta))
  p <- plogis(eta)
  m <- calibrationMetrics(y, p)
  g <- glm(y ~ eta, family = binomial())
  expect_equal(m$calibrationSlope, unname(coef(g)[2]), tolerance = 1e-6)
  expect_equal(m$calibrationIntercept, unname(coef(g)[1]), tolerance = 1e-6)
  large <- glm(y ~ 1 + offset(eta), family = binomial())
  expect_equal(m$calibrationInTheLarge, unname(coef(large)[1]), tolerance = 1e-6)
  expect_equal(m$expectedOverObserved, sum(p) / sum(y))
  # DeLong's variance against its O(n1 x n0) definition
  s <- sample(n, 400); ys <- y[s]; ps <- p[s]
  pos <- ps[ys == 1]; neg <- ps[ys == 0]
  psi <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
  byHand <- sqrt(var(rowMeans(psi)) / length(pos) + var(colMeans(psi)) / length(neg))
  expect_equal(aucDeLongSe(ys, ps), byHand, tolerance = 1e-10)
  expect_true(all(is.na(unlist(calibrationMetrics(rep(0, 10), rep(0.1, 10))[c("calibrationSlope", "calibrationInTheLarge")]))))
  constant <- calibrationMetrics(c(0, 1, 0, 1), rep(0.5, 4))
  expect_true(is.na(constant$calibrationSlope)); expect_equal(constant$calibrationInTheLarge, 0, tolerance = 1e-8)
})

test_that("the follow-up metrics use only rows with the full time at risk", {
  population <- data.frame(rowId = 1:6, timeAtRisk = c(365, 365, 365, 100, 365, 10))
  test <- data.frame(rowId = 1:6, outcomeCount = c(1, 0, 1, 1, 0, 0), value = c(0.8, 0.2, 0.6, 0.1, 0.4, 0.9))
  m <- followUpMetrics(test, population, list(riskWindowStart = 1, riskWindowEnd = 365), c(0, 1))
  expect_equal(m$nTestFullFollowUp, 4); expect_equal(m$eventsTestFullFollowUp, 2); expect_equal(m$aucFullFollowUp, 1)
  none <- followUpMetrics(test[1, ], population, list(riskWindowStart = 1, riskWindowEnd = 365), c(0, 1))
  expect_true(is.na(none$aucFullFollowUp))
})

test_that("every method of the final study runs end to end, saves, and gives the new metrics", {
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
  folder <- tempfile("sar_final_")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  methods <- finalStudyMethods()
  expect_setequal(names(methods), c("SAR", "LASSO", "Ridge", "PLP_LASSO", "ElasticNet", "SAR_fixed", "SAR_no_c3"))
  task$budgets <- list(100)
  results <- execute(source, task, folder, mode = "models", orderings = 1, methods = methods, workers = 1)
  expect_equal(nrow(results), length(methods))
  expect_setequal(results$method, names(methods))
  expect_true(all(results$auc > 0.5 & results$auc < 1))
  expect_true(all(results$featureSet == "defaults"))
  for (column in c("aucSe", "calibrationSlope", "calibrationInTheLarge", "expectedOverObserved", "nTestFullFollowUp", "aucFullFollowUp")) expect_false(anyNA(results[[column]]), info = column)
  expect_true(all(results$nTestFullFollowUp <= results$nTest))
  for (column in c("modelStatus", "foldFitsFailed", "foldFitStatus", "fitWarnings", "finalIterations", "finalAlgorithm", "finalNewtonIterations",
                   "finalGapEstimate", "foldFitIterationsTotal", "foldFitIterationsMax", "foldNewtonIterations", "foldFitAlgorithms")) expect_true(column %in% names(results), info = column)
  # time inside Cyclops is part of the fit, which is part of the run; PLP's LASSO does not report it
  study <- results[results$method != "PLP_LASSO", ]
  expect_true(all(study$cyclopsFits > 0 & study$cyclopsFitSeconds > 0 & study$cyclopsDataSeconds >= 0))
  expect_true(all(study$cyclopsFitSeconds + study$cyclopsDataSeconds <= study$fitSeconds + 1))
  expect_true(all(results$predictSeconds > 0 & results$saveSeconds > 0))
  expect_true(all(results$fitSeconds + results$predictSeconds + results$saveSeconds <= results$runSeconds + 1))
  expect_gt(study$cyclopsFits[study$method == "SAR"], study$cyclopsFits[study$method == "LASSO"])
  zipFile <- exportResults(folder, minCellCount = 5)
  expect_true(any(grepl("^steps_SIM\\.csv$", utils::unzip(zipFile, list = TRUE)$Name)))
  expect_true(all(results$modelStatus %in% c("SUCCESS", "OK", "INTERCEPT_ONLY", "POOR_BLR_STEP")))
  expect_true(all(results$foldFitsFailed[results$method != "PLP_LASSO"] >= 0))
  # the same training set for every method
  expect_equal(length(unique(results$nTrain)), 1)
  # each saved model predicts the same as when it was fitted
  elastic <- results[results$method == "ElasticNet", ]
  expect_gt(elastic$nonZero, 0)
  expect_equal(elastic$search, "auto")
  expect_lte(elastic$kkt, 1e-3)
  expect_equal(results$c1[results$method == "SAR_fixed"], 0)
  expect_equal(results$c2[results$method == "SAR_fixed"], 1)
  expect_equal(results$c3[results$method == "SAR_fixed"], -1)
  expect_equal(results$c3[results$method == "SAR_no_c3"], 0)
  # a restart skips finished runs
  again <- execute(source, task, folder, mode = "models", orderings = 1, methods = methods, workers = 1)
  expect_identical(sort(again$auc), sort(results$auc))
})

test_that("the elastic net's fit and search agree with glmnet and its definition", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  d <- simulateDesign()
  alpha <- ELASTIC_NET_ALPHA
  rows <- rep(TRUE, length(d$y))
  problem <- createFitProblem(d$X, d$y, rows)
  sigma2 <- 0.05
  v <- rep(sigma2, length(problem$ids))
  fit <- fitProblem(problem, "elasticnet", v, rep(TRUE, length(v)), NULL, Cyclops::createControl(fitter = "super", convergenceType = "kkt", maxIterations = 100000, threads = 1, tolerance = 1e-6))
  expect_true(fit$converged)
  lambda <- varianceToLambda("elasticnet", sigma2, problem$n)
  ref <- glmnet::glmnet(as.matrix(d$X[, problem$ids]), d$y, family = "binomial", alpha = alpha, lambda = lambda, standardize = FALSE,
                        control = list(thresh = 1e-14, maxit = 1e7))
  expect_equal(fit$coefficients[-1], as.numeric(ref$beta), tolerance = 1e-3)
  expect_equal(fit$coefficients[1], as.numeric(ref$a0), tolerance = 1e-3)
  expect_lte(kktResidual(d$X, d$y, problem, fit$coefficients, v, rep(TRUE, length(v)), "elasticnet"), 1e-3)
  # at the null-model bound every penalised coefficient is 0, and just above it some is not
  setup <- problemWeights(d$X, d$y, problem, c(0, 0, 0))
  control <- Cyclops::createControl(fitter = "super", maxIterations = 100000, threads = 1)
  bound <- nullModelVariance(d$X, d$y, problem, setup, control, alpha = alpha)
  at <- fitProblem(problem, "elasticnet", rep(bound * 0.999, length(v)), rep(TRUE, length(v)), NULL, control)
  above <- fitProblem(problem, "elasticnet", rep(bound * 1.5, length(v)), rep(TRUE, length(v)), NULL, control)
  expect_true(all(at$coefficients[-1] == 0))
  expect_true(any(above$coefficients[-1] != 0))
  # the whole search, on the study's folds
  cv <- signalAdaptiveCv(d$X, d$y, d$folds, penalty = "elasticnet", shapes = c(0, 0, 0))
  expect_equal(cv$penalty, "elasticnet")
  expect_gt(sum(cv$beta != 0), 0)
  expect_lte(cv$kkt, 1e-3)
  expect_equal(setElasticNet()$settings$penalty, "elasticnet")
})

test_that("the phenotype features and benchmark population run through the real extraction (Eunomia)", {
  skip_if_not_installed("Eunomia")
  skip_on_cran()
  details <- tryCatch(Eunomia::getEunomiaConnectionDetails(databaseFile = tempfile(fileext = ".duckdb"), dbms = "duckdb"),
                      error = function(e) skip(paste("Eunomia not available:", conditionMessage(e))))
  outputFolder <- tempfile("sar_eunomia_final_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  features <- benchmarkCohortDefinitionSet()
  features <- features[features$use == "feature", ]
  tasks <- benchmarkTasks()[1, ]
  tasks$taskId <- "eunomia_task"
  tasks$targetId <- features$cohortId[grepl("Osteoarthritis", features$cohortName)][1]
  tasks$outcomeId <- features$cohortId[grepl("gastrointestinal bleeding", features$cohortName, ignore.case = TRUE)][1]
  source <- createCdmDataSource(details, cdmDatabaseSchema = "main", cohortDatabaseSchema = "main", databaseName = "EUNOMIA",
                                covariateSettings = finalCovariateSettings("main", "sar_study_cohort"), cohortDefinitionSet = features,
                                incremental = TRUE, combinedCohorts = combinedFeatureCohorts())
  # the combined cohort is the union of the 14 antibiotic classes, each episode once
  invisible(suppressWarnings(source$loadPlpData(expandFeatureSets(tasks, "phenotypes"), outputFolder)))
  connection <- DatabaseConnector::connect(details)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
  q <- function(sql) DatabaseConnector::querySql(connection, sql)[[1]]
  union <- q(sprintf("SELECT COUNT(*) FROM (SELECT DISTINCT subject_id, cohort_start_date, cohort_end_date FROM main.sar_study_cohort WHERE cohort_definition_id IN (%s)) x",
                     paste(100000 + 1201:1214, collapse = ",")))
  combined <- q("SELECT COUNT(*) FROM main.sar_study_cohort WHERE cohort_definition_id = 199001")
  expect_gt(union, 0)
  expect_equal(combined, union)
  for (set in c("defaults", "phenotypes")) {
    task <- expandFeatureSets(tasks, set)
    plpData <- suppressWarnings(source$loadPlpData(task, outputFolder))
    ref <- dplyr::collect(plpData$covariateData$covariateRef)
    values <- dplyr::collect(plpData$covariateData$covariates)
    expect_true(all(values$covariateValue == 1), info = set) # every covariate binary
    expect_true(any(grepl("age group", ref$covariateName)), info = set)
    expect_false(any(ref$covariateId == 1002), info = set) # no continuous age
    expect_false(any(grepl("Charlson|DCSI|CHADS", ref$covariateName, ignore.case = TRUE)), info = set)
    if (set == "phenotypes") {
      cohorts <- ref$covariateId[ref$analysisId == 456]
      expect_gt(length(cohorts), 0)
      ids <- (cohorts - 456) / 1000
      expect_true(all(ids %in% c(100000 + setdiff(1152:1215, 1201:1214), 199001)))
      expect_true(all(ref$analysisId %in% c(1, 3, 456)))
    } else {
      expect_gt(length(unique(ref$analysisId)), 3)
    }
  }
})

test_that("the two fitting algorithms reach the same fit, and report how", {
  skip_on_cran()
  d <- simulateDesign()
  ccd <- signalAdaptiveCv(d$X, d$y, d$folds, penalty = "lasso", shapes = c(0, 0, 0), searchSettings = list(fitAlgorithm = "ccd"))
  auto <- signalAdaptiveCv(d$X, d$y, d$folds, penalty = "lasso", shapes = c(0, 0, 0), searchSettings = list(fitAlgorithm = "auto"))
  expect_equal(ccd$convergence$finalAlgorithm, "ccd")
  expect_true(auto$convergence$finalAlgorithm %in% c("newton", "ccd", "auto"))
  expect_gt(ccd$convergence$foldFitIterationsTotal, 0)
  expect_gte(ccd$convergence$foldFitIterationsMax, 1)
  expect_equal(ccd$convergence$foldNewtonIterations, 0)
  expect_true(grepl("ccd", ccd$convergence$foldFitAlgorithms))
  expect_lt(abs(ccd$selected$cvLoss - auto$selected$cvLoss), 1e-3)
  expect_error(signalAdaptiveCv(d$X, d$y, d$folds, penalty = "lasso", shapes = c(0, 0, 0), searchSettings = list(fitAlgorithm = "x")), "unknown fit algorithm")
  methods <- fitAlgorithmComparisonMethods()
  expect_length(methods, 10)
  expect_equal(methods$LASSO_ccd$param$searchSettings$fitAlgorithm, "ccd")
  expect_equal(methods$SAR_fixed_auto$param$searchSettings$fitAlgorithm, "auto")
})

test_that("failed fits are counted, summarised and reported", {
  expect_equal(statusSummary(c("MAX_ITERATIONS", "MAX_ITERATIONS", "ILLCONDITIONED")), "ILLCONDITIONED x 1; MAX_ITERATIONS x 2")
  expect_equal(statusSummary(c("MAX_ITERATIONS x 2", "MAX_ITERATIONS", "")), "MAX_ITERATIONS x 3")
  expect_equal(statusSummary(character()), "")
  # a fit stopped by the iteration limit is not converged, and says why
  d <- simulateDesign()
  problem <- createFitProblem(d$X, d$y, rep(TRUE, length(d$y)))
  fit <- fitProblem(problem, "lasso", rep(0.5, length(problem$ids)), rep(TRUE, length(problem$ids)), NULL,
                    Cyclops::createControl(fitter = "super", convergenceType = "kkt", tolerance = 1e-12, maxIterations = 1, threads = 1))
  expect_false(fit$converged)
  expect_equal(fit$status, "MAX_ITERATIONS")
  results <- data.frame(modelStatus = c("SUCCESS", "REACHED MAXIMUM NUMBER OF ITERATIONS, CANNOT FIT"), foldFitsFailed = c(2, NA), fitWarnings = c("", ""))
  expect_message(convergenceReport(results), "1 run\\(s\\) with a model status other than success, 1 with failed fold fits, 0 with fit warnings")
})
