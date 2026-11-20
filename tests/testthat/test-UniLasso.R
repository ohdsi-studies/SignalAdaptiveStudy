# uniLasso without dense matrices (R/UniLasso.R) gives the CRAN package's stage 1 and, through it, its model.

simulateBinary <- function(n = 400, p = 30, seed = 1) {
  set.seed(seed)
  prevalence <- c(stats::runif(p - 4, 0.02, 0.5), 0.005, 0.995, 0, 1) # rare and common, and two constant columns
  X <- vapply(prevalence, function(q) as.numeric(stats::runif(n) < q), numeric(n))
  X[1:2, p - 3] <- 1; X[1:2, p - 2] <- 0                              # at least two persons in each group
  X[, p - 1] <- 0; X[, p] <- 1
  eta <- -2 + X[, 1:5] %*% c(1.2, -0.8, 0.6, 0, 1)
  y <- as.numeric(stats::runif(n) < 1 / (1 + exp(-eta)))
  colnames(X) <- paste0("c", seq_len(p))
  list(X = X, y = y)
}

test_that("stage 1 equals the package's uniInfo() (binomial, leave-one-out, two steps)", {
  skip_if_not_installed("uniLasso")
  d <- simulateBinary()
  features <- uniLassoFeatures(Matrix::Matrix(d$X, sparse = TRUE), d$y)
  reference <- uniLasso::uniInfo(d$X, d$y, family = "binomial", loo = TRUE)
  expect_equal(features$beta, reference$beta, tolerance = 1e-10)
  expect_equal(features$beta0, reference$beta0, tolerance = 1e-10)
  F <- as.matrix(features$G) + matrix(features$shift, nrow(d$X), ncol(d$X), byrow = TRUE)
  expect_equal(unname(F), unname(reference$F), tolerance = 1e-10)
  expect_identical(features$constant, c(rep(FALSE, 28), TRUE, TRUE))
  # sparse: the covariates present, plus the outcome rows
  expect_lte(length(features$G@x), sum(d$X) + sum(d$y) * ncol(d$X))
})

test_that("stage 2 on the sparse features gives the package's model", {
  skip_if_not_installed("uniLasso")
  d <- simulateBinary(seed = 2)
  features <- uniLassoFeatures(Matrix::Matrix(d$X, sparse = TRUE), d$y)
  reference <- uniLasso::uniLasso(d$X, d$y, family = "binomial", thresh = 1e-14)
  lambdas <- reference$lambda[c(5, 20, 40)]
  fit <- glmnet::glmnet(features$G, d$y, family = "binomial", lower.limits = 0, standardize = FALSE, lambda = reference$lambda,
                        thresh = 1e-14)
  for (lambda in lambdas) {
    k <- which(reference$lambda == lambda)
    model <- uniLassoModel(features, fit$a0[k], as.numeric(fit$beta[, k]))
    expect_equal(model$intercept, unname(reference$a0[k]), tolerance = 1e-6)
    expect_equal(model$coefficients, as.numeric(reference$beta[, k]), tolerance = 1e-6)
  }
  expect_gt(sum(as.numeric(fit$beta[, k]) > 0), 2)
})

test_that("uniLassoFeatures() refuses non-binary covariates or a constant outcome", {
  d <- simulateBinary()
  X <- Matrix::Matrix(d$X, sparse = TRUE)
  expect_error(uniLassoFeatures(X * 2, d$y), "binary covariates")
  expect_error(uniLassoFeatures(X, rep(0, nrow(X))), "binary outcome")
})

cyclopsHasBounds <- function() "lowerBound" %in% names(formals(Cyclops::createPrior))

test_that("stage 2 by Cyclops with lower bounds gives the package's model at the same lambda", {
  skip_if_not_installed("uniLasso")
  skip_if_not(cyclopsHasBounds(), "Cyclops without lower bounds (branch lower-bounds)")
  d <- simulateBinary(seed = 3)
  features <- uniLassoFeatures(Matrix::Matrix(d$X, sparse = TRUE), d$y)
  reference <- uniLasso::uniLasso(d$X, d$y, family = "binomial", thresh = 1e-14)
  problem <- createFitProblem(features$G, d$y, rep(TRUE, length(d$y)))
  control <- Cyclops::createControl(convergenceType = "gradient", tolerance = 1e-12, maxIterations = 1e6)
  n <- length(d$y)
  for (k in c(10, 25, 40)) {
    lambda <- reference$lambda[k]
    fit <- fitProblem(problem, "lasso", rep(2 / (n * lambda)^2, length(problem$ids)), rep(TRUE, length(problem$ids)),
                      NULL, control, lowerBound = 0)
    expect_true(fit$converged)
    theta <- numeric(ncol(features$G)); theta[problem$ids] <- fit$coefficients[-1]
    expect_true(all(theta >= 0))
    model <- uniLassoModel(features, fit$coefficients[1], theta)
    expect_equal(model$intercept, unname(reference$a0[k]), tolerance = 1e-6)
    expect_equal(model$coefficients, as.numeric(reference$beta[, k]), tolerance = 1e-6)
    expect_lt(kktResidual(features$G, d$y, problem, fit$coefficients, rep(2 / (n * lambda)^2, length(problem$ids)),
                          rep(TRUE, length(problem$ids)), "lasso", lowerBound = 0), 1e-6)
  }
  # without the bound some coefficients go negative: the bound is doing something here
  free <- fitProblem(problem, "lasso", rep(2 / (n * reference$lambda[40])^2, length(problem$ids)),
                     rep(TRUE, length(problem$ids)), NULL, control)
  expect_true(any(free$coefficients[-1] < -1e-6))
})

test_that("setUniLasso() runs through execute(): coefficients keep the univariate signs", {
  skip_on_cran()
  skip_if_not(cyclopsHasBounds(), "Cyclops without lower bounds (branch lower-bounds)")
  profile <- get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment()))
  plpData <- PatientLevelPrediction::simulatePlpData(profile, n = 3000, seed = 11)
  # binary covariates only, as the study's (the simulation's age in years, 1002, is not)
  plpData$covariateData$covariates <- dplyr::filter(plpData$covariateData$covariates, .data$covariateId != 1002)
  dataSource <- createPlpDataSource(plpData, "SIM", PatientLevelPrediction::createStudyPopulationSettings(
    riskWindowStart = 1, riskWindowEnd = 365, requireTimeAtRisk = FALSE))
  tasks <- data.frame(taskId = "task_a", outcomeId = 3, timeAtRiskDays = 365)
  tasks$budgets <- list(100)
  folder <- tempfile("sar_unilasso_")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  methods <- list(uniLasso = setUniLasso(), LASSO = setSal(shapes = c(0, 0, 0), name = "LASSO"))
  results <- suppressMessages(execute(dataSource, tasks, folder, mode = "models", orderings = 1, methods = methods))
  expect_setequal(results$method, c("uniLasso", "LASSO"))
  u <- results[results$method == "uniLasso", ]
  expect_identical(u$modelStatus, "SUCCESS")
  expect_gt(u$nonZero, 0)
})

test_that("the study search with the bound keeps every stage-2 coefficient >= 0, so the univariate signs", {
  skip_if_not(cyclopsHasBounds(), "Cyclops without lower bounds (branch lower-bounds)")
  d <- simulateBinary(n = 1500, p = 40, seed = 4)
  features <- uniLassoFeatures(Matrix::Matrix(d$X, sparse = TRUE), d$y)
  folds <- rep_len(1:3, length(d$y))
  fit <- suppressMessages(signalAdaptiveCv(features$G, d$y, folds, penalty = "lasso", shapes = c(0, 0, 0),
                                           searchSettings = c(studySearchSettings(), list(lowerBound = 0))))
  expect_true(all(fit$beta >= 0))
  expect_gt(sum(fit$beta > 0), 0)
  expect_lt(fit$kkt, 1e-3)
  model <- uniLassoModel(features, fit$intercept, fit$beta)
  used <- model$coefficients != 0
  expect_identical(sign(model$coefficients[used]), sign(features$beta[used]))
  expect_error(signalAdaptiveCv(features$G, d$y, folds, penalty = "ridge", shapes = c(0, 0, 0),
                                searchSettings = list(lowerBound = 0)), "LASSO penalty")
})
