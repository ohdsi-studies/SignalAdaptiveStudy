# Tests 2-4 of PLAN.md.

# Test 2: with c = (0, 0, 0), SAR is ridge and SAL is LASSO (glmnet as the independent reference, same lambda).
test_that("SAR and SAL with c = (0, 0, 0) equal ridge and LASSO at the same lambda", {
  skip_if_not_installed("glmnet")
  d <- simulateDesign()
  tight <- Cyclops::createControl(fitter = "super", convergenceType = "kkt", tolerance = 1e-9, maxIterations = 100000)
  for (penalty in c("ridge", "lasso")) {
    # one fit at a fixed lambda (the lambda <-> prior variance conversion)
    problem <- createFitProblem(d$X, d$y, rep(TRUE, nrow(d$X)))
    w <- problemWeights(d$X, d$y, problem, c(0, 0, 0))
    lambda <- if (penalty == "ridge") 0.01 else 0.005
    fit <- fitProblem(problem, penalty, priorVariances(penalty, w$weights, w$keep, lambda, "lambda", problem$n), w$keep,
                      NULL, tight)
    reference <- glmnetFit(d$X, d$y, lambda, penalty)
    beta <- numeric(ncol(d$X)); beta[problem$ids] <- fit$coefficients[-1]
    expect_lt(abs(fit$coefficients[1] - reference$intercept), 1e-6)
    expect_lt(max(abs(beta - reference$beta)), 1e-6)

    # the cross-validated pipeline: the final model is ridge / LASSO at the selected lambda
    cv <- signalAdaptiveCv(d$X, d$y, d$folds, penalty = penalty, shapes = c(0, 0, 0))
    expect_false(cv$selected$interceptOnly)
    reference <- glmnetFit(d$X, d$y, cv$selected$lambda, penalty)
    expect_lt(abs(cv$intercept - reference$intercept), 1e-3) # Cyclops "super" final-fit tolerance (kkt 1e-4)
    expect_lt(max(abs(cv$beta - reference$beta)), 1e-3)
  }
})

# Test 3: changing the held-out fold's labels does not change the weights (or the fit) of the training fold.
test_that("held-out labels do not change the training fold's weights or fit", {
  d <- simulateDesign(seed = 2)
  heldOut <- d$folds == 1
  yChanged <- d$y
  yChanged[heldOut] <- 1 - yChanged[heldOut]
  problem <- createFitProblem(d$X, d$y, !heldOut)
  problemChanged <- createFitProblem(d$X, yChanged, !heldOut)
  for (shape in list(c(1, 1, -1), c(0, 1, -0.5), c(1, 0, 0))) {
    w <- problemWeights(d$X, d$y, problem, shape)
    wChanged <- problemWeights(d$X, yChanged, problemChanged, shape)
    expect_identical(wChanged$weights, w$weights)
    expect_identical(wChanged$keep, w$keep)
    for (penalty in c("ridge", "lasso")) {
      v <- priorVariances(penalty, w$weights, w$keep, 0.01, "lambda", problem$n)
      control <- Cyclops::createControl(fitter = "super", convergenceType = "kkt", tolerance = 1e-6)
      fit <- fitProblem(problem, penalty, v, w$keep, NULL, control)
      fitChanged <- fitProblem(problemChanged, penalty, v, wChanged$keep, NULL, control)
      expect_identical(fitChanged$coefficients, fit$coefficients)
    }
  }
  # and the weights do depend on the training rows' labels (so the test can fail)
  yTrainChanged <- d$y
  yTrainChanged[!heldOut][1:50] <- 1 - yTrainChanged[!heldOut][1:50]
  expect_false(identical(problemWeights(d$X, yTrainChanged, problem, c(1, 1, -1))$weights,
                         problemWeights(d$X, d$y, problem, c(1, 1, -1))$weights))
})

# Test 4: save -> load -> predict gives identical predictions.
test_that("a saved and reloaded model predicts identically", {
  plpData <- PatientLevelPrediction::simulatePlpData(
    get(utils::data("simulationProfile", package = "PatientLevelPrediction", envir = environment())),
    n = 1500, seed = 4
  )
  saveDirectory <- tempfile()
  on.exit(unlink(saveDirectory, recursive = TRUE), add = TRUE)
  result <- PatientLevelPrediction::runPlp(
    plpData = plpData, outcomeId = 3, analysisId = "saveLoad",
    populationSettings = PatientLevelPrediction::createStudyPopulationSettings(riskWindowStart = 1, riskWindowEnd = 365,
                                                                               requireTimeAtRisk = FALSE),
    splitSettings = PatientLevelPrediction::createDefaultSplitSetting(splitSeed = 1, nfold = 5),
    preprocessSettings = PatientLevelPrediction::createPreprocessSettings(minFraction = 0, normalize = TRUE,
                                                                          removeRedundancy = FALSE),
    modelSettings = setSal(shapes = rbind(c(0, 0, 0), c(1, 1, -1))),
    logSettings = PatientLevelPrediction::createLogSettings(verbosity = "FATAL"),
    executeSettings = PatientLevelPrediction::createExecuteSettings(
      runSplitData = TRUE, runSampleData = FALSE, runFeatureEngineering = FALSE, runPreprocessData = TRUE,
      runModelDevelopment = TRUE, runCovariateSummary = FALSE),
    saveDirectory = saveDirectory
  )
  modelDirectory <- file.path(saveDirectory, "model")
  PatientLevelPrediction::savePlpModel(result$model, modelDirectory)
  reloaded <- PatientLevelPrediction::loadPlpModel(modelDirectory)
  expect_identical(attr(reloaded, "predictionFunction"), "SignalAdaptiveStudy::predictSignalAdaptive")

  population <- result$prediction[result$prediction$evaluationType == "Test", c("rowId", "outcomeCount")]
  before <- PatientLevelPrediction::predictPlp(result$model, plpData, population)
  after <- PatientLevelPrediction::predictPlp(reloaded, plpData, population)
  expect_identical(after$value, before$value)
  expect_gt(stats::sd(before$value), 0)
})

test_that("unpenalised covariates get a flat prior and stay out of the weights' centring", {
  d <- simulateDesign()
  age <- round(stats::runif(nrow(d$X), 0.2, 1), 2) # a continuous column, as PLP's normalised age
  X <- cbind(d$X, age)
  y <- stats::rbinom(nrow(X), 1, stats::plogis(-4 + 3 * age + as.numeric(d$X[, 1:5] %*% rep(0.7, 5))))
  free <- c(rep(FALSE, ncol(d$X)), TRUE)
  problem <- createFitProblem(X, y, rep(TRUE, nrow(X)))
  shape <- c(1, 1, -1)
  w <- problemWeights(X, y, problem, shape, free)
  reference <- computeSarWeights(d$X, y, shape)
  expect_equal(w$full$weights[!free], reference$weights)
  expect_identical(w$full$weights[free], 0)
  expect_true(w$full$keep[free])

  # at a fixed lambda: the gradient of the log likelihood is 0 for the unpenalised column, lambda * w * beta otherwise
  tight <- Cyclops::createControl(fitter = "super", convergenceType = "kkt", tolerance = 1e-10, maxIterations = 100000)
  lambda <- 0.05
  v <- priorVariances("ridge", w$weights, w$keep, lambda, "lambda", problem$n)
  expect_true(is.infinite(v[problem$ids == ncol(X)]))
  fit <- fitProblem(problem, "ridge", v, w$keep, NULL, tight)
  beta <- numeric(ncol(X)); beta[problem$ids] <- fit$coefficients[-1]
  gradient <- as.numeric(Matrix::crossprod(X, y - stats::plogis(fit$coefficients[1] + as.numeric(X %*% beta)))) / nrow(X)
  expect_lt(abs(gradient[free]), 1e-8)
  kept <- which(w$full$keep & !free)
  expect_lt(max(abs(gradient[kept] - lambda * w$full$weights[kept] * beta[kept])), 1e-8)

  # through the cross-validated pipeline: penalised with the others, age is shrunk; unpenalised, it is not
  for (penalty in c("ridge", "lasso")) {
    shrunk <- signalAdaptiveCv(X, y, d$folds, penalty = penalty, shapes = c(0, 0, 0))
    unshrunk <- signalAdaptiveCv(X, y, d$folds, penalty = penalty, shapes = c(0, 0, 0), unpenalized = free)
    expect_gt(unshrunk$beta[ncol(X)], shrunk$beta[ncol(X)])
    expect_identical(unshrunk$weights[ncol(X)], 0)
  }
  expect_error(signalAdaptiveCv(X, y, d$folds, unpenalized = TRUE), "one TRUE/FALSE per column")
})
