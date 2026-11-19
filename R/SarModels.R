# SAR and SAL as PatientLevelPrediction models. PLP supplies the data, split, folds and evaluation; the fit function
# runs signalAdaptiveCv() on PLP's training rows and folds (option B: our own cross-validation).

#' Signal-adaptive ridge (SAR) model settings
#'
#' Logistic regression with the penalty \eqn{(\lambda/2) \sum_j w_j \beta_j^2}, where the weights come from
#' [computeSarWeights()] on the training rows of each fit. With `shapes = c(0, 0, 0)` this is ridge regression
#' through the same folds, search and code.
#'
#' The penalty is chosen for each shape by Cyclops's auto-search over the prior variance, on PLP's folds, with the
#' weights recomputed from each training fold's rows; c is chosen by the pooled cross-validated log loss.
#'
#' @param shapes Shapes c = (c1, c2, c3) to choose from: a matrix with one row per shape, or one vector.
#'   Default: the 12 shapes of [sarShapes()], the original grid (the final study uses [finalSarShapes()]).
#' @param interceptOnly Allow the intercept-only model to win the cross-validation, as the original does.
#' @param name Model name shown by PLP.
#' @param fitAlgorithm `"auto"`: Cyclops's coordinate descent that moves to Newton or prox-Newton steps where they pay.
#'   `"ccd"`: coordinate descent. Both stop at the same criteria (KKT 3e-3 in folds, Cyclops's default for the final fit),
#'   used for fold and final fits alike.
#' @param sharedShapes `NULL`, or a matrix of shapes that contains `shapes` (for an ablation: [sarShapes()]): the shape
#'   searches are then SAR's own, walked in this order with each shape's fold fits started from the shapes before it, and
#'   taken from SAR's fit on the same training set when it ran first (execute() runs SAR before its ablations); c is
#'   chosen among `shapes`.
#'
#' @return An object of class `modelSettings`.
#' @export
setSar <- function(shapes = sarShapes(), interceptOnly = TRUE, name = "SAR", fitAlgorithm = "auto", sharedShapes = NULL) {
  signalAdaptiveSettings("ridge", shapes, interceptOnly, name, searchSettings = studySearchSettings(fitAlgorithm, sharedShapes))
}

#' Signal-adaptive LASSO (SAL) model settings
#'
#' Logistic regression with the penalty \eqn{\lambda \sum_j \sqrt{w_j} |\beta_j|}, the same weights as SAR. With
#' `shapes = c(0, 0, 0)` this is LASSO through the same folds, search and code.
#'
#' @inheritParams setSar
#' @return An object of class `modelSettings`.
#' @export
setSal <- function(shapes = sarShapes(), interceptOnly = TRUE, name = "SAL", fitAlgorithm = "auto") {
  signalAdaptiveSettings("lasso", shapes, interceptOnly, name, searchSettings = studySearchSettings(fitAlgorithm))
}

#' Posterior-moment ridge
#'
#' Ridge with each coefficient's prior variance proportional to its posterior second moment under a N(0, 1/kappa) pilot
#' prior and the diagonal quadratic approximation of the likelihood at the intercept-only model: weight
#' (h + kappa)^2 / (h + kappa + s^2), with SAR's score s and information h, centred like SAR's weights. The overall
#' variance is searched as for SAR. A secondary method of the study, added after a pilot on the development tasks.
#'
#' @param kappa The pilot prior precision; 0.5 matches SAR's stabiliser.
#' @inheritParams setSar
#' @return PLP model settings.
#' @export
setPosteriorMomentRidge <- function(kappa = 0.5, interceptOnly = TRUE, name = "PM_ridge", fitAlgorithm = "auto") {
  signalAdaptiveSettings("ridge", c(0, 0, 0), interceptOnly, name, searchSettings = studySearchSettings(fitAlgorithm),
                         weightRule = list(type = "posteriorMoment", kappa = kappa))
}

# The study's variance search (autoStrategy()): Cyclops's auto-search started at the "priorScale" variance (c = 1 ridge
# and elastic net, 3 LASSO; elastic net 1 since 2026-11-18 after the CCAE trace), stepping up
# by sqrt(10) while bracketing, each SAR shape's first fold fits started from the nearest shape done (investigator,
# 2026-11-18, after the CCAE fit traces). sharedShapes: see setSar() (SAR's ablations, investigator 2026-11-19).
studySearchSettings <- function(fitAlgorithm = "auto", sharedShapes = NULL) {
  settings <- list(fitAlgorithm = fitAlgorithm, start = "priorScale", upStep = sqrt(10), warmShapes = TRUE)
  if (!is.null(sharedShapes)) settings$sharedShapes <- normalizeShapes(sharedShapes) # only then: other fingerprints as before
  settings
}

# search, cSelection, searchSettings and unpenalizedCovariateIds are for development comparisons; the study uses the
# defaults (Cyclops's auto-search, c chosen by the observed CV loss, every covariate penalised).
signalAdaptiveSettings <- function(penalty, shapes, interceptOnly, name, search = "auto", cSelection = "observed",
                                   searchSettings = list(), unpenalizedCovariateIds = c(), weightRule = NULL) {
  shapes <- normalizeShapes(shapes)
  weightRule <- checkWeightRule(weightRule)
  if (!is.null(weightRule) && (nrow(shapes) != 1 || any(shapes != 0))) stop("a weight rule takes the one shape c(0, 0, 0)")
  if (any(!is.finite(shapes))) stop("shapes must be finite")
  if (!is.null(unpenalizedCovariateIds) && (!is.numeric(unpenalizedCovariateIds) || anyNA(unpenalizedCovariateIds))) {
    stop("unpenalizedCovariateIds must be numeric covariate IDs")
  }
  if (!is.logical(interceptOnly) || length(interceptOnly) != 1) stop("interceptOnly must be TRUE or FALSE")
  structure(list(
    fitFunction = "SignalAdaptiveStudy::fitSignalAdaptive",
    param = list(shapes = shapes, interceptOnly = interceptOnly, search = search, cSelection = cSelection,
                 searchSettings = searchSettings, unpenalizedCovariateIds = unpenalizedCovariateIds,
                 weightRule = weightRule)[c(rep(TRUE, 6), !is.null(weightRule))], # SAR's settings, so their fingerprints, as before
    settings = list(modelName = name, penalty = penalty, modelType = "binary", saveType = "RtoJson")
  ), class = "modelSettings")
}

#' Fit a SAR or SAL model (called by PatientLevelPrediction)
#'
#' PatientLevelPrediction calls this through `modelSettings$fitFunction`; it is exported for that reason only. It
#' runs the cross-validation on the training rows with PLP's folds and returns a `plpModel`.
#'
#' @param trainData PLP training data, with `folds`.
#' @param modelSettings From [setSar()] or [setSal()].
#' @param hyperparameterSettings Not used: the search is part of the model.
#' @param analysisId,analysisPath From PLP.
#' @param ... Not used.
#' @return A `plpModel`.
#' @keywords internal
#' @export
fitSignalAdaptive <- function(trainData, modelSettings, hyperparameterSettings = NULL, analysisId = NULL,
                              analysisPath = NULL, ...) {
  start <- Sys.time()
  param <- modelSettings$param
  settings <- modelSettings$settings
  # The sparse matrix and the fit problems of a training set, shared by the methods that fit on it in turn (execute()
  # runs a unit's methods on the same trainData): same covariates object, rows and folds.
  cached <- trainingSetCache(trainData)
  sparse <- cached$sparse %||% PatientLevelPrediction::toSparseM(trainData, trainData$labels)
  labels <- sparse$labels
  folds <- trainData$folds$index[match(labels$originalRowId, trainData$folds$rowId)]
  if (anyNA(folds) || any(folds < 1)) stop("every training row needs a fold")
  covariateIds <- sparse$covariateMap$covariateId[match(seq_len(ncol(sparse$dataMatrix)), sparse$covariateMap$columnId)]

  fit <- signalAdaptiveCv(sparse$dataMatrix, as.numeric(labels$outcomeCount > 0), folds,
                          penalty = settings$penalty, shapes = param$shapes, interceptOnly = param$interceptOnly,
                          search = param$search %||% "auto", cSelection = param$cSelection %||% "observed",
                          searchSettings = param$searchSettings %||% list(),
                          unpenalized = covariateIds %in% (param$unpenalizedCovariateIds %||% c()),
                          weightRule = param$weightRule, prepared = cached$prepared)
  rememberTrainingSet(trainData, sparse, fit$prepared)
  fit$prepared <- NULL

  ids <- c("(Intercept)", as.character(covariateIds))
  betas <- c(fit$intercept, fit$beta)
  model <- list(
    coefficients = data.frame(betas = betas, covariateIds = ids, stringsAsFactors = FALSE),
    coefficientsExact = data.frame(betas = sprintf("%.17g", betas), covariateIds = ids, stringsAsFactors = FALSE),
    modelType = "logistic",
    penalty = settings$penalty, search = fit$search, cSelection = fit$cSelection,
    selected = fit$selected,
    weights = data.frame(covariateId = covariateIds, weight = fit$weights)
  )
  plpModel <- list(model = model)
  attr(plpModel, "modelType") <- "binary"
  # Training predictions (PLP's runPlp() evaluates them; the study's runs drop them) from the matrix in hand rather than
  # predictSignalAdaptive()'s join over the covariates in DuckDB.
  prediction <- trainData$labels
  prediction$value <- stats::plogis(fit$intercept + as.numeric(sparse$dataMatrix %*% fit$beta)[match(prediction$rowId, labels$originalRowId)])
  prediction$evaluationType <- "Train"
  cvPrediction <- merge(trainData$labels, data.frame(rowId = labels$originalRowId, value = fit$cvPrediction), by = "rowId")
  cvPrediction$evaluationType <- "CV"
  prediction <- rbind(prediction, cvPrediction[, colnames(prediction)])

  covariateImportance <- as.data.frame(trainData$covariateData$covariateRef)
  covariateImportance$covariateValue <- 0
  i <- match(covariateImportance$covariateId, covariateIds)
  covariateImportance$covariateValue[!is.na(i)] <- fit$beta[i[!is.na(i)]]
  covariateImportance$included <- as.integer(!is.na(i))

  delta <- Sys.time() - start
  result <- list(
    model = model,
    preprocessing = list(
      featureEngineering = attr(trainData$covariateData, "metaData")$featureEngineering,
      tidyCovariates = attr(trainData$covariateData, "metaData")$tidyCovariateDataSettings,
      requiresDenseMatrix = FALSE
    ),
    prediction = prediction,
    modelDesign = PatientLevelPrediction::createModelDesign(
      targetId = attr(trainData, "metaData")$targetId,
      outcomeId = attr(trainData, "metaData")$outcomeId,
      restrictPlpDataSettings = attr(trainData, "metaData")$restrictPlpDataSettings,
      covariateSettings = attr(trainData, "metaData")$covariateSettings,
      populationSettings = attr(trainData, "metaData")$populationSettings,
      featureEngineeringSettings = attr(trainData, "metaData")$featureEngineeringSettings,
      preprocessSettings = attr(trainData$covariateData, "metaData")$preprocessSettings,
      modelSettings = modelSettings,
      splitSettings = attr(trainData, "metaData")$splitSettings,
      sampleSettings = attr(trainData, "metaData")$sampleSettings
    ),
    trainDetails = list(
      analysisId = analysisId,
      analysisSource = "",
      developmentDatabase = attr(trainData, "metaData")$cdmDatabaseName,
      developmentDatabaseSchema = attr(trainData, "metaData")$cdmDatabaseSchema,
      attrition = attr(trainData, "metaData")$attrition,
      trainingTime = paste(as.character(abs(delta)), attr(delta, "units")),
      trainingDate = Sys.Date(),
      modelName = settings$modelName,
      finalModelParameters = c(fit$selected, list(penalty = fit$penalty, search = fit$search, weightRule = fit$weightRule,
                                                  kappa = param$weightRule$kappa %||% NA_real_, kkt = fit$kkt,
                                                  nonZero = fit$nonZero, interceptCvLoss = fit$interceptCvLoss,
                                                  cvSeconds = fit$seconds), fit$convergence, fit$timing),
      hyperParamSearch = fit$records,
      searchTrace = fit$trace
    ),
    covariateImportance = covariateImportance
  )
  class(result) <- "plpModel"
  attr(result, "predictionFunction") <- "SignalAdaptiveStudy::predictSignalAdaptive"
  attr(result, "modelType") <- "binary"
  attr(result, "saveType") <- "RtoJson"
  result
}

#' Predict with a SAR or SAL model (called by PatientLevelPrediction)
#'
#' \eqn{\mathrm{plogis}(\beta_0 + x^\top \beta)}, with the coefficients stored at full precision
#' (`coefficientsExact`): PLP's JSON model files keep only about 4 significant digits of numbers. PLP calls this by
#' name; it is exported for that reason only.
#'
#' @param plpModel A `plpModel` from [fitSignalAdaptive()].
#' @param data PLP data with `covariateData`.
#' @param cohort The rows to predict (`rowId` plus cohort columns).
#' @return `cohort` with the predicted risk in `value`.
#' @keywords internal
#' @export
predictSignalAdaptive <- function(plpModel, data, cohort) {
  co <- plpModel$model$coefficientsExact
  if (is.null(co)) co <- plpModel$model$coefficients
  co <- as.data.frame(co, stringsAsFactors = FALSE)
  betas <- as.numeric(co$betas)
  isIntercept <- co$covariateIds == "(Intercept)"
  intercept <- if (any(isIntercept)) betas[isIntercept] else 0
  coefficients <- data.frame(covariateId = as.numeric(co$covariateIds[!isIntercept]), beta = betas[!isIntercept])
  coefficients <- coefficients[coefficients$beta != 0, , drop = FALSE]
  linear <- numeric(nrow(cohort))
  if (nrow(coefficients) > 0) {
    covariateData <- data$covariateData
    covariateData$sarCoefficients <- coefficients
    covariateData$sarRows <- data.frame(rowId = cohort$rowId)
    on.exit({ covariateData$sarCoefficients <- NULL; covariateData$sarRows <- NULL }, add = TRUE)
    terms <- covariateData$covariates %>%
      dplyr::inner_join(covariateData$sarRows, by = "rowId") %>%
      dplyr::inner_join(covariateData$sarCoefficients, by = "covariateId") %>%
      dplyr::transmute(rowId = .data$rowId, covariateId = .data$covariateId,
                       term = .data$covariateValue * .data$beta) %>%
      dplyr::arrange(.data$rowId, .data$covariateId) %>%
      dplyr::collect()
    # Summed in R in a fixed order, so predictions are reproducible to the last bit.
    sums <- rowsum(terms$term, terms$rowId, reorder = FALSE)
    i <- match(cohort$rowId, as.numeric(rownames(sums)))
    linear[!is.na(i)] <- sums[i[!is.na(i)], 1]
  }
  prediction <- cohort
  prediction$value <- stats::plogis(intercept + linear)
  attr(prediction, "metaData") <- list(modelType = "binary")
  prediction
}

# The last training set's sparse matrix and fit problems (fitSignalAdaptive()): one at a time, held until another
# training set comes, or clearTrainingSetCache() (execute() at the end of each unit, to free the memory).
.trainingSetCache <- new.env()
trainingSetKey <- function(trainData) {
  list(covariateData = trainData$covariateData, rows = trainData$labels$rowId, folds = trainData$folds)
}
trainingSetCache <- function(trainData) {
  if (!is.null(.trainingSetCache$key) && identical(.trainingSetCache$key, trainingSetKey(trainData))) {
    list(sparse = .trainingSetCache$sparse, prepared = .trainingSetCache$prepared)
  } else {
    clearTrainingSetCache()
    list()
  }
}
rememberTrainingSet <- function(trainData, sparse, prepared) {
  .trainingSetCache$key <- trainingSetKey(trainData)
  .trainingSetCache$sparse <- sparse
  .trainingSetCache$prepared <- prepared
}
clearTrainingSetCache <- function() {
  rm(list = ls(.trainingSetCache, all.names = TRUE), envir = .trainingSetCache)
  invisible(NULL)
}
