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
#'   Default: the 12 shapes of [sarShapes()].
#' @param interceptOnly Allow the intercept-only model to win the cross-validation, as the original does.
#' @param name Model name shown by PLP.
#'
#' @return An object of class `modelSettings`.
#' @export
setSar <- function(shapes = sarShapes(), interceptOnly = TRUE, name = "SAR") {
  signalAdaptiveSettings("ridge", shapes, interceptOnly, name)
}

#' Signal-adaptive LASSO (SAL) model settings
#'
#' Logistic regression with the penalty \eqn{\lambda \sum_j \sqrt{w_j} |\beta_j|}, the same weights as SAR. With
#' `shapes = c(0, 0, 0)` this is LASSO through the same folds, search and code.
#'
#' @inheritParams setSar
#' @return An object of class `modelSettings`.
#' @export
setSal <- function(shapes = sarShapes(), interceptOnly = TRUE, name = "SAL") {
  signalAdaptiveSettings("lasso", shapes, interceptOnly, name)
}

# search, cSelection, searchSettings and unpenalizedCovariateIds are for development comparisons; the study uses the
# defaults (Cyclops's auto-search, c chosen by the observed CV loss, every covariate penalised).
signalAdaptiveSettings <- function(penalty, shapes, interceptOnly, name, search = "auto", cSelection = "observed",
                                   searchSettings = list(), unpenalizedCovariateIds = c()) {
  shapes <- normalizeShapes(shapes)
  if (any(!is.finite(shapes))) stop("shapes must be finite")
  if (!is.null(unpenalizedCovariateIds) && (!is.numeric(unpenalizedCovariateIds) || anyNA(unpenalizedCovariateIds))) {
    stop("unpenalizedCovariateIds must be numeric covariate IDs")
  }
  if (!is.logical(interceptOnly) || length(interceptOnly) != 1) stop("interceptOnly must be TRUE or FALSE")
  structure(list(
    fitFunction = "SignalAdaptiveStudy::fitSignalAdaptive",
    param = list(shapes = shapes, interceptOnly = interceptOnly, search = search, cSelection = cSelection,
                 searchSettings = searchSettings, unpenalizedCovariateIds = unpenalizedCovariateIds),
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
  sparse <- PatientLevelPrediction::toSparseM(trainData, trainData$labels)
  labels <- sparse$labels
  folds <- trainData$folds$index[match(labels$originalRowId, trainData$folds$rowId)]
  if (anyNA(folds) || any(folds < 1)) stop("every training row needs a fold")
  covariateIds <- sparse$covariateMap$covariateId[match(seq_len(ncol(sparse$dataMatrix)), sparse$covariateMap$columnId)]

  fit <- signalAdaptiveCv(sparse$dataMatrix, as.numeric(labels$outcomeCount > 0), folds,
                          penalty = settings$penalty, shapes = param$shapes, interceptOnly = param$interceptOnly,
                          search = param$search %||% "auto", cSelection = param$cSelection %||% "observed",
                          searchSettings = param$searchSettings %||% list(),
                          unpenalized = covariateIds %in% (param$unpenalizedCovariateIds %||% c()))

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
  prediction <- predictSignalAdaptive(plpModel, trainData, trainData$labels)
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
      finalModelParameters = c(fit$selected, list(penalty = fit$penalty, search = fit$search, kkt = fit$kkt,
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
