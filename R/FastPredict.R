# Test-set predictions of the study's own models without PLP's predictPlp(), which for every model copies the test
# covariates through the model's preprocessing in DuckDB, then joins, sorts and collects every term. The test set is the
# same for every model of a task (the study split fixes it), so its covariates are read once per task in each worker
# (testDesign()) and each model's predictions are computed from them in R, with the same arithmetic in the same order
# as predictPlp() with predictSignalAdaptive(): the value normalised as PLP does, (1 * value) / maxValue, times the
# coefficient, summed per person in covariate ID order from 0, terms with a zero coefficient left out. So the
# predictions are identical to the last bit (tests/testthat/test-FastPredict.R).

# The test covariates of the persons in test$labels, by person in covariate ID order: covariate indices into ids, the
# number of entries per person, and the values (NULL when all are 1, as for the study's binary covariates).
testDesign <- function(test) {
  covariateData <- test$covariateData
  covariateData$sarTestRows <- data.frame(rowId = test$labels$rowId)
  on.exit(covariateData$sarTestRows <- NULL, add = TRUE)
  entries <- covariateData$covariates %>%
    dplyr::inner_join(covariateData$sarTestRows, by = "rowId") %>%
    dplyr::select("rowId", "covariateId", "covariateValue") %>%
    dplyr::arrange(.data$rowId, .data$covariateId) %>%
    dplyr::collect()
  ids <- sort(unique(entries$covariateId))
  runs <- rle(entries$rowId)
  list(rowIds = test$labels$rowId, persons = runs$values, lengths = runs$lengths, ids = ids,
       covariate = match(entries$covariateId, ids),
       value = if (all(entries$covariateValue == 1)) NULL else as.numeric(entries$covariateValue))
}

# The design of the test set in data, from the worker's cache when the test persons are the same.
cachedTestDesign <- function(cache, test) {
  if (is.null(cache$testDesign) || !identical(cache$testDesign$rowIds, test$labels$rowId)) {
    cache$testDesign <- NULL
    cache$testDesign <- testDesign(test)
  }
  cache$testDesign
}

# The linear predictor (without the intercept) of each person of the design, for coefficients (covariateId, beta) and
# the model's tidyCovariates preprocessing (PLP's applyTidyCovariateData(): covariates removed in training dropped,
# values divided by the training maximum).
designLinearPredictor <- function(design, coefficients, tidyCovariates, chunk = 50000) {
  beta <- numeric(length(design$ids))
  at <- match(coefficients$covariateId, design$ids)
  beta[at[!is.na(at)]] <- coefficients$beta[!is.na(at)]
  deleted <- c(tidyCovariates$deletedRedundantCovariateIds, tidyCovariates$deletedInfrequentCovariateIds)
  beta[design$ids %in% deleted] <- 0
  maxs <- tidyCovariates$normFactors
  maxValue <- NULL
  if (!is.null(maxs)) {
    if ("bins" %in% colnames(maxs)) maxs <- data.frame(covariateId = maxs$bins, maxValue = maxs$maxs)
    maxValue <- maxs$maxValue[match(design$ids, maxs$covariateId)]
    beta[is.na(maxValue)] <- 0 # PLP's inner join drops covariates without a normalisation factor
  }
  linear <- numeric(length(design$persons))
  ends <- cumsum(as.numeric(design$lengths)); starts <- ends - design$lengths
  for (first in seq(1, length(design$persons), by = chunk)) {
    people <- first:min(first + chunk - 1, length(design$persons))
    range <- if (ends[people[length(people)]] > starts[people[1]]) (starts[people[1]] + 1):ends[people[length(people)]] else integer()
    j <- design$covariate[range]
    use <- beta[j] != 0
    if (!any(use)) next
    value <- if (is.null(design$value)) rep(1, length(range)) else design$value[range]
    normalised <- if (is.null(maxValue)) value[use] else 1 * value[use] / maxValue[j[use]]
    person <- rep(people, design$lengths[people])[use]
    sums <- rowsum(normalised * beta[j[use]], person, reorder = FALSE)
    linear[as.integer(rownames(sums))] <- sums[, 1]
  }
  linear
}

# predictPlp(model, data$Test, data$Test$labels) for the study's own models, from the design of data$Test.
predictFromDesign <- function(plpModel, design, labels) {
  co <- plpModel$model$coefficientsExact %||% plpModel$model$coefficients
  co <- as.data.frame(co, stringsAsFactors = FALSE)
  betas <- as.numeric(co$betas)
  isIntercept <- co$covariateIds == "(Intercept)"
  intercept <- if (any(isIntercept)) betas[isIntercept] else 0
  coefficients <- data.frame(covariateId = as.numeric(co$covariateIds[!isIntercept]), beta = betas[!isIntercept])
  coefficients <- coefficients[coefficients$beta != 0, , drop = FALSE]
  linear <- designLinearPredictor(design, coefficients, plpModel$preprocessing$tidyCovariates)
  i <- match(labels$rowId, design$persons)
  prediction <- labels
  prediction$value <- stats::plogis(intercept + ifelse(is.na(i), 0, linear[i]))
  attr(prediction, "metaData") <- list(modelType = "binary")
  prediction
}

isStudyModel <- function(modelSettings) identical(modelSettings$fitFunction, "SignalAdaptiveStudy::fitSignalAdaptive")
