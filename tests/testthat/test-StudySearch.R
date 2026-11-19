test_that("the study's search settings and preprocessing are as decided", {
  expect_equal(studySearchSettings("auto"), list(fitAlgorithm = "auto", start = "priorScale", upStep = sqrt(10), warmShapes = TRUE))
  expect_equal(setSar()$param$searchSettings, studySearchSettings("auto"))
  expect_equal(setElasticNet(fitAlgorithm = "ccd")$param$searchSettings, studySearchSettings("ccd"))
  expect_equal(studyPreprocessSettings()$minFraction, 0.001)
  expect_true(studyPreprocessSettings()$removeRedundancy)
})

test_that("priorScaleStart is c n / sum ||x_j||^2 / w_j over the penalised columns kept", {
  squares <- c(10, 20, 40, 5)
  weights <- list(weights = c(1, 2, 0, NA), keep = c(TRUE, TRUE, TRUE, FALSE)) # column 3 unpenalised, 4 not kept
  expect_equal(priorScaleStart(squares, weights, 100, "ridge"), 100 / (10 + 10))
  expect_equal(priorScaleStart(squares, weights, 100, "lasso"), 3 * 100 / 20)
  expect_equal(priorScaleStart(squares, weights, 100, "elasticnet"), 100 / 20)
  # binary columns, equal weights: 1 / (mean non-zeros per row)
  X <- Matrix::rsparsematrix(200, 50, 0.1, rand.x = function(n) rep(1, n))
  expect_equal(priorScaleStart(Matrix::colSums(X^2), list(weights = rep(1, 50), keep = rep(TRUE, 50)), 200, "ridge"),
               1 / mean(Matrix::rowSums(X)))
})

test_that("autoSearch steps up by upStep and down by 10", {
  peak <- log(0.5)
  evaluate <- function(v) -(log(v) - peak)^2 + c(0, 0.01, -0.01) # three folds
  cyclops <- autoSearch(evaluate, start = 0.01)
  study <- autoSearch(evaluate, start = 0.01, upStep = sqrt(10))
  expect_equal(cyclops$tried[1:3], c(0.01, 0.001, 0.1))
  expect_equal(study$tried[1:3], c(0.01, 0.001, 0.01 * sqrt(10)))
  expect_equal(study$variance, 0.5, tolerance = 0.15)
  expect_equal(cyclops$variance, 0.5, tolerance = 0.15)
})

test_that("shapeSeeds take the nearest shape's fold fits at its selected variance", {
  state <- function(k) lapply(1:2, function(f) list(variances = c(0.1, 1), coefficients = list(c(k, f, 0.1), c(k, f, 1))))
  done <- list(list(shape = c(0, 0, 0), state = state(1), variance = 1), list(shape = c(1, 0, -1), state = state(2), variance = 0.1))
  seeds <- shapeSeeds(c(1, 0, -0.5), done, 2)
  expect_equal(seeds, list(c(2, 1, 0.1), c(2, 2, 0.1)))
  expect_equal(shapeSeeds(c(0, 1, 0), done, 2), list(c(1, 1, 1), c(1, 2, 1)))
  expect_equal(shapeSeeds(c(0, 0, 0), list(), 2), list(NULL, NULL))
})

test_that("the study's search selects as Cyclops's start and steps do, within the refinement's resolution", {
  set.seed(5)
  n <- 3000; p <- 60
  X <- Matrix::rsparsematrix(n, p, 0.15, rand.x = function(k) rep(1, k))
  y <- stats::rbinom(n, 1, stats::plogis(-2 + as.numeric(X[, 1:6] %*% rep(0.8, 6))))
  folds <- rep(1:3, length.out = n)
  shapes <- rbind(c(0, 0, 0), c(0, 1, -1), c(1, 0, -1))
  old <- signalAdaptiveCv(X, y, folds, "ridge", shapes = shapes, searchSettings = list(fitAlgorithm = "ccd"))
  new <- signalAdaptiveCv(X, y, folds, "ridge", shapes = shapes, searchSettings = studySearchSettings("ccd"))
  expect_equal(new$selected$shape, old$selected$shape)
  expect_lt(abs(log(new$selected$variance) - log(old$selected$variance)), log(1.5))
  expect_equal(new$selected$cvLoss, old$selected$cvLoss, tolerance = 1e-3)
})

test_that("shape searches are shared on a training set: ridge and LASSO as alone, the ablations as SAR's", {
  set.seed(7)
  n <- 2000; p <- 50
  X <- Matrix::rsparsematrix(n, p, 0.15, rand.x = function(k) rep(1, k))
  y <- stats::rbinom(n, 1, stats::plogis(-2 + as.numeric(X[, 1:6] %*% rep(0.8, 6))))
  folds <- rep(1:3, length.out = n)
  all <- sarShapes()[c(1, 2, 3, 7, 9), ] # (0,0,0) first, as in sarShapes()
  fit <- function(penalty, shapes, settings = studySearchSettings("ccd"), prepared = NULL) {
    signalAdaptiveCv(X, y, folds, penalty, shapes = shapes, searchSettings = settings, prepared = prepared)
  }
  same <- function(a, b) {
    expect_identical(a$beta, b$beta); expect_identical(a$intercept, b$intercept)
    expect_identical(a$selected$shape, b$selected$shape); expect_identical(a$selected$variance, b$selected$variance)
  }
  sar <- fit("ridge", all)
  expect_equal(sar$reusedShapes, 0L)
  # ridge after SAR: SAR's first shape, started from no other: the same as ridge alone
  ridgeAlone <- fit("ridge", c(0, 0, 0))
  ridgeAfter <- fit("ridge", c(0, 0, 0), prepared = sar$prepared)
  expect_equal(ridgeAfter$reusedShapes, 1L)
  same(ridgeAfter, ridgeAlone)
  # an ablation after SAR takes all of SAR's searches; without SAR first it walks SAR's shapes itself: the same model
  sub <- all[all[, "c3"] == 0, , drop = FALSE]
  shared <- studySearchSettings("ccd", sharedShapes = all)
  after <- fit("ridge", sub, shared, prepared = sar$prepared)
  expect_equal(after$reusedShapes, nrow(all))
  same(after, fit("ridge", sub, shared))
  expect_equal(after$records, sar$records[sar$records$c3 == 0, ], ignore_attr = TRUE)
  # a search started from no other shape is not taken by a method that would start that shape from others
  fixedAlone <- fit("ridge", c(0, 1, -1))
  sarAfterFixed <- fit("ridge", all, prepared = fixedAlone$prepared)
  same(sarAfterFixed, sar)
  # the penalty is part of what is shared: LASSO does not take ridge's search
  lassoAfter <- fit("lasso", c(0, 0, 0), prepared = sar$prepared)
  expect_equal(lassoAfter$reusedShapes, 0L)
  same(lassoAfter, fit("lasso", c(0, 0, 0)))
  expect_error(fit("ridge", rbind(c(2, 0, 0)), shared), "must contain every shape")
  # other methods' settings, so their fingerprints, unchanged
  expect_null(setSar()$param$searchSettings$sharedShapes)
  expect_equal(finalStudyMethods()$SAR_no_c3$param$searchSettings$sharedShapes, finalSarShapes(), ignore_attr = TRUE)
})

test_that("with sarShareSearches FALSE every method searches its shapes itself, to the same model", {
  set.seed(7)
  n <- 1500; p <- 40
  X <- Matrix::rsparsematrix(n, p, 0.15, rand.x = function(k) rep(1, k))
  y <- stats::rbinom(n, 1, stats::plogis(-2 + as.numeric(X[, 1:5] %*% rep(0.8, 5))))
  folds <- rep(1:3, length.out = n)
  settings <- studySearchSettings("ccd")
  sar <- signalAdaptiveCv(X, y, folds, "ridge", shapes = sarShapes()[1:3, ], searchSettings = settings)
  old <- options(sarShareSearches = FALSE)
  on.exit(options(old))
  ridge <- signalAdaptiveCv(X, y, folds, "ridge", shapes = c(0, 0, 0), searchSettings = settings, prepared = sar$prepared)
  expect_equal(ridge$reusedShapes, 0L)
  expect_identical(ridge$beta, signalAdaptiveCv(X, y, folds, "ridge", shapes = c(0, 0, 0), searchSettings = settings)$beta)
})
