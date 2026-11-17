# Equal-weight models (ridge, LASSO, elastic net) use one Cyclops prior with the intercept excluded; the same fits as one
# prior per column.

perColumnFit <- function(problem, penalty, variances, keep, control) {
  type <- switch(penalty, ridge = "normal", lasso = "laplace", elasticnet = "elasticnet")
  prior <- Cyclops::createPrior(c("none", rep(type, length(variances))), variance = c(0, variances),
                                laplaceFraction = if (penalty == "elasticnet") ELASTIC_NET_ALPHA)
  fixed <- if (any(!keep)) c(FALSE, !keep)
  fit <- Cyclops::fitCyclopsModel(problem$data, prior = prior, control = control, fixedCoefficients = fixed, warnings = FALSE,
                                  startingCoefficients = if (any(!keep)) numeric(length(variances) + 1))
  as.numeric(stats::coef(fit, ignoreConvergence = TRUE))
}

test_that("one global prior gives the same coefficients as one prior per column, for every penalty and algorithm", {
  d <- simulateDesign(n = 2000, p = 60, seed = 3)
  X <- cbind(d$X, 0) # a column that is constant in the rows: not kept, fixed at 0
  problem <- createFitProblem(X, d$y, rep(TRUE, length(d$y)))
  keep <- c(rep(TRUE, length(problem$ids) - 1), FALSE)
  for (algorithm in c("ccd", "auto")) {
    control <- if (algorithm == "ccd") {
      Cyclops::createControl(fitter = "super", algorithm = "ccd", convergenceType = "kkt", tolerance = 1e-10, maxIterations = 1e5, threads = 1)
    } else {
      Cyclops::createControl(fitter = "super", algorithm = "auto", convergenceType = "gap", tolerance = 1e-10, maxIterations = 1e5, threads = 1)
    }
    for (penalty in c("ridge", "lasso", "elasticnet")) {
      v <- rep(0.05, length(problem$ids))
      global <- fitProblem(problem, penalty, v, keep, NULL, control)
      expect_true(global$converged, info = paste(algorithm, penalty))
      expect_equal(global$coefficients, perColumnFit(problem, penalty, v, keep, control), tolerance = 1e-8,
                   info = paste(algorithm, penalty))
      expect_equal(global$coefficients[length(global$coefficients)], 0)
    }
  }
})

test_that("unequal variances or an unpenalised column keep one prior per column", {
  d <- simulateDesign(n = 1500, p = 30, seed = 4)
  problem <- createFitProblem(d$X, d$y, rep(TRUE, length(d$y)))
  keep <- rep(TRUE, length(problem$ids))
  control <- Cyclops::createControl(fitter = "super", algorithm = "ccd", convergenceType = "kkt", tolerance = 1e-10, maxIterations = 1e5, threads = 1)
  v <- seq(0.01, 0.1, length.out = length(problem$ids))
  expect_equal(fitProblem(problem, "ridge", v, keep, NULL, control)$coefficients, perColumnFit(problem, "ridge", v, keep, control),
               tolerance = 1e-8)
  v <- rep(0.05, length(problem$ids)); v[1] <- Inf
  fit <- fitProblem(problem, "lasso", v, keep, NULL, control)
  expect_true(fit$converged)
  expect_true(fit$coefficients[2] != 0) # unpenalised
})
