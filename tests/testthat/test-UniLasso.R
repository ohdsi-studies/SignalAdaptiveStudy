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
