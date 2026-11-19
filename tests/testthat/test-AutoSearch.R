# The R port of Cyclops's auto-search on synthetic CV curves. (Its agreement with native Cyclops on real fold scores
# was checked during development.)

test_that("the search finds the maximum of a concave curve in log variance", {
  evaluate <- function(v) -(log(v) - log(0.3))^2 + c(-0.1, 0, 0.1) # three folds, maximum at 0.3
  s <- autoSearch(evaluate, start = 1)
  expect_lt(abs(log(s$variance) - log(0.3)), log(1.15))
})

test_that("a flat CV curve ends the search at a tried value, without error", {
  s <- autoSearch(function(v) c(-100, -101, -99), start = 0.05)
  expect_true(s$variance %in% s$tried)
  expect_lte(length(s$tried), 50)
})

test_that("a monotone curve walks to the edge and stops at the step limit", {
  s <- autoSearch(function(v) rep(log(v), 3), start = 1, maxSteps = 6, refine = FALSE)
  expect_equal(s$steps, 6)
  expect_equal(max(s$tried), 1e4) # steps down first (1, 0.1), then up
})

test_that("failed fold fits are left out as in Cyclops (common folds), and a value without any stops the search", {
  scores <- list(c(-1, -2, -3), c(-1, NaN, -3), c(-2, -2, -2))
  e <- foldEstimates(scores)
  expect_equal(e$m, c(-2, -2, -2)) # fold 2 failed somewhere, so it is left out everywhere
  evaluate <- function(v) if (v < 0.01) rep(NaN, 3) else -(log(v) - log(0.001))^2 + c(0, 0.1, 0.2)
  s <- autoSearch(evaluate, start = 0.1)
  expect_true(s$variance >= 0.01) # stopped at the best value tried where folds converged
  expect_error(autoSearch(function(v) rep(NaN, 3), start = 1), "no fold fit converged")
})

test_that("a flat region with rounding noise is a tie: the search moves up, as Cyclops does", {
  # all-zero LASSO models below variance 0.01: equal scores up to the last bit (lower variance 1 ulp "better")
  evaluate <- function(v) {
    if (v < 0.01) return(rep(-115.82961810307479, 3) + (v < 0.001) * 1.4e-14)
    -(log(v) - log(0.3))^2 + c(-0.1, 0, 0.1) - 100
  }
  s <- autoSearch(evaluate, start = 0.0032)
  expect_lt(abs(log(s$variance) - log(0.3)), log(1.5))
})

test_that("the default starting variance does not overflow on large designs", {
  X <- Matrix::sparseMatrix(i = c(1L, 282852L), j = c(1L, 9184L), x = 1, dims = c(282852L, 9184L))
  expect_equal(defaultStartVariance(X, ncol(X)), 2 / (282852 * 9184))
})

test_that("the null-model bound keeps the LASSO search off the all-zero plateau (Cyclops PR #101)", {
  # CV curve in log variance: flat (all-zero) below v0 = 1e-3, a slight dip just above it, the maximum at 0.1
  v0 <- 1e-3
  curve <- function(v) if (v <= v0) -100 else -100.05 + 3 * exp(-(log(v) - log(0.1))^2 / 2) - 0.05 * (v > 0.2)
  evaluate <- function(v) curve(v) + c(-0.1, 0, 0.1)
  start <- 1.5e-3 # just above v0: its score is below the plateau, so stepping down "improves"
  plain <- autoSearch(evaluate, start = start)
  expect_lte(plain$variance, v0) # the trap: the search ends on the plateau
  bounded <- autoSearch(evaluate, start = start, nullVariance = v0)
  expect_gt(bounded$variance, 0.03)
  expect_true(all(bounded$tried >= v0 * (1 - 1e-12)))
  # a start below v0 moves up to it
  expect_identical(autoSearch(evaluate, start = 1e-5, nullVariance = v0)$tried[1], v0)
  # when the null model really is best, the search looks up to 100 x v0 and stops there
  flat <- function(v) if (v <= v0) c(-100, -100, -100) else c(-100, -100, -100) - 1 - log(v / v0)
  stops <- autoSearch(flat, start = start, nullVariance = v0)
  expect_equal(stops$variance, v0)
  expect_lt(max(stops$tried), 1000 * v0)
})

test_that("the null-model variance is where the weighted LASSO becomes all zero", {
  d <- simulateDesign()
  age <- round(stats::runif(nrow(d$X), 0.2, 1), 2)
  X <- cbind(d$X, age)
  free <- c(rep(FALSE, ncol(d$X)), TRUE)
  problem <- createFitProblem(X, d$y, rep(TRUE, nrow(X)))
  tight <- Cyclops::createControl(fitter = "super", convergenceType = "kkt", tolerance = 1e-10, maxIterations = 100000)
  for (shape in list(c(0, 0, 0), c(1, 1, -1))) for (unpenalized in list(rep(FALSE, ncol(X)), free)) {
    w <- problemWeights(X, d$y, problem, shape, unpenalized)
    v0 <- nullModelVariance(X, d$y, problem, w, tight)
    nonZero <- function(v) {
      fit <- fitProblem(problem, "lasso", priorVariances("lasso", w$weights, w$keep, v, "variance", problem$n), w$keep,
                        NULL, tight)
      sum(fit$coefficients[-1][!(w$weights == 0 & w$keep)] != 0)
    }
    expect_identical(nonZero(v0 * 0.98), 0L)
    expect_gt(nonZero(v0 * 1.05), 0L)
  }
})

test_that("the local fit keeps the bracketed step between the best value's neighbours", {
  # a CV curve that still rises almost linearly in log variance, its maximum just bracketed (best 3, then 9.5 a hair
  # lower), tried as the study's search does: the fit through every value is nearly linear and its maximum far out
  x <- c(0.03, 0.003, 0.09, 0.3, 0.9, 3, 9.5)
  m <- -0.0522 + 1e-5 * log(x)
  m[7] <- m[6] - 1e-8
  s <- rep(0.001, length(x))
  global <- uniModalStep(x, m, s, upStep = sqrt(10))
  local <- uniModalStep(x, m, s, upStep = sqrt(10), localFit = TRUE)
  sorted <- sort(x)
  b <- which(sorted == x[which.max(m)])
  expect_true(local$next_ > sorted[b - 1] && local$next_ < sorted[b + 1])
  expect_gt(global$next_, 1000 * max(x)) # the fit through every value lands far outside
  # fewer than three values, or the best at an end: no fit, the same step either way
  expect_equal(uniModalStep(x[1:2], m[1:2], s[1:2], localFit = TRUE), uniModalStep(x[1:2], m[1:2], s[1:2]))
  expect_equal(uniModalStep(x[2:4], m[2:4], s[2:4], localFit = TRUE), uniModalStep(x[2:4], m[2:4], s[2:4]))
})

test_that("the search comparison adds every study model with the local fit, and only those change", {
  methods <- searchComparisonMethods()
  base <- finalStudyMethods()
  local <- methods[grepl("_local$", names(methods))]
  expect_setequal(names(local), paste0(setdiff(names(base), "PLP_LASSO"), "_local"))
  expect_true(all(vapply(local, function(m) isTRUE(m$param$searchSettings$localFit), logical(1))))
  expect_identical(methods[names(base)], base)
  expect_true(all(vapply(base[names(base) != "PLP_LASSO"], function(m) is.null(m$param$searchSettings$localFit), logical(1))))
  expect_identical(local$SAR_no_c3_local$param$searchSettings$sharedShapes, base$SAR_no_c3$param$searchSettings$sharedShapes)
})
