# Posterior-moment ridge (exploratory): algebra, the link to SAR's c = (0, 1, -1), centring, and leakage.

pmStatistics <- function(info, score) list(n = 100, mu = rep(0.1, length(info)), info = info, score = score)
pmRule <- function(kappa = 0.5) list(type = "posteriorMoment", kappa = kappa)

test_that("the weight is the inverse of the pilot posterior second moment", {
  info <- c(0, 0.01, 0.3, 2, 40, 900); score <- c(0, -0.2, 1.5, -3, 12, 40); kappa <- 0.7
  d <- info + kappa
  secondMoment <- (score / d)^2 + 1 / d
  expect_equal(exp(posteriorMomentLogWeights(pmStatistics(info, score), kappa)), 1 / secondMoment, tolerance = 1e-12)
})

test_that("limiting cases: s = 0 gives h + kappa, h = s = 0 gives kappa, the sign of s does not matter", {
  info <- c(0, 0.5, 3, 50); kappa <- 0.5
  expect_equal(exp(posteriorMomentLogWeights(pmStatistics(info, 0 * info), kappa)), info + kappa)
  expect_equal(exp(posteriorMomentLogWeights(pmStatistics(0, 0), 2)), 2)
  score <- c(0.3, -2, 5, -9)
  expect_identical(posteriorMomentLogWeights(pmStatistics(info, score), kappa),
                   posteriorMomentLogWeights(pmStatistics(info, -score), kappa))
})

test_that("at fixed h the weight falls as |s| grows", {
  u <- posteriorMomentLogWeights(pmStatistics(rep(4, 5), c(0, 0.5, 1, 3, 10)), 0.5)
  expect_true(all(diff(u) < 0))
})

test_that("extreme but finite statistics give finite weights; a bad kappa is refused", {
  u <- posteriorMomentLogWeights(pmStatistics(c(0, 1e-300, 1e12, 1e6), c(1e150, 0, 1e9, 1e-200)), 0.5)
  expect_true(all(is.finite(u)))
  for (k in list(0, -1, Inf, NA_real_, c(1, 2), "a")) expect_error(posteriorMomentLogWeights(pmStatistics(1, 1), k))
  expect_error(checkWeightRule(list(type = "other", kappa = 1)), "unknown weight rule")
})

test_that("kappa = 0.5: the ratio to SAR's c = (0, 1, -1) before centring is (h + 0.5) / (h + 1)", {
  fixture <- simulateDesign(n = 800, p = 30, seed = 4)
  s <- sarStatistics(fixture$X, fixture$y)
  sarLog <- log1p(s$info) - log1p(s$score^2 / (s$info + 0.5))
  pmLog <- posteriorMomentLogWeights(s, 0.5)
  expect_equal(pmLog - sarLog, log((s$info + 0.5) / (s$info + 1)), tolerance = 1e-12)
  # after centring the two differ by that ratio up to one common factor
  a <- weightsFromStatistics(s, c(0, 1, -1)); b <- weightsFromStatistics(s, c(0, 0, 0), pmRule())
  k <- a$keep
  r <- b$logWeights[k] - a$logWeights[k] - log((s$info[k] + 0.5) / (s$info[k] + 1))
  expect_lt(max(r) - min(r), 1e-12)
})

test_that("centred weights have geometric mean one over the kept columns, the same columns as SAR", {
  fixture <- simulateDesign(n = 600, p = 25, seed = 2)
  X <- cbind(fixture$X, 0, 1)
  s <- sarStatistics(X, fixture$y)
  w <- weightsFromStatistics(s, c(0, 0, 0), pmRule())
  expect_identical(w$keep, weightsFromStatistics(s, c(1, 1, -1))$keep)
  expect_equal(mean(w$logWeights[w$keep]), 0, tolerance = 1e-12)
  expect_true(all(is.na(w$logWeights[!w$keep])))
})

test_that("a large kappa approaches ridge's equal weights", {
  fixture <- simulateDesign(n = 600, p = 25, seed = 2)
  s <- sarStatistics(fixture$X, fixture$y)
  w <- weightsFromStatistics(s, c(0, 0, 0), pmRule(1e9))
  expect_lt(max(abs(w$logWeights[w$keep])), 1e-6)
})

test_that("fold weights use only the fold's training rows", {
  fixture <- simulateDesign(n = 900, p = 30, seed = 7)
  rows <- fixture$folds != 1
  problem <- list(rows = rows, ids = seq_len(ncol(fixture$X)))
  before <- problemWeights(fixture$X, fixture$y, problem, c(0, 0, 0), weightRule = pmRule())
  y <- fixture$y; y[!rows] <- 1 - y[!rows] # change every held-out outcome
  after <- problemWeights(fixture$X, y, problem, c(0, 0, 0), weightRule = pmRule())
  expect_identical(after$weights, before$weights)
})

test_that("the CV fit uses the weight rule, leaves the intercept unpenalised and refuses other shapes", {
  fixture <- simulateDesign(n = 1500, p = 40, seed = 1)
  pm <- signalAdaptiveCv(fixture$X, fixture$y, fixture$folds, "ridge", shapes = c(0, 0, 0), weightRule = pmRule())
  sar <- signalAdaptiveCv(fixture$X, fixture$y, fixture$folds, "ridge", shapes = c(0, 1, -1))
  expect_identical(pm$weightRule, "posteriorMoment(kappa=0.5)")
  expect_identical(sar$weightRule, "SAR")
  expected <- weightsFromStatistics(sarStatistics(fixture$X, fixture$y), c(0, 0, 0), pmRule())$weights
  expect_equal(pm$weights, expected)
  expect_false(isTRUE(all.equal(pm$beta, sar$beta)))
  expect_error(signalAdaptiveCv(fixture$X, fixture$y, fixture$folds, "ridge", shapes = sarShapes(), weightRule = pmRule()),
               "one shape")
  expect_error(setPosteriorMomentRidge(kappa = -1), "kappa")
})

test_that("SAR settings without a weight rule are unchanged, so their fingerprints are too", {
  expect_false("weightRule" %in% names(setSar()$param))
  expect_identical(setPosteriorMomentRidge()$param$weightRule, pmRule())
})

