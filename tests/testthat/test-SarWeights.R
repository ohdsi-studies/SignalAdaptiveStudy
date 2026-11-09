# Test 1 (PLAN.md): the weights match the original code's on a small fixture (tolerance 1e-8).
# Fixture: the original implementation's log_weights (control_cv_439.py) on synthetic data.

readFixture <- function() {
  folder <- testthat::test_path("fixtures")
  triplets <- read.csv(file.path(folder, "x.csv"))
  X <- Matrix::sparseMatrix(i = triplets$row, j = triplets$col, x = 1,
                            dims = c(triplets$n[1], triplets$p[1]))
  list(X = X, y = read.csv(file.path(folder, "y.csv"))$y,
       expected = read.csv(file.path(folder, "expected.csv"), colClasses = c(log_weight = "numeric")))
}

test_that("weights match the original code for all 12 shapes", {
  fixture <- readFixture()
  shapes <- unique(fixture$expected[, c("c1", "c2", "c3")])
  expect_equal(nrow(shapes), 12)
  for (i in seq_len(nrow(shapes))) {
    shape <- unlist(shapes[i, ])
    expected <- merge(fixture$expected, shapes[i, ])
    expected <- expected[order(expected$col), ]
    result <- computeSarWeights(fixture$X, fixture$y, shape)
    expect_identical(result$keep, expected$keep == 1)
    kept <- result$keep
    expect_lt(max(abs(result$logWeights[kept] - expected$log_weight[kept])), 1e-8)
    expect_true(all(is.na(result$logWeights[!kept])))
  }
})

test_that("c = (0, 0, 0) gives unit weights on the kept columns", {
  fixture <- readFixture()
  result <- computeSarWeights(fixture$X, fixture$y, c(0, 0, 0))
  expect_equal(unname(result$weights[result$keep]), rep(1, sum(result$keep)))
})

test_that("a dense matrix gives the same weights as a sparse one", {
  fixture <- readFixture()
  sparse <- computeSarWeights(fixture$X, fixture$y, c(1, 1, -1))
  dense <- computeSarWeights(as.matrix(fixture$X), fixture$y, c(1, 1, -1))
  expect_equal(dense, sparse)
})

test_that("non-binary columns are kept when they vary, whatever their mean", {
  set.seed(3)
  n <- 200
  binary <- Matrix::Matrix(stats::rbinom(n * 3, 1, 0.3), n, 3, sparse = TRUE)
  X <- cbind(binary, age = stats::runif(n, 18, 90), constantHalf = 0.5, constantOne = 1, allZero = 0)
  y <- stats::rbinom(n, 1, 0.2)
  keep <- computeSarWeights(X, y, c(1, 1, -1))$keep
  expect_identical(unname(keep), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  # binary columns: the original rule (0 < mean < 1), including a column of all ones and one with a single 1
  B <- Matrix::Matrix(cbind(stats::rbinom(n, 1, 0.5), 1, c(1, rep(0, n - 1)), 0), sparse = TRUE)
  expect_identical(computeSarWeights(B, y, c(0, 0, 0))$keep, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("the weights are the same whatever the values of a sparse column, by the mean-square shortcut or X * X", {
  set.seed(3)
  X <- Matrix::rsparsematrix(800, 120, 0.1)
  X@x <- abs(X@x) * stats::runif(length(X@x), 0.5, 40)
  X <- methods::as(X, "CsparseMatrix")
  y <- stats::rbinom(800, 1, 0.15)
  for (shape in list(c(0, 0, 0), c(1, 1, -1), c(1, 0, -0.5), c(0, 1, -1))) {
    expect_identical(computeSarWeights(X, y, shape), computeSarWeights(methods::as(X, "TsparseMatrix"), y, shape))
  }
})
