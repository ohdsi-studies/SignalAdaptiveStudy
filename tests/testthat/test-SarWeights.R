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
