# The study split: PLP's default 75/25 test set, seeded orderings cut at the k-th event, 3 stratified R folds.

test_that("the budget prefix ends at the k-th event", {
  y <- c(0, 0, 1, 0, 1, 1, 0, 0, 1)
  expect_identical(budgetPrefixLength(y, 1), 3L)
  expect_identical(budgetPrefixLength(y, 3), 6L)
  expect_identical(budgetPrefixLength(y, 4), 9L)
  expect_error(budgetPrefixLength(y, 5), "infeasible")
})

test_that("R folds are reproducible and leave the global random stream alone", {
  set.seed(1); before <- .Random.seed
  y <- rep(c(0, 1), c(90, 10))
  f1 <- stratifiedFolds(y, 5, 1729)
  expect_identical(.Random.seed, before)
  expect_identical(stratifiedFolds(y, 5, 1729), f1)
  expect_identical(as.vector(table(f1[y == 1])), rep(2L, 5))
})

test_that("study split: PLP test set fixed across budgets and orderings, prefix to the k-th event", {
  population <- simulatePopulation()
  set.seed(1); before <- .Random.seed
  splits <- list()
  for (ordering in 1:3) for (k in c(25, 100)) splits[[paste(ordering, k)]] <- studySplit(population, list(), k, ordering)
  expect_identical(.Random.seed, before)
  test <- expectSplitProperties(splits, population, c(25, 100), nfold = 3)
  expect_equal(length(test) / nrow(population), 0.25, tolerance = 0.01)
  # the orderings differ
  first <- function(s) s$rowId[s$index > 0]
  expect_false(identical(first(splits[["1 25"]]), first(splits[["2 25"]])))
  # the whole pool when there is no budget
  all <- studySplit(population, list(), NULL, 1)
  expect_identical(attr(all, "trainingN"), nrow(population) - length(test))
})
