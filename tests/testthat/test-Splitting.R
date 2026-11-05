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

test_that("repeated entries of a person stay on one side of the split, in the ordering and in the folds", {
  set.seed(1)
  n <- 3000
  population <- data.frame(rowId = seq_len(n), subjectId = c(seq_len(n - 300), sample(n - 300, 300)), cohortStartDate = as.Date("2020-01-01"),
                           outcomeCount = rbinom(n, 1, 0.1))
  s <- studySplit(population, list(taskId = "t"), 100, 1)
  subjects <- population$subjectId[match(s$rowId, population$rowId)]
  expect_length(intersect(subjects[s$index < 0], subjects[s$index > 0]), 0)
  training <- s$index > 0
  expect_true(all(tapply(s$index[training], subjects[training], function(x) length(unique(x))) == 1))
  expect_equal(attr(s, "trainingEvents"), sum(population$outcomeCount[match(s$rowId[training], population$rowId)]))
  expect_gte(attr(s, "trainingEvents"), 100)
  again <- studySplit(population, list(taskId = "t"), 100, 1)
  expect_identical(s$rowId, again$rowId)
  expect_identical(s$index, again$index)
  # smaller budgets stay prefixes of larger ones
  small <- studySplit(population, list(taskId = "t"), 50, 1)
  expect_true(all(small$rowId[small$index > 0] %in% s$rowId[s$index > 0]))
})

test_that("a population with one entry per person keeps its split", {
  set.seed(2)
  n <- 3000
  population <- data.frame(rowId = seq_len(n), subjectId = seq_len(n), cohortStartDate = as.Date("2020-01-01"), outcomeCount = rbinom(n, 1, 0.1))
  noSubject <- population[, c("rowId", "cohortStartDate", "outcomeCount")]
  expect_false(hasRepeatedPersons(population))
  expect_identical(studySplit(population, list(taskId = "t"), 100, 2)$rowId, studySplit(noSubject, list(taskId = "t"), 100, 2)$rowId)
})
