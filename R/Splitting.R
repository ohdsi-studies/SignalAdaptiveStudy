# Training subsets for the study: PLP's default 75/25 split fixes the test set for every budget, ordering and
# method; each ordering is a seeded random order of the training pool, cut at the k-th event (so smaller budgets are
# prefixes of larger ones); 3 folds (PLP's default) are stratified within each training subset. 3 against 5 folds
# was checked on the replication tasks.

STUDY_SEEDS <- list(orderings = c(2026101603, 2026101604, 2026101605), folds = 1729, testSplit = 1729)

# Runs expr with its own random seed and leaves the global random stream as it was.
withSeed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv())) {
    saved <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

# The number of leading rows that end at the k-th event.
budgetPrefixLength <- function(orderedOutcome, events) {
  hits <- which(orderedOutcome > 0)
  if (length(hits) < events) stop("infeasible: fewer than ", events, " events in the training pool")
  hits[events]
}

# Stratified folds: within each outcome class, rows in random order dealt to the folds in turn.
stratifiedFolds <- function(y, nSplits, seed) {
  withSeed(seed, {
    folds <- integer(length(y))
    for (cls in sort(unique(y))) {
      i <- which(y == cls)
      folds[i] <- sample(rep_len(seq_len(nSplits), length(i)))
    }
    folds
  })
}

# PLP splitIds (test -1, training rows their fold; training-pool rows beyond the budget left out) and counts.
assembleSplit <- function(population, training, trainingFolds, isTest) {
  outcome <- as.integer(population$outcomeCount > 0)
  structure(data.frame(rowId = c(population$rowId[training], population$rowId[isTest]),
                       index = c(as.integer(trainingFolds), rep(-1L, sum(isTest)))),
            trainingN = length(training), trainingEvents = sum(outcome[training]), testN = sum(isTest),
            testEvents = sum(outcome[isTest]), poolN = sum(!isTest), poolEvents = sum(outcome[!isTest]))
}

trainingFolds <- function(y, nfold, seed) {
  nSplits <- min(nfold, sum(y), length(y) - sum(y))
  if (nSplits < 2) stop("fewer than 2 rows in the minority class of the training set")
  stratifiedFolds(y, nSplits, seed)
}

# Several entries per person (PLP's firstExposureOnly = FALSE, as in the COVID benchmarks): the rows of one person stay
# together, in the test split, in the ordering and in the folds, so that no person is on both sides of a split.
hasRepeatedPersons <- function(population) "subjectId" %in% names(population) && anyDuplicated(population$subjectId) > 0

# The pool's rows in the order of a seeded random order of its persons, a person's rows together.
orderByPerson <- function(pool, subjects, seed) {
  persons <- unique(subjects[pool])
  shuffled <- persons[withSeed(seed, sample.int(length(persons)))]
  pool[order(match(subjects[pool], shuffled))]
}

# Folds stratified on whether the person has any event; every row of a person in the person's fold.
personFolds <- function(y, subjects, nfold, seed) {
  personY <- tapply(y, subjects, max)
  folds <- trainingFolds(as.numeric(personY), nfold, seed)
  folds[match(as.character(subjects), names(personY))]
}

# The study split of one budget (NULL: the whole training pool) and ordering (1, 2, 3).
studySplit <- function(population, task, budget, ordering, nfold = 3) {
  repeated <- hasRepeatedPersons(population)
  settings <- PatientLevelPrediction::createDefaultSplitSetting(testFraction = 0.25, trainFraction = 0.75, nfold = 3,
                                                                splitSeed = STUDY_SEEDS$testSplit,
                                                                type = if (repeated) "subject" else "stratified")
  splitter <- utils::getFromNamespace(attr(settings, "fun"), "PatientLevelPrediction")
  plpSplit <- withSeed(STUDY_SEEDS$testSplit, splitter(population = population, splitSettings = settings))
  isTest <- population$rowId %in% plpSplit$rowId[plpSplit$index < 0]
  outcome <- as.integer(population$outcomeCount > 0)
  pool <- which(!isTest)
  ordered <- if (repeated) orderByPerson(pool, population$subjectId, STUDY_SEEDS$orderings[ordering])
             else pool[withSeed(STUDY_SEEDS$orderings[ordering], sample.int(length(pool)))]
  size <- if (is.null(budget)) length(ordered) else budgetPrefixLength(outcome[ordered], budget)
  training <- ordered[seq_len(size)]
  folds <- if (repeated) personFolds(outcome[training], population$subjectId[training], nfold, STUDY_SEEDS$folds)
           else trainingFolds(outcome[training], nfold, STUDY_SEEDS$folds)
  assembleSplit(population, training, folds, isTest)
}
