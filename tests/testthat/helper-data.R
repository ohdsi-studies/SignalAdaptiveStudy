# Shared by the core and replication tests.

simulateDesign <- function(n = 1500, p = 40, seed = 1) {
  set.seed(seed)
  X <- Matrix::rsparsematrix(n, p, 0.1)
  X@x[] <- 1
  X <- methods::as(X, "CsparseMatrix")
  y <- stats::rbinom(n, 1, stats::plogis(-2 + as.numeric(X %*% c(stats::rnorm(10, 0, 0.8), rep(0, p - 10)))))
  list(X = X, y = y, folds = sample(rep_len(1:5, n)))
}

glmnetFit <- function(X, y, lambda, penalty) {
  keep <- Matrix::colSums(X) > 0 & Matrix::colSums(X) < nrow(X)
  fit <- glmnet::glmnet(as.matrix(X[, keep]), y, family = "binomial", alpha = if (penalty == "ridge") 0 else 1,
                        lambda = lambda, standardize = FALSE, control = list(thresh = 1e-14, maxit = 1e7))
  beta <- numeric(ncol(X))
  beta[keep] <- as.numeric(fit$beta)
  list(intercept = as.numeric(fit$a0), beta = beta)
}

simulatePopulation <- function(n = 6000, rate = 0.05, seed = 5) {
  set.seed(seed)
  data.frame(rowId = seq_len(n) + 100, subjectId = sample(1e6, n), outcomeCount = stats::rbinom(n, 1, rate))
}

# Properties every split must have: k events ending at the k-th event, stratified folds 1..5, a test set that is the
# same for every budget and ordering, and smaller budgets as prefixes of larger ones. splits: named "ordering budget".
expectSplitProperties <- function(splits, population, budgets, nfold = 5) {
  for (name in names(splits)) {
    s <- splits[[name]]
    k <- as.integer(strsplit(name, " ")[[1]][2])
    training <- s$rowId[s$index > 0]
    y <- population$outcomeCount[match(training, population$rowId)]
    expect_identical(attr(s, "trainingEvents"), k)
    expect_identical(sum(y), k)
    expect_equal(utils::tail(y, 1), 1)
    expect_setequal(unique(s$index[s$index > 0]), seq_len(nfold))
    perFold <- tapply(y, s$index[s$index > 0], sum)
    expect_lte(max(perFold) - min(perFold), 1)
    expect_false(anyDuplicated(s$rowId) > 0)
  }
  tests <- lapply(splits, function(s) sort(s$rowId[s$index == -1]))
  for (t in tests) expect_identical(t, tests[[1]])
  for (o in unique(sapply(strsplit(names(splits), " "), `[`, 1))) {
    a <- splits[[paste(o, min(budgets))]]; b <- splits[[paste(o, max(budgets))]]
    expect_identical(a$rowId[a$index > 0], b$rowId[b$index > 0][seq_len(sum(a$index > 0))])
  }
  invisible(tests[[1]])
}
