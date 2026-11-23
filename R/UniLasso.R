# uniLasso (Chatterjee, Hastie & Tibshirani 2025, arXiv 2501.18360; CRAN package uniLasso 2.11), for binary covariates
# and a binary outcome, without dense matrices.
#
# Stage 1 (the package's uniInfo(family = "binomial", loo = TRUE, nit = 2, ridge = 0), unweighted): each covariate's
# univariate logistic fit, two IRLS steps from mu = (y + 0.5) / 2, and its leave-one-out ("prevalidated") fitted values F.
# The package builds several dense n x p matrices; with a binary covariate every quantity depends only on the covariate's
# 2 x 2 table (covariate present or not, outcome or not), so each column's fits, leverages and F take four values.
#
# Stage 2 (the package's glmnet(F, y, lower.limits = 0, standardize = FALSE)): a LASSO of y on F with coefficients >= 0.
# F_ij = A_j + G_ij, with A_j the value for a person without the covariate and without the outcome; A is absorbed by the
# unpenalised intercept, and G is sparse: non-zero where the covariate is present and in the outcome rows. The final
# model is linear in x (uniLassoModel()).

# Stage 1 for a binary sparse matrix X (n x p, values 1) and y in {0, 1}. Returns beta and beta0 (the univariate fits,
# 0 for a constant column, as the package), shift (A_j above) and G (n x p, sparse), with F = G + 1 shift'. A covariate in
# exactly one person (or in all but one) has leverage 1 in its group, where the package's F is not finite; there F is
# the fitted value instead (as the package does where F is NaN).
uniLassoFeatures <- function(X, y) {
  X <- methods::as(X, "CsparseMatrix")
  if (!all(X@x == 1)) stop("uniLassoFeatures() needs binary covariates (values 1)")
  y <- as.numeric(y)
  if (!all(y %in% c(0, 1)) || length(unique(y)) < 2) stop("uniLassoFeatures() needs a binary outcome with both values")
  n <- nrow(X); p <- ncol(X)
  n1 <- unname(Matrix::colSums(X)); e1 <- as.numeric(Matrix::crossprod(X, y)); n0 <- n - n1; e0 <- sum(y) - e1
  constant <- n1 == 0 | n1 == n
  # first IRLS step: working response +-(log 3 + 4/3) at mu = 0.75 / 0.25, equal weights: the group means
  zOne <- log(3) + 4 / 3
  groupMean <- function(e, m) ifelse(m > 0, zOne * (2 * e - m) / pmax(m, 1), 0)
  eta1 <- groupMean(e1, n1); eta0 <- groupMean(e0, n0)
  # second step: weights mu (1 - mu) and working response per group, then the group means of the working response
  step <- function(eta, e, m) {
    mu <- 1 / (1 + exp(-eta)); w <- mu * (1 - mu)
    list(w = w, z0 = eta - mu / w, z1 = eta + (1 - mu) / w, eta = ifelse(m > 0, eta + (e / pmax(m, 1) - mu) / w, eta))
  }
  g1 <- step(eta1, e1, n1); g0 <- step(eta0, e0, n0)
  beta <- ifelse(constant, 0, g1$eta - g0$eta)
  beta0 <- ifelse(constant, 0, g0$eta)
  # leave-one-out fits: leverage per group (weighted least squares on the indicator), then (eta - h z) / (1 - h)
  s1 <- n1 * g1$w / n; s0 <- n0 * g0$w / n; sW <- s1 + s0
  xbar <- s1 / sW
  sWxx <- s1 * (1 - xbar)^2 + s0 * xbar^2
  h1 <- g1$w / n / sW + g1$w / n * (1 - xbar)^2 / sWxx
  h0 <- g0$w / n / sW + g0$w / n * xbar^2 / sWxx
  loo <- function(g, h, z) {
    f <- (g$eta - h * z) / (1 - h)
    ifelse(is.finite(f) & abs(1 - h) > 1e-12, f, g$eta)
  }
  f00 <- loo(g0, h0, g0$z0); f01 <- loo(g0, h0, g0$z1); f10 <- loo(g1, h1, g1$z0); f11 <- loo(g1, h1, g1$z1)
  zero <- function(v) ifelse(constant, 0, v)
  shift <- zero(f00)
  # G: where the covariate is present, F - shift for the person's outcome; in the outcome rows, F01 - shift elsewhere.
  # Built on X's own column-compressed structure (no triplet copies).
  j <- rep.int(seq_len(p), diff(X@p)); i <- X@i + 1L
  event <- y[i] == 1
  present <- X
  present@x <- ifelse(event, zero(f11)[j] - zero(f01)[j], zero(f10)[j] - shift[j])
  rm(j, i, event)
  events <- which(y == 1)
  absentInEvents <- Matrix::sparseMatrix(i = rep(events, p), j = rep(seq_len(p), each = length(events)),
                                         x = rep(zero(f01) - shift, each = length(events)), dims = c(n, p))
  G <- present + absentInEvents
  rm(present, absentInEvents)
  dimnames(G) <- dimnames(X)
  list(beta = beta, beta0 = beta0, shift = shift, G = methods::as(G, "CsparseMatrix"), constant = constant)
}

# The final model from stage 2's fit on G (intercept on G's scale, coefficients theta >= 0): as the package's, an
# intercept and coefficients on x, gamma_j = theta_j beta_j, intercept theta0 + sum_j theta_j (beta0_j - shift_j).
uniLassoModel <- function(features, intercept, theta) {
  list(intercept = unname(intercept) + sum(theta * (features$beta0 - features$shift)), coefficients = unname(theta * features$beta))
}

# The fit problems of uniLasso's cross-validation, with stage 1 inside each fold (investigator, 2026-11-23): the package's
# cv.uniLasso computes the stage-1 features once on all rows and then cross-validates stage 2, so a held-out row's
# feature already holds a trace of its own outcome (the leave-one-out of two IRLS steps is approximate: the first step
# includes the person), which with many rare binary covariates makes the CV loss optimistic. Here each fold's stage 1
# uses the fold's training rows only (the package's leave-one-out features for them), and its held-out rows get the
# features of a new person, beta0 + beta x from that fit, as at test time. The full problem is stage 1 on all rows, as the
# package's final model. Returns list(G, features, prepared) for signalAdaptiveCv(G, ..., prepared = prepared).
uniLassoProblems <- function(X, y, folds) {
  X <- methods::as(X, "CsparseMatrix")
  y <- as.numeric(y)
  problem <- function(rows) {
    f <- uniLassoFeatures(X[rows, , drop = FALSE], y[rows])
    pr <- list(rows = rows, n = sum(rows), data = makeCyclopsData(f$G, y[rows]),
               statistics = sarStatistics(f$G, y[rows]))
    pr$ids <- as.integer(Cyclops::getCovariateIds(pr$data)[-1])
    if (any(!rows)) {
      # new-person features, beta0 + beta x, relative to the training rows' shift (absorbed by the intercept there)
      pr$heldOutX <- X[!rows, pr$ids, drop = FALSE] %*% Matrix::Diagonal(x = f$beta[pr$ids])
      pr$heldOutShift <- (f$beta0 - f$shift)[pr$ids]
    }
    list(problem = pr, features = f)
  }
  problems <- lapply(sort(unique(folds)), function(k) problem(folds != k)$problem)
  full <- problem(rep(TRUE, length(y)))
  list(G = full$features$G, features = full$features,
       prepared = list(problems = problems, full = full$problem, unpenalized = rep(FALSE, ncol(X)), shapeResults = list()))
}
