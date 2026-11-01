# Cross-validation for SAR and SAL (option B: our own cross-validation): our folds, weights recomputed from each
# training fold's rows, Cyclops as the solver of one fit at a time.
#
# Penalties (mean logistic loss + penalty, unpenalized intercept, raw columns; columns constant in the fit's rows
# are left out with coefficient 0):
#   SAR (ridge)  (lambda / 2) * sum_j w_j beta_j^2
#   SAL (lasso)  lambda * sum_j sqrt(w_j) |beta_j|
# Cyclops minimises -sum loglik + prior, so with n rows in the fit the per-covariate prior variance is
#   ridge  v_j = 1 / (n lambda w_j)          (normal: beta^2 / (2 v))
#   lasso  v_j = 2 / (n^2 lambda^2 w_j)      (Laplace: rate sqrt(2 / v))
# and, in the Cyclops convention, v_j = sigma^2 / w_j for both.
#
# The penalty is chosen per shape by a search strategy. The study uses "auto": Cyclops's auto-search over sigma^2
# (R/AutoSearch.R), with the same sigma^2 on every fold and on the final fit, as in Cyclops.

# Cyclops data for the rows of X (binary columns become indicator columns).
makeCyclopsData <- function(X, y) {
  triplets <- Matrix::summary(methods::as(X, "TsparseMatrix"))
  triplets <- triplets[order(triplets$i, triplets$j), ]
  Cyclops::convertToCyclopsData(
    outcomes = data.frame(rowId = seq_len(nrow(X)), y = y),
    covariates = data.frame(rowId = triplets$i, covariateId = triplets$j, covariateValue = triplets$x),
    modelType = "lr", addIntercept = TRUE, quiet = TRUE
  )
}

# One fit problem: the Cyclops data of the fit rows and the columns it holds (indices into X).
createFitProblem <- function(X, y, rows) {
  data <- makeCyclopsData(X[rows, , drop = FALSE], y[rows])
  list(rows = rows, n = sum(rows), data = data, ids = as.integer(Cyclops::getCovariateIds(data)[-1]))
}

# Weights of one shape for a fit problem, from that problem's rows only, restricted to the columns it holds. The SAR
# weights (and their centring) cover the penalised columns; an unpenalised column gets weight 0, i.e. an infinite
# prior variance (a flat prior), and is kept when it varies in these rows.
problemWeights <- function(X, y, problem, shape, unpenalized = rep(FALSE, ncol(X))) {
  penalized <- !unpenalized
  w <- computeSarWeights(X[problem$rows, penalized, drop = FALSE], y[problem$rows], shape)
  weights <- rep(NA_real_, ncol(X)); keep <- logical(ncol(X))
  weights[penalized] <- w$weights; keep[penalized] <- w$keep
  if (any(unpenalized)) {
    free <- X[problem$rows, unpenalized, drop = FALSE]
    varies <- Matrix::colMeans(free * free) - Matrix::colMeans(free)^2 > 0
    keep[unpenalized] <- varies
    weights[unpenalized] <- ifelse(varies, 0, NA)
  }
  list(weights = weights[problem$ids], keep = keep[problem$ids], full = list(weights = weights, keep = keep))
}

priorVariances <- function(penalty, weights, keep, scale, scaleType, n) {
  base <- switch(scaleType,
    variance = scale,
    lambda = if (penalty == "ridge") 1 / (n * scale) else 2 / (n^2 * scale^2)
  )
  v <- base / weights
  v[!keep] <- 1 # these coefficients are fixed at 0
  v
}

varianceToLambda <- function(penalty, variance, n) if (penalty == "ridge") 1 / (n * variance) else sqrt(2 / variance) / n

# Mean-loss KKT residual on the original scale (the original code's certificate).
kktResidual <- function(X, y, problem, coefficients, variances, keep, penalty) {
  Xs <- X[problem$rows, problem$ids, drop = FALSE]
  b <- coefficients[-1]
  r <- stats::plogis(coefficients[1] + as.numeric(Xs %*% b)) - y[problem$rows]
  g <- as.numeric(Matrix::crossprod(Xs, r)) / problem$n
  if (penalty == "ridge") {
    v <- abs(g + b / (problem$n * variances))
  } else {
    rate <- sqrt(2 / variances) / problem$n
    v <- ifelse(abs(b) > 1e-8, abs(g + rate * sign(b)), pmax(abs(g) - rate, 0))
  }
  max(abs(mean(r)), v[keep], 0)
}

fitProblem <- function(problem, penalty, variances, keep, start, control) {
  type <- if (penalty == "ridge") "normal" else "laplace"
  free <- is.infinite(variances) # unpenalised (weight 0)
  prior <- Cyclops::createPrior(c("none", ifelse(free, "none", type)), variance = c(0, ifelse(free, 0, variances)))
  fixed <- NULL
  if (any(!keep)) {
    fixed <- c(FALSE, !keep)
    if (is.null(start)) start <- numeric(length(problem$ids) + 1)
    start[fixed] <- 0
  }
  fit <- Cyclops::fitCyclopsModel(problem$data, prior = prior, control = control, startingCoefficients = start,
                                  fixedCoefficients = fixed, warnings = FALSE)
  coefficients <- if (is.null(fit$estimation)) NULL else as.numeric(stats::coef(fit, ignoreConvergence = TRUE))
  # converged as Cyclops's own CV counts it: SUCCESS, or POOR_BLR_STEP under the kkt criterion
  list(coefficients = coefficients, status = fit$return_flag,
       converged = !is.null(coefficients) && fit$return_flag %in% c("SUCCESS", "POOR_BLR_STEP"))
}

softplus <- function(z) ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))

heldOutLogLikelihood <- function(X, y, problem, coefficients) {
  heldOut <- !problem$rows
  eta <- coefficients[1] + as.numeric(X[heldOut, problem$ids, drop = FALSE] %*% coefficients[-1])
  -sum((1 - y[heldOut]) * softplus(eta) + y[heldOut] * softplus(-eta))
}

smoothedIntercept <- function(y) log((sum(y) + 0.5) / (length(y) - sum(y) + 0.5))

# The variance below which every penalised coefficient of the LASSO is 0 (Cyclops PR #101's null-model variance,
# 2 / lambdaMax^2 without weights). At the fit of the unpenalised coefficients alone (the intercept and any unpenalised
# columns) on the problem's rows, coefficient j stays 0 while |g_j| <= sqrt(2 w_j / sigma^2), g the log-likelihood
# gradient, so all are 0 for sigma^2 <= min_j 2 w_j / g_j^2. NA when no penalised column has a gradient.
nullModelVariance <- function(X, y, problem, setup, control) {
  free <- setup$keep & setup$weights == 0
  penalised <- setup$keep & !free
  if (!any(penalised)) return(NA_real_)
  fit <- fitProblem(problem, "lasso", priorVariances("lasso", setup$weights, free, 1, "variance", problem$n), free, NULL,
                    control)
  if (!fit$converged) return(NA_real_)
  Xs <- X[problem$rows, problem$ids, drop = FALSE]
  g <- as.numeric(Matrix::crossprod(Xs, y[problem$rows] - stats::plogis(fit$coefficients[1] +
                                                                           as.numeric(Xs %*% fit$coefficients[-1]))))
  bound <- 2 * setup$weights[penalised] / g[penalised]^2
  if (!any(is.finite(bound))) NA_real_ else min(bound[is.finite(bound)])
}

# Cyclops's default starting variance (without weights): non-zeros / (columns x rows), in double precision (the
# integer product overflows above 2^31 cells, e.g. 9,184 x 282,852).
defaultStartVariance <- function(X, columns) length(X@x) / (as.numeric(columns) * nrow(X))

shapeKey <- function(shape) c(sum(shape != 0), sum(abs(shape)))

# ---------------------------------------------------------------- search strategies
# A strategy has: searchShape(X, y, problems, setups, shape, penalty, context) -> list(record, trace, state);
# penaltyScale(records) (larger = stronger penalty, for ties); finalFit(X, y, full, penalty, w, record) ->
# list(fit, variances, lambda, variance); foldCoefficients(X, y, problem, f, penalty, setup, record, state).
searchStrategy <- function(search, settings = list()) {
  switch(search,
    auto = autoStrategy(),
    stop("unknown search: ", search)
  )
}

# Cyclops's auto-search over sigma^2 (fold fits as Cyclops's own CV under fitter = "super"), final fit at the same
# sigma^2 on all rows.
# One thread per fit (Cyclops's default, explicit): execute() runs fits in parallel in worker processes instead.
autoStrategy <- function() {
  foldControl <- Cyclops::createControl(fitter = "super", convergenceType = "kkt", tolerance = 3e-3, maxIterations = 5000,
                                        threads = 1)
  # Cyclops's final-fit settings (super: kkt, 1e-4) with more iterations than its default 5000, which a weighted
  # LASSO on 283k rows exceeded; a fit that converges sooner stops at the same point either way.
  finalControl <- Cyclops::createControl(fitter = "super", maxIterations = 100000, threads = 1)
  list(
    searchShape = function(X, y, problems, setups, shape, penalty, context) {
      n <- length(y)
      cache <- lapply(problems, function(p) list(variances = numeric(), coefficients = list()))
      evaluate <- function(sigma2) {
        vapply(seq_along(problems), function(f) {
          p <- problems[[f]]; s <- setups[[f]]; ca <- cache[[f]]
          warm <- if (length(ca$variances)) ca$coefficients[[which.min(abs(log(ca$variances) - log(sigma2)))]] else NULL
          fit <- fitProblem(p, penalty, priorVariances(penalty, s$weights, s$keep, sigma2, "variance", p$n), s$keep,
                            warm, foldControl)
          if (!fit$converged) return(NaN) # as Cyclops: not scored and not used as a warm start
          cache[[f]]$variances <<- c(ca$variances, sigma2)
          cache[[f]]$coefficients[[length(ca$variances) + 1]] <<- fit$coefficients
          heldOutLogLikelihood(X, y, p, fit$coefficients)
        }, numeric(1))
      }
      nullVariance <- if (penalty == "lasso") {
        nullModelVariance(X, y, context$full, problemWeights(X, y, context$full, shape, context$unpenalized), finalControl)
      } else NA_real_
      search <- autoSearch(evaluate, start = context$startVariance, firstCut = 0, refine = TRUE,
                           nullVariance = nullVariance)
      chosen <- search$variance
      i <- which(abs(log(search$tried) - log(chosen)) < 1e-12)
      observed <- if (length(i)) search$scores[[i[1]]] else evaluate(chosen) # observed score at the chosen variance
      if (anyNA(observed)) stop("a fold fit did not converge at the selected variance")
      list(
        record = data.frame(c1 = shape[1], c2 = shape[2], c3 = shape[3], lambda = NA_real_, variance = chosen,
                            cvLoss = -sum(observed) / n, ordinateLoss = -search$ordinate * length(problems) / n,
                            nullVariance = nullVariance),
        trace = data.frame(c1 = shape[1], c2 = shape[2], c3 = shape[3], step = seq_along(search$tried),
                           variance = search$tried, meanFoldLogLikelihood = vapply(search$scores, mean, numeric(1))),
        state = cache
      )
    },
    penaltyScale = function(records) -records$variance,
    finalFit = function(X, y, full, penalty, w, record) {
      v <- priorVariances(penalty, w$weights, w$keep, record$variance, "variance", full$n)
      fit <- fitProblem(full, penalty, v, w$keep, NULL, finalControl)
      if (!fit$converged) stop("the final fit did not converge (", fit$status, ")")
      list(fit = fit, variances = v, variance = record$variance, lambda = varianceToLambda(penalty, record$variance, full$n))
    },
    foldCoefficients = function(X, y, problem, f, penalty, setup, record, state) {
      state[[f]]$coefficients[[which.min(abs(log(state[[f]]$variances) - log(record$variance)))]]
    }
  )
}

# ---------------------------------------------------------------- driver
# Cross-validated SAR / SAL fit: selects c and the penalty on the given folds, with weights recomputed from each
# training fold's rows, and fits the selected model on all rows. Returns the intercept and coefficients (original
# scale, one per column of X), the selection, the full-training weights, the CV records, the search trace and the
# held-out CV predictions of the selected model.
#   cSelection: choose c by the observed pooled CV log loss ("observed") or by the search's predicted optimum
#   ("ordinate"; auto only). Both are recorded.
#   interceptOnly: the intercept-only model wins if its CV loss is not worse (the original's rule).
#   unpenalized: logical, one per column of X: columns left out of the penalty (and of the weights' centring).
signalAdaptiveCv <- function(X, y, folds, penalty = c("ridge", "lasso"), shapes = sarShapes(), interceptOnly = TRUE,
                             search = "auto", cSelection = c("observed", "ordinate"), searchSettings = list(),
                             unpenalized = rep(FALSE, ncol(X))) {
  penalty <- match.arg(penalty); cSelection <- match.arg(cSelection)
  strategy <- searchStrategy(search, searchSettings)
  start <- Sys.time()
  X <- methods::as(methods::as(X, "CsparseMatrix"), "generalMatrix")
  y <- as.numeric(y)
  shapes <- normalizeShapes(shapes)
  if (min(sum(y), length(y) - sum(y)) < 2) stop("fewer than 2 rows in the minority class")
  if (length(unpenalized) != ncol(X) || anyNA(unpenalized)) stop("unpenalized needs one TRUE/FALSE per column of X")
  folds <- as.integer(folds)
  problems <- lapply(sort(unique(folds)), function(f) createFitProblem(X, y, folds != f))
  full <- createFitProblem(X, y, rep(TRUE, length(y)))
  interceptLoss <- sum(vapply(problems, function(p) {
    a <- smoothedIntercept(y[p$rows]); ho <- !p$rows
    sum((1 - y[ho]) * softplus(a) + y[ho] * softplus(-a))
  }, numeric(1))) / length(y)
  context <- list(startVariance = defaultStartVariance(X, length(full$ids)), full = full, unpenalized = unpenalized)

  # progress in the log (PLP's plpLog.txt of the run): one line per shape
  ParallelLogger::logInfo(sprintf("%s cross-validation: %d rows, %d events, %d columns, %d folds, %d shape(s)",
                                  penalty, length(y), sum(y), ncol(X), length(problems), nrow(shapes)))
  records <- list(); traces <- list(); states <- list()
  for (s in seq_len(nrow(shapes))) {
    shape <- unname(shapes[s, ])
    setups <- lapply(problems, function(p) problemWeights(X, y, p, shape, unpenalized))
    result <- strategy$searchShape(X, y, problems, setups, shape, penalty, context)
    records[[s]] <- result$record; traces[s] <- list(result$trace); states[s] <- list(result$state) # keep NULL slots
    best <- result$record[which.min(result$record$cvLoss), ] # one row (auto) or one per lambda (grid)
    chosen <- if (is.na(best$lambda)) sprintf("variance %.4g", best$variance) else sprintf("lambda %.4g", best$lambda)
    ParallelLogger::logInfo(sprintf("shape %d/%d (c = %s): %s, CV loss %.6f; %.1f min so far", s, nrow(shapes),
                                    paste(shape, collapse = ", "), chosen, best$cvLoss,
                                    as.numeric(difftime(Sys.time(), start, units = "mins"))))
  }
  records <- do.call(rbind, records)
  score <- if (cSelection == "ordinate") records$ordinateLoss else records$cvLoss
  if (anyNA(score)) stop("this search records no ", cSelection, " score")
  best <- min(score)
  selected <- NULL
  if (!(interceptOnly && interceptLoss <= best + 1e-10)) {
    tied <- which(score <= best + 1e-10)
    scale <- strategy$penaltyScale(records)
    keys <- t(vapply(tied, function(i) c(shapeKey(unlist(records[i, c("c1", "c2", "c3")])), -scale[i]), numeric(3)))
    selected <- tied[do.call(order, as.data.frame(keys))[1]]
  }

  p <- ncol(X); beta <- numeric(p); weights <- rep(NA_real_, p); kkt <- NA_real_; cvPrediction <- numeric(length(y))
  if (is.null(selected)) {
    intercept <- smoothedIntercept(y)
    for (pr in problems) cvPrediction[!pr$rows] <- stats::plogis(smoothedIntercept(y[pr$rows]))
    shape <- c(NA, NA, NA); lambda <- NA_real_; variance <- NA_real_
  } else {
    record <- records[selected, ]
    shape <- unname(unlist(record[c("c1", "c2", "c3")]))
    w <- problemWeights(X, y, full, shape, unpenalized)
    final <- strategy$finalFit(X, y, full, penalty, w, record)
    lambda <- final$lambda; variance <- final$variance
    kkt <- kktResidual(X, y, full, final$fit$coefficients, final$variances, w$keep, penalty)
    intercept <- final$fit$coefficients[1]
    beta[full$ids] <- final$fit$coefficients[-1]
    weights <- w$full$weights
    state <- states[[which(apply(shapes, 1, function(r) all(r == shape)))[1]]]
    for (f in seq_along(problems)) {
      pr <- problems[[f]]
      co <- strategy$foldCoefficients(X, y, pr, f, penalty, problemWeights(X, y, pr, shape, unpenalized), record,
                                      state)
      cvPrediction[!pr$rows] <- stats::plogis(co[1] + as.numeric(X[!pr$rows, pr$ids, drop = FALSE] %*% co[-1]))
    }
  }
  ParallelLogger::logInfo(sprintf("final fit done: %s; %.1f min in all",
                                  if (is.null(selected)) "intercept only" else paste("c =", paste(shape, collapse = ", ")),
                                  as.numeric(difftime(Sys.time(), start, units = "mins"))))
  list(
    intercept = unname(intercept), beta = beta, penalty = penalty, search = search, cSelection = cSelection,
    selected = list(shape = shape, lambda = lambda, variance = variance, interceptOnly = is.null(selected),
                    cvLoss = if (is.null(selected)) interceptLoss else records$cvLoss[selected],
                    ordinateLoss = if (is.null(selected)) NA_real_ else records$ordinateLoss[selected]),
    interceptCvLoss = interceptLoss, weights = weights, kkt = kkt, nonZero = sum(beta != 0),
    records = records, trace = if (length(Filter(Negate(is.null), traces))) do.call(rbind, traces) else NULL,
    cvPrediction = cvPrediction, folds = folds,
    seconds = as.numeric(difftime(Sys.time(), start, units = "secs"))
  )
}

#' The 12 SAR/SAL shapes
#'
#' c1, c2 in {0, 1} and c3 in {0, -0.5, -1}, as in the original code.
#' @return A matrix with columns c1, c2, c3.
#' @export
sarShapes <- function() {
  g <- expand.grid(c3 = c(0, -0.5, -1), c2 = c(0, 1), c1 = c(0, 1))
  as.matrix(g[, c("c1", "c2", "c3")])
}

normalizeShapes <- function(shapes) {
  if (is.list(shapes) && !is.data.frame(shapes)) shapes <- do.call(rbind, shapes)
  matrix(as.numeric(as.matrix(shapes)), ncol = 3, dimnames = list(NULL, c("c1", "c2", "c3")))
}
