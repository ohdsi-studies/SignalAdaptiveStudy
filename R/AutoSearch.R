# R port of Cyclops's one-dimensional auto-search (Cyclops develop 91e0e49b):
#   src/utils/HParSearch.cpp                        UniModalSearch::step, QuadrLogFit
#   src/cyclops/drivers/AutoSearchCrossValidationDriver.cpp   doCrossValidationLoop (nDim = 1), refineSearch
#   src/cyclops/drivers/AbstractCrossValidationDriver.cpp     computeCommonFoldEstimates (no failed folds here)
# The caller supplies evaluate(variance) -> per-fold held-out log-likelihoods. Verified against native Cyclops by
# replaying native fold scores.
# Null-model bound (Cyclops PR #101, cvNullModelBound, on by default for LASSO under fitter = "super"): with
# nullVariance, the variance below which every penalised coefficient is 0, the search starts at no less than it and
# never evaluates below it (see autoSearch).
# C++ semantics kept: a comparison with NaN is false (isTRUE() here), and a singular quadratic fit gives NaN instead
# of an error (Eigen's solver does not throw), so a flat CV curve ends the search at the best value tried.
# Ties: Cyclops compares fold means exactly, and in a flat region (e.g. LASSO models that are all zero) its means tie
# exactly, so it moves right (up). Here the means come from separate fits and can differ in the last bit, so means
# within rounding (relative 1e-12) count as ties; ties then break as in Cyclops (2 points: move right; best: the first
# value tried, or the leftmost in the refinement).

TIE <- 1e-12
above <- function(a, b) a > b + TIE * max(abs(a), abs(b))   # a beats b by more than rounding
firstBest <- function(m) which(m >= max(m) - TIE * abs(max(m)))[1]

quadraticLogFit <- function(x, m, s) {
  positive <- !is.na(s) & s > 0
  zeroAdjust <- if (any(positive)) 100 / min(s[positive]) else 1
  weight <- ifelse(positive, 1 / s, zeroAdjust)
  design <- cbind(1, log(x), log(x)^2)
  coefficients <- tryCatch(solve(crossprod(design * weight, design), crossprod(design * weight, m)),
                           error = function(e) rep(NaN, 3))
  list(c0 = coefficients[1], c1 = coefficients[2], c2 = coefficients[3])
}

# x, m, s: tried values with fold mean and sd, in the order tried; best = first maximum in that order.
uniModalStep <- function(x, m, s, stdStep = 10, stopByY = 0.01, stopByX = log(1.5), firstCut = 0) {
  best <- x[firstBest(m)]
  bestM <- m[firstBest(m)]
  ord <- order(x)
  x <- x[ord]; m <- m[ord]; s <- s[ord]
  result <- list(more = TRUE, next_ = 0, expected = 0)
  if (length(x) == 1) {
    result$next_ <- if (x < firstCut) x * stdStep else x / stdStep
  } else if (length(x) == 2) {
    result$next_ <- if (above(m[1], m[2])) x[1] / stdStep else x[2] * stdStep
  } else if (x[1] == best) {
    result$next_ <- best / stdStep
    if (!(best > .Machine$double.xmin)) { result$more <- FALSE; result$expected <- bestM }
  } else if (x[length(x)] == best) {
    result$next_ <- best * stdStep
    if (!is.finite(best)) { result$more <- FALSE; result$expected <- bestM }
  } else {
    q <- quadraticLogFit(x, m, s)
    logArgmax <- -q$c1 / q$c2 / 2
    expectedMax <- -q$c1^2 / q$c2 / 4 + q$c0
    if (isTRUE(bestM == 0)) result$more <- FALSE
    else if (isTRUE((expectedMax - bestM) / abs(bestM) < stopByY)) result$more <- FALSE
    else if (isTRUE(abs(logArgmax - log(best)) < stopByX)) result$more <- FALSE
    result$next_ <- exp(logArgmax)
    result$expected <- expectedMax
  }
  result
}

# computeCommonFoldEstimates: a fold that failed (NaN) at any tried value is left out at every value, as long as some
# fold succeeded everywhere; means and standard deviations ignore the remaining NaNs.
foldEstimates <- function(scores) {
  failed <- Reduce(`|`, lapply(scores, is.nan))
  if (!all(failed)) scores <- lapply(scores, function(x) { x[failed] <- NaN; x })
  m <- vapply(scores, function(x) mean(x[!is.nan(x)]), numeric(1))
  s <- vapply(seq_along(scores), function(i) { x <- scores[[i]][!is.nan(scores[[i]])]; sqrt(mean(x^2) - m[i]^2) }, numeric(1))
  list(m = m, s = s)
}

# evaluate(variance) returns the per-fold held-out log-likelihoods at that variance.
# nullVariance (LASSO; NA otherwise): start at max(start, nullVariance); when a step falls below it, try nullVariance
# if it has not been tried; else, if it is the best so far and every value tried is below 100 x nullVariance, try 10x
# the largest value tried; else stop at the best value.
autoSearch <- function(evaluate, start, firstCut = 0, refine = TRUE, maxSteps = 50, nullVariance = NA) {
  bounded <- isTRUE(nullVariance > 0 && is.finite(nullVariance))
  if (bounded && start < nullVariance) start <- nullVariance
  tried <- numeric(); scores <- list(); current <- start; value <- NA_real_; step <- 0
  record <- function(v) {
    i <- match(v, tried)
    if (is.na(i)) { tried <<- c(tried, v); scores[[length(scores) + 1]] <<- evaluate(v) }
    else scores[[i]] <<- evaluate(v)
  }
  failed <- FALSE
  repeat {
    record(current)
    est <- foldEstimates(scores)
    scored <- !is.nan(est$m)
    step <- step + 1
    if (!scored[match(current, tried)]) { # no fold fit converged here: stop at the best value tried
      if (!any(scored)) stop("cross-validation failed: no fold fit converged")
      current <- tried[scored][firstBest(est$m[scored])]; value <- max(est$m[scored]); failed <- TRUE
      break
    }
    nxt <- uniModalStep(tried[scored], est$m[scored], est$s[scored], firstCut = firstCut)
    if (bounded && isTRUE(nxt$next_ < nullVariance)) {
      best <- tried[scored][firstBest(est$m[scored])]
      same <- function(a, b) abs(log(a) - log(b)) < 1e-9
      if (!any(same(tried, nullVariance))) {
        nxt <- list(more = TRUE, next_ = nullVariance, expected = 0)
      } else if (same(best, nullVariance) && max(tried) < 100 * nullVariance) { # the null model is best: look further up
        nxt <- list(more = TRUE, next_ = 10 * max(tried), expected = 0)
      } else {
        nxt <- list(more = FALSE, next_ = best, expected = max(est$m[scored]))
      }
    }
    if (isTRUE(nxt$next_ > 0 && is.finite(nxt$next_))) {
      current <- nxt$next_; value <- nxt$expected
      if (!nxt$more) break
    } else {
      current <- tried[scored][firstBest(est$m[scored])]; value <- max(est$m[scored]); break
    }
    if (step >= maxSteps) break
  }
  if (refine && !failed) {
    for (extra in 1:3) {
      if (any(abs(log(tried) - log(current)) < 1e-9)) break
      record(current); step <- step + 1
      est <- foldEstimates(scores)
      keep <- !is.nan(est$m)
      if (!any(keep)) break
      ord <- order(tried[keep]); x <- tried[keep][ord]; m <- est$m[keep][ord]; s <- est$s[keep][ord]
      b <- firstBest(m)
      if (b == 1 || b == length(x)) { current <- x[b]; value <- m[b]; break }
      local <- (b - 1):(b + 1)
      q <- quadraticLogFit(x[local], m[local], s[local])
      if (!isTRUE(q$c2 < 0)) { current <- x[b]; value <- m[b]; break }
      nextValue <- min(max(exp(-q$c1 / q$c2 / 2), x[b - 1]), x[b + 1])
      value <- -q$c1^2 / q$c2 / 4 + q$c0
      done <- abs(log(nextValue) - log(x[b])) < log(1.15)
      current <- nextValue
      if (done) break
    }
  }
  list(variance = current, ordinate = value, tried = tried, scores = scores, steps = step)
}
