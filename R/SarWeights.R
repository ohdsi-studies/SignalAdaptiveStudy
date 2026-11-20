# Ported from the original SAR code (signal-adaptive/evaluation/control_cv_439.py, log_weights, sha256 53acf5df...,
# and control_cv_438.py, weights_for). Port only: the formula must not change. The kept-column rule is the original's
# for binary columns (0 < mean < 1, i.e. not constant); for other columns it is "not constant" too (variance > 0;
# investigator, 2026-11-01), where the original's rule would drop any column with mean >= 1, such as age in years.

#' Signal-adaptive penalty weights
#'
#' Computes the SAR/SAL per-feature penalty weights from the rows given, which must be the training rows of the
#' current fit (a CV training fold, or the whole training set for the final fit).
#'
#' With n rows, e events, \eqn{\mu_j} the column mean, \eqn{t = (e + 0.5)/(n + 1)},
#' \eqn{h_j = n t (1 - t) \mathrm{var}(x_j)} and \eqn{s_j = x_j^\top (y - \bar{y})}:
#' \deqn{\log w_j = c_1 \frac{\log(1 + n\mu_j)}{\log(1 + n)} + c_2 \log(1 + h_j) + c_3 \log\left(1 + \frac{s_j^2}{h_j + 0.5}\right)}
#' minus its mean over the kept columns. Kept columns are those that are not constant in these rows (variance > 0;
#' for a binary column \eqn{0 < \mu_j < 1}); the others are left out of the fit and get coefficient 0.
#'
#' @param X A numeric matrix or `Matrix` (rows = training rows, columns = features), binary in the study.
#' @param y Outcome vector of 0/1, one per row of `X`.
#' @param c Numeric vector `c(c1, c2, c3)`.
#'
#' @return A list with `logWeights` and `weights` (length `ncol(X)`, `NA` for columns not kept) and the logical
#'   vector `keep`.
#' @export
computeSarWeights <- function(X, y, c) {
  if (length(c) != 3) stop("c must be c(c1, c2, c3)")
  if (nrow(X) != length(y)) stop("X and y have different numbers of rows")
  weightsFromStatistics(sarStatistics(X, y), c)
}

# What the weights of every shape share, computed once per fit problem: n, the column means, variances and scores. The
# mean squares come from the stored values of a sparse column matrix, which is far cheaper than X * X.
sarStatistics <- function(X, y) {
  n <- nrow(X)
  y <- as.numeric(y)
  target <- (sum(y) + 0.5) / (n + 1)
  mu <- Matrix::colMeans(X)
  square <- X
  if (methods::is(square, "dgCMatrix")) square@x <- square@x^2 else square <- X * X
  meanSquare <- Matrix::colMeans(square)
  variance <- pmax(meanSquare - mu * mu, 0)
  list(n = n, mu = mu, variance = variance, meanSquare = meanSquare, score = as.numeric(Matrix::crossprod(X, y - mean(y))),
       info = n * target * (1 - target) * variance)
}

weightsFromStatistics <- function(s, c, weightRule = NULL) {
  logWeights <- if (is.null(weightRule)) {
    # a term with c = 0 is left out rather than multiplied by 0 (the same weights; on columns that are not binary, e.g.
    # uniLasso's features, a term may be NaN)
    term <- function(cj, value) if (cj == 0) numeric(length(s$mu)) else cj * value
    term(c[1], log1p(s$n * s$mu) / log1p(s$n)) + term(c[2], log1p(s$info)) + term(c[3], log1p(s$score^2 / (s$info + 0.5)))
  } else {
    posteriorMomentLogWeights(s, weightRule$kappa)
  }
  # not constant in these rows; relative to the mean square, so that rounding in a constant column (about 1e-16 of it)
  # does not count, while a binary column's variance / mean square is 1 - mu >= 1/n
  keep <- s$variance > 1e-10 * s$meanSquare
  logWeights <- logWeights - mean(logWeights[keep])
  logWeights[!keep] <- NA
  list(logWeights = logWeights, weights = exp(logWeights), keep = keep)
}

# Posterior-moment ridge (PM-ridge; a secondary method, added 2026-11-13 after a pilot on development tasks). Under a N(0, 1/kappa) prior and the
# diagonal quadratic approximation of the log-likelihood at the intercept-only model (score s_j, information h_j, the
# same statistics as SAR), the posterior of beta_j has mean s_j / D_j and variance 1 / D_j, D_j = h_j + kappa. The
# weight is the inverse of the posterior second moment, u_j = D_j^2 / (D_j + s_j^2): one EM update of beta_j's prior
# variance. Centred in the log over the same columns as SAR. With kappa = 0.5 its ratio to SAR's c = (0, 1, -1) before
# centring is (h_j + 0.5) / (h_j + 1).
posteriorMomentLogWeights <- function(s, kappa) {
  if (!is.numeric(kappa) || length(kappa) != 1 || !is.finite(kappa) || kappa <= 0) stop("kappa must be positive and finite")
  d <- s$info + kappa
  log(d) - log1p(s$score^2 / d)
}

# NULL (SAR's formula) or list(type = "posteriorMoment", kappa); checked once, where the settings are made.
checkWeightRule <- function(weightRule) {
  if (is.null(weightRule)) return(NULL)
  if (!is.list(weightRule) || !identical(weightRule$type, "posteriorMoment")) stop("unknown weight rule")
  posteriorMomentLogWeights(list(info = 0, score = 0), weightRule$kappa)
  list(type = "posteriorMoment", kappa = as.numeric(weightRule$kappa))
}

weightRuleLabel <- function(weightRule) {
  if (is.null(weightRule)) "SAR" else sprintf("posteriorMoment(kappa=%g)", weightRule$kappa)
}
