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
  n <- nrow(X)
  y <- as.numeric(y)
  target <- (sum(y) + 0.5) / (n + 1)
  mu <- Matrix::colMeans(X)
  variance <- pmax(Matrix::colMeans(X * X) - mu * mu, 0)
  score <- as.numeric(Matrix::crossprod(X, y - mean(y)))
  info <- n * target * (1 - target) * variance
  logWeights <- c[1] * log1p(n * mu) / log1p(n) + c[2] * log1p(info) + c[3] * log1p(score^2 / (info + 0.5))
  # not constant in these rows; relative to the mean square, so that rounding in a constant column (about 1e-16 of it)
  # does not count, while a binary column's variance / mean square is 1 - mu >= 1/n
  keep <- variance > 1e-10 * Matrix::colMeans(X * X)
  logWeights <- logWeights - mean(logWeights[keep])
  logWeights[!keep] <- NA
  list(logWeights = logWeights, weights = exp(logWeights), keep = keep)
}
