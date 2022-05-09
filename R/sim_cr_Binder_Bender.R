# block_corr_mat <- function(n = 100, p = 15, corr = c(0.5, 0.35, 0.05)) {
#   block_len <- p/length(corr)
#
#   corr_blocks <- lapply(corr, function(x) {
#     m <- matrix(rep(x, block_len^2), ncol = block_len)
#     diag(m) <- rep(1, block_len)
#     m
#   })
#
#   # Assemble into sparse block matrix
#   Sigma <- Matrix::bdiag(corr_blocks)
#
#   # generate observations
#   MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
#
# }

#' Generate block-correlated data as in Binder et. al. (2008) page 10.
#' @references https://www.fdm.uni-freiburg.de/publications-preprints/preprints/papers/pre100.pdf
#' "This results in correlations of about
#' 0.50 for j ≤ 0.05p,
#' 0.35 for 0.05p < j ≤ 0.1p,
#' 0.05 for 0.1p < j ≤ 0.2p,
#' 0.32 for 0.2p < j ≤ 0.3p,
#' and no correlation otherwise."
#' @param n Number of obs.
#' @param p Number of predictors
block_corr_binder <- function(n = 100, p = 100) {

  X <- matrix(0, nrow = n, ncol = p)

  # This is a shamelessly dumb copypasta-reimplementation of what's in the paper, sorry
  # using row(X) as i and col(X) as j
  cond1 <- which((row(X) <= 0.5 * n) & (col(X) <= 0.05 * p))
  X[cond1] <- -1 * rnorm(length(cond1))

  cond2 <- which((row(X) > 0.5 * n) & (col(X) <= 0.05 * p))
  X[cond2] <- 1 * rnorm(length(cond2))

  cond3 <- which((0.05 * p < col(X)) & (col(X) <= 0.1 * p))
  X[cond3] <- 1.5 * as.integer((runif(length(cond3))) < 0.4) + rnorm(length(cond3))

  cond4 <- which((0.1 * p < col(X)) & (col(X) <= 0.2 * p))
  X[cond4] <- 0.5 * as.integer((runif(length(cond4))) < 0.7) + rnorm(length(cond4))

  cond5 <- which((0.2 * p < col(X)) & (col(X) <= 0.3 * p))
  X[cond5] <- 1.5 * as.integer((runif(length(cond5))) < 0.3) + rnorm(length(cond5))

  cond6 <- which(col(X) > 0.3 * p)
  X[cond6] <- rnorm(length(cond6))


  X
}

if (FALSE) {
  X <- block_corr_binder(n = 1000, p = 1000)
  Xcorr <- cor(X)

  # First couple features should have correlations around 0.5
  Xcorr[1:5, 1:5]
}

