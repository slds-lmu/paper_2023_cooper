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
#' @param n Number of obs. (must be even).
#' @param p Number of predictors.
block_corr_binder <- function(n = 50, p = 1000) {
  stopifnot("n must be even" = n %% 2 == 0)
  stopifnot("p should be > 20" = 0.05 * p > 1)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  ui1 <- runif(n)
  ui2 <- runif(n)
  ui3 <- runif(n)
  j_seq <- seq_len(p)

  block1 <- which(j_seq <= 0.05 * p)
  X[seq_len(n/2), block1] <- -1 + X[seq_len(n/2), block1]
  X[-seq_len(n/2), block1] <- 1 + X[-seq_len(n/2), block1]


  block2 <- which((j_seq > (0.05 * p)) & (j_seq <= (0.1 * p)))
  X[, block2] <- 1.5 * (ui1 < 0.4) + X[, block2]


  block3 <- which((0.1 * p < j_seq) & (j_seq <= 0.2 * p))
  X[, block3] <- 0.5 * (ui2 < 0.7) + X[, block3]

  block4 <- which((0.2 * p < j_seq) & (j_seq <= 0.3 * p))
  X[, block4] <- 1.5 * (ui3 < 0.3) + X[, block4]

  X
}

#' Take a covariate matrix and simulate a survival outcome, based on the same paper.
#' @param n,p Passed to block_corr_binder().
#' @param ce Effect constant. Used as magnitude of effect.
#' @param lambda,lambda_c Constants for Cox-exponential model for survival & censoring times.
#'
#' Binder et al use ce in {0.05, 0.075, 0.1} for the censored survival setting.
#' @note Binder et. al. describe U = runif(n), but log(runif(n)) is used to ensure results as stated
#' with mean baseline survival time of around 10, rather than negative survival times.
#' This is also consistent with other descriptions of the Cox exponential model.
sim_surv_binder <- function(n = 50, p = 1000, ce = 0.05, lambda = 0.1, lambda_c = 0.1) {
  checkmate::assert_integerish(n, lower = 10, len = 1)
  checkmate::assert_true(n %% 2 == 0)
  checkmate::assert_integerish(p, lower = 400, len = 1)
  checkmate::assert_numeric(ce, lower = 0.01, len = 1)
  checkmate::assert_numeric(lambda, lower = 0.01, len = 1)
  checkmate::assert_numeric(lambda_c, lower = 0.01, len = 1)

  X <- block_corr_binder(n = n, p = p)
  beta <- rep(0, p)

  # j * 200/p: if odd -> ce, if even -> -ce, else 0
  tmp <- ((seq_len(p) * 200) / p) %% 2
  j_odd <- which(tmp == 1)
  j_even <- which(tmp == 0)

  checkmate::assert_integer(j_odd, null.ok = FALSE, min.len = 1)
  checkmate::assert_integer(j_even, null.ok = FALSE, min.len = 1)

  beta[j_odd] <- ce
  beta[j_even] <- -ce
  nonzero <- sum(beta > 0)
  checkmate::assert_count(nonzero, positive = TRUE, null.ok = FALSE)

  lp <- X %*% beta

  # Survival and censoring times
  Ti <- -log(runif(n)) / (lambda * exp(lp))
  Ci <- -log(runif(n)) / lambda_c

  ti <- pmin(Ti, Ci)
  di <- as.integer(Ti <= Ci)

  X <- as.data.frame(X)
  names(X) <- paste0("x", seq_len(p))
  res <- data.frame(time = ti, status = di, X)

  list(data = res, beta = beta, nonzero = nonzero)

}


# Debugging and sanity checking -----------------------------------------------------------------------------------

if (FALSE) {
  xdat <- sim_surv_binder(n = 100, p = 600)
}

if (FALSE) {
  n <- 1000
  p <- 100
  X <- block_corr_binder(n = n, p = p)
  Xcorr <- cor(X)

  # First couple features should have correlations around 0.5
  Xcorr[1:5, 1:5]

  # Peek at empirical correlation structure
  Xcorr <- as.data.frame(Xcorr)
  names(Xcorr) <- paste0("x", seq_len(p))
  Xcorr$var <- paste0("x", seq_len(p))

  Xcorr |>
    tidyr::pivot_longer(cols =  tidyr::starts_with("x")) |>
    dplyr::mutate(
      var = factor(var, paste0("x", seq_len(p))),
      name = factor(name, paste0("x", seq_len(p)))
    ) |>
    ggplot2::ggplot(ggplot2::aes(var, name, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_binned(type = "viridis", breaks = seq(0, 1, .25))
}

if (FALSE) {
  survdat <- sim_surv_binder(n = 200, p = 200)
}

# Original/naive/"safe" implementation
# block_corr_binder <- function(n = 50, p = 1000) {
#   stopifnot("n must be even" = n %% 2 == 0)
#
#   X <- matrix(NA_real_, nrow = n, ncol = p)
#   ui1 <- runif(n)
#   ui2 <- runif(n)
#   ui3 <- runif(n)
#
#   for (j in seq_len(ncol(X))) {
#     if (j <= 0.05 * p) {
#       X[seq_len(n/2), j] <- -1 + rnorm(n/2)
#       X[-seq_len(n/2), j] <- 1 + rnorm(n/2)
#     }
#
#     if ((j > (0.05 * p)) & (j <= (0.1 * p))) {
#       X[, j] <- 1.5 * (ui1 < 0.4) + rnorm(n)
#     }
#
#     if ((0.1 * p < j) & (j <= 0.2 * p)) {
#       X[, j] <- 0.5 * (ui2 < 0.7) + rnorm(n)
#     }
#
#     if ((0.2 * p < j) & (j <= 0.3 * p)) {
#       X[, j] <- 1.5 * (ui3 < 0.3) + rnorm(n)
#     }
#
#     if (j > 0.3 * p) {
#       X[, j] <- rnorm(n)
#     }
#   }
#
#   X <- as.data.frame(X)
#   names(X) <- paste0("x", seq_len(p))
#
#   X
# }
