# Simulation wrapper functions for high dimension data simulation
# Batchtools main script: 2-run-batchtools.R


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
#'
#' Block 4 is not filled with correlated data to adhere to setup in Binder (2009), p. 893
#' @param n Number of obs. (must be even).
#' @param p Number of predictors.
block_corr_binder <- function(n = 400, p = 5000) {
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

  # Binder (2009) only use 3 blocks of correlated variables, we still include it as part of the setup though.
  block4 <- which((0.2 * p < j_seq) & (j_seq <= 0.3 * p))
  X[, block4] <- 1.5 * (ui3 < 0.3) + X[, block4]

  X
}

#' Take a covariate matrix and simulate a survival outcome, based on the same paper.
#' @param job,data For batchtools use only.
#' @param n_train,n_test Number of observations in train- and test sets. Training data will be in element `$train`,
#'   and test data in element `$test`.
#' @param p Passed to block_corr_binder().
#' @param ce Effect constant. Used as magnitude of effect. Binder et al 2009 use 0.5.
#' @param lambda,lambda_c Baseline-hazard constants for Cox-exponential model for survival & censoring times.
#' @param lambda1,lambda2 Baseline-hazard constant for cause 1 and 2 respectively. By default, both will be set to
#'   the value of `lambda`, resulting in roughly equal event numbers for each cause. If set, these values supersede
#'   `lambda`.
#'
#' Binder 2008 use ce in {0.05, 0.075, 0.1} for the censored survival setting.
#' Target Binder 2009 paper uses ce = +/-0.5.
#' @note Binder et. al. describe U = runif(n), but log(runif(n)) is used to ensure results as stated
#' with mean baseline survival time of around 10, rather than negative survival times.
#' This is also consistent with other descriptions of the Cox exponential model.
sim_surv_binder <- function(job, data,
                            n_train = 50, n_test = 0,
                            p = 1000,
                            ce = 0.5,
                            lambda = 0.1,
                            lambda1 = lambda, lambda2 = lambda,
                            lambda_c = 0.1) {
  checkmate::assert_integerish(n_train, lower = 10, len = 1)
  checkmate::assert_integerish(n_test, lower = 0, len = 1)
  checkmate::assert_true((n_train + n_test) %% 2 == 0)
  checkmate::assert_integerish(p, lower = 400, len = 1)
  checkmate::assert_numeric(ce, lower = 0.01, len = 1)
  checkmate::assert_numeric(lambda, lower = 0.01, len = 1)
  checkmate::assert_numeric(lambda1, lower = 0.01, len = 1)
  checkmate::assert_numeric(lambda2, lower = 0.01, len = 1)
  checkmate::assert_numeric(lambda_c, lower = 0.01, len = 1)

  n <- n_train + n_test

  X <- block_corr_binder(n = n, p = p)
  beta1 <- numeric(p)
  beta2 <- beta1

  # Original Binder 2008 setup differs from target setup in Binder 2009
  # j * 200/p: if odd -> ce, if even -> -ce, else 0
  # tmp <- ((seq_len(p) * 200) / p) %% 2
  # j_odd <- which(tmp == 1)
  # j_even <- which(tmp == 0)
  # beta[j_odd] <- ce
  # beta[j_even] <- -ce

  # From Binder 2009: (increasing/decreasing means ce = +/- 0.5 )
  # In the first block, four covariates have an increasing effect on both hazards.
  # The second block has four informative covariates with an
  #   increasing effect on the cause-specific hazard for type 1 events, and a
  #   decreasing effect on the competing cause-specific hazard.
  # In the third block, there are
  #   - four covariates that have a decreasing effect on the event type 1 hazard only, and
  #   - four other covariates that have an increasing effect on the event type 2 hazard.

  j_seq <- seq_len(p)

  # first block: just use coefs 1-4, since block starts at j = 1
  j_block1 <- which(j_seq <= 0.05 * p)
  beta1[j_block1[1:4]] <- ce
  beta2[j_block1[1:4]] <- ce

  # second block:
  j_block2 <- which((j_seq > (0.05 * p)) & (j_seq <= (0.1 * p)))
  beta1[j_block2[1:4]] <- ce
  beta2[j_block2[1:4]] <- -ce

  # third block
  j_block3 <- which((0.1 * p < j_seq) & (j_seq <= 0.2 * p))
  beta1[j_block3[1:4]] <- -ce
  beta2[j_block3[5:8]] <- ce # offset by 4, b/c "four *other* covariates..."

  # Unused for true effects but tracked to see if false positives pop up here
  j_block4 <- which((0.2 * p < j_seq) & (j_seq <= 0.3 * p))

  # Noise variables, uncorrelated, also no true effects
  j_noise <- which(j_seq > 0.3 * p)

  # Save indices of effect variables grouped by their covar blocks for later
  # since effects/direction on causes depends on blocks, and they have different correlations
  covar_true_effect <- list(
    block1 = j_block1[1:4],
    block2 = j_block2[1:4],
    block31 = j_block3[1:4],
    block32 = j_block3[5:8],
    block4 = integer(0), # no true effects, so empty set (of integers for later assertions)
    noise = integer(0)
  )

  # Blocks have unequal length so don't use a data.frame to avoid recycling
  covar_blocks <- list(
    block1 = j_block1,
    block2 = j_block2,
    block31 = j_block3,
    block32 = j_block3,
    block4 = j_block4,
    noise = j_noise
  )

  # Keep track of number of nonzero effects, total count is identical for both causes in setup above
  nonzero1 <- sum(beta1 > 0)
  nonzero2 <- sum(beta2 > 0)
  checkmate::assert_count(nonzero1, positive = TRUE, null.ok = FALSE)
  checkmate::assert_count(nonzero2, positive = TRUE, null.ok = FALSE)
  # Triple sanity check we have 16 unique informative covariates
  # checkmate::assert_true(length(c(1:4, j_block2[1:4], j_block3[1:4], j_block3[5:8])) == 16)
  checkmate::assert_true(
    length(c(
      covar_true_effect$block1,
      covar_true_effect$block2,
      covar_true_effect$block31,
      covar_true_effect$block32)) == 16
  )


  # Linear predictors
  lp1 <- X %*% beta1
  lp2 <- X %*% beta2

  # Survival and censoring times
  Ti1 <- -log(runif(n)) / (lambda1 * exp(lp1))
  Ti2 <- -log(runif(n)) / (lambda2 * exp(lp2))
  # Default lambda_c = 0.1 should yield roughly 36% censored events.
  Ci <- -log(runif(n)) / lambda_c

  ti <- pmin(Ti1, Ti2, Ci)
  # Censored where neither event occurs before censoring time
  # if event before censoring time -> TRUE -> integer 1
  di <- as.integer(Ti1 <= Ci | Ti2 <= Ci)
  # Set status == 2 if Ti2 is observed
  di[which(Ti2 <= Ti1 & Ti2 <= Ci)] <- 2

  # Glimpse for debugging because thinking is hard
  # data.frame(Ti1, Ti2, Ci, ti, di)

  X <- as.data.frame(X)
  names(X) <- paste0("x", seq_len(p))
  res <- data.frame(time = ti, status = di, X)

  if (n_test > 0) {
    id_train <- sample.int(n, size = n_train)
  } else {
    id_train = seq_len(n_train)
  }

  list(
    train = res[id_train, ],
    test = res[-c(id_train), ],
    covar_blocks = covar_blocks,
    covar_true_effect = covar_true_effect,
    true_model = list(
      beta1 = beta1, beta2 = beta2,
      lp1 = lp1, lp2 = lp2
    )
  )
}

# Debugging and sanity checking -----------------------------------------------------------------------------------

if (FALSE) {
  # Check event times match expectations
  library(dplyr)
  n <- 400

  xdat <- sim_surv_binder(n_train = n, p = 5000, lambda1 = 0.01, lambda2 = 0.1, lambda_c = 0.1)
  table(xdat$train$status)

  xsum <- purrr::map_df(1:100, ~{
    status <- sim_surv_binder(n_train = n, p = 5000, lambda1 = 0.1, lambda2 = 0.1)$train$status
    data.frame(rep = .x, table(status))
  })

  xsum |>
    group_by(status) |>
    summarize(
      freq_min = min(Freq),
      freq_mean = mean(Freq),
      freq_max = max(Freq),
      prop_min = min(Freq/n),
      prop_mean = mean(Freq/n),
      prop_max = max(Freq/n),
      .groups = "keep"
    ) |>
    mutate(across(starts_with("prop"), \(x) scales::label_percent(accuracy = .1)(x))) |>
    transmute(
      n = glue::glue("{freq_mean} ({freq_min} - {freq_max})"),
      prop = glue::glue("{prop_mean} ({prop_min} - {prop_max})")
    )
}

# Original/naive/"safe" implementation for reference
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
