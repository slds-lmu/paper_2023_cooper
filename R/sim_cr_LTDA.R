# Simulate competing risk data with different baseline hazards
# Vased on exercise from Lifetime Data Analysis course via AB

sim_compRisks_fullData <- function(w_param, X, beta_, max.int = 5000) {
  #browser()
  temp <- function(t_, u, x) {
    (w_param[2] * t_)^w_param[1] * exp(sum(x * beta_[, 1])) +
      1.8 * log(0.5 * t_ + 1) * exp(sum(x * beta_[, 2])) + log(1 - u)
  }

  ev.time <- function(u, x) {
    if (temp(0, u, x) * temp(max.int, u, x) < 0) {
      return(uniroot(temp, c(0, max.int), tol = 0.0001, u = u, x = x)$root)
    } else {
      return(NA)
    }
  }

  time <- apply(X, 1, function(x) ev.time(u = runif(1), x = x))

  probs <- ((w_param[1] * w_param[2]^w_param[1]) * time^(w_param[1] - 1) *

  exp(X %*% beta_[, 1])) /
  ((w_param[1] * w_param[2]^w_param[1]) * time^(w_param[1] - 1) * exp(X %*% beta_[, 1]) +
    (1.8 / (time + 2)) * exp(X %*% beta_[, 2]))

  to <- rbinom(dim(X)[1], 1, prob = probs)

  to <- ifelse(to == 1, 1, 2)

  #from <- 0
  id <- seq_along(time)

  x_df <- data.frame(X)
  names(x_df) <- paste0("x", seq_len(ncol(X)))

  data_cr <- data.frame(id, time, to, x_df)

  ## post hoc censoring
  cens_time <- rweibull(n, 0.8, 1/0.005)
  delta <- data_cr$time < cens_time
  data_cr$time <- data_cr$time * delta + cens_time * (1 - delta)
  data_cr$delta <- as.numeric(delta)
  data_cr$delta1 <- as.numeric((data_cr$to == 1) & delta)
  data_cr$delta2 <- as.numeric((data_cr$to == 2) & delta)

  data_cr
}

if (FALSE) {
  set.seed(123)

  n <- 200
  w_param <- c(3, .5)
  X <- replicate(120, rnorm(n))
  beta_ <- -cbind(
    cause1 = c(seq(0.1, 1, .1), seq(1, .1, -.1), rep(0, 100)),
    cause2 = c(seq(0.1, 1, .1), seq(1, .1, -.1), rep(0, 100))
  )

  data_cr <- sim_compRisks_fullData(w_param = w_param, X = X, beta_ = beta_)

  data_cr
}
