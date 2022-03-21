source(here::here("R/sim_cr.R"))

# Simulation settings -----------------------------------------------------
true_effects <- tibble::tribble(
  ~problem, ~beta,    ~x,  ~truth,
  "sim_a", "beta1", "x1", 1,
  "sim_a", "beta2", "x1", 1,

  "sim_b", "beta1", "x1", 1,
  "sim_b", "beta2", "x2", 1,

  "sim_c", "beta1", "x1", 1,
  "sim_c", "beta2", "x2", 0.25,

  "sim_d", "beta1", "x1", 1,
  "sim_d", "beta1", "x2", 0.75,
  "sim_d", "beta1", "x3", -0.5,
  "sim_d", "beta2", "x1", 1,
  "sim_d", "beta2", "x2", 0.75,
  "sim_d", "beta2", "x3", -0.5,
)

# A: x1 has the same (large) effect in both causes, both causes have equal prevalence
# Censoring about as likely as event 1 or 2
sim_a <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 |
        -2 + 2*dgamma(t, 8, 2) + 1 * x1,
    n = n
  )

  list(
    data = xdat
  )
}

# B: x1 has effect on cause 1, x2 effect auf cause 2, causes have same prevalence
# Censoring about as likely as event 1 or 2
sim_b <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 |
        -2 + 2*dgamma(t, 8, 2) + 1 * x2,
    n = n
  )

  list(
    data = xdat
  )
}

# C: x1 has effect on cause 1 and cause 2, smaller effect on cause 2, cause 2 less prevalent
# Censoring about less likely than cause 1, cause 2 less likely than both
sim_c <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 |
        -4 + 2*dgamma(t, 8, 2) + 0.25 * x2,
    n = n
  )

  list(
    data = xdat
  )
}

# D: Various effects but identical in both causes, cause 2 less likely
# Censoring about less likely than cause 1, cause 2 less likely than both
sim_d <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 + 0.75 * x2 - 0.5 * x3 |
        -4 + 2*dgamma(t, 8, 2) + 1 * x1 + 0.75 * x2 - 0.5 * x3,
    n = n
  )

  list(
    data = xdat
  )
}
