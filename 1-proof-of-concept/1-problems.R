# Simulation settings -----------------------------------------------------
# true_effects <- tibble::tribble(
#   ~problem, ~beta,    ~x,  ~truth,
#   "sim_a", "beta1", "x1", 1,
#   "sim_a", "beta2", "x1", 1,
#
#   "sim_b", "beta1", "x1", 1,
#   "sim_b", "beta2", "x2", 1,
#
#   "sim_c", "beta1", "x1", 1,
#   "sim_c", "beta2", "x1", 0.25,
#
#   "sim_d", "beta1", "x1", 1,
#   "sim_d", "beta1", "x2", 0.75,
#   "sim_d", "beta1", "x3", -0.5,
#   "sim_d", "beta2", "x1", 1,
#   "sim_d", "beta2", "x2", 0.75,
#   "sim_d", "beta2", "x3", -0.5,
# )

# sim_labels <- c(
#   "sim_a" = "A: x1 has equal effect on both causes / same prevalence",
#   "sim_b" = "B: x1 has effect on cause 1, x2 has equal effect on cause 2 / same prevalence",
#   "sim_c" = "C: x1 has effect on cause 1 and smaller effect on cause2 / cause 2 less prevalent",
#   "sim_d" = "D: x{1,2,3} have equal effects in both causes / cause 2 less prevalent"
# )

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
    data = xdat,
    effects = tibble::tribble(
      ~variable, ~truth, ~cause,
      "x1", 1, 1,
      "x1", 1, 1
    )
  )
}

# B: x1 has effect on cause 1, x2 effect on cause 2, causes have same prevalence
# Censoring about as likely as event 1 or 2
sim_b <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 |
        -2 + 2*dgamma(t, 8, 2) + 1 * x2,
    n = n
  )

  list(
    data = xdat,
    effects = tibble::tribble(
      ~variable, ~truth, ~cause,
      "x1", 1, 1,
      "x2", 1, 1
    )
  )
}

# C: x1 has effect on cause 1 and cause 2, smaller effect on cause 2, cause 2 less prevalent
# Censoring about less likely than cause 1, cause 2 less likely than both
sim_c <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 |
        -4 + 2*dgamma(t, 8, 2) + 0.25 * x1,
    n = n
  )

  list(
    data = xdat,
    effects = tibble::tribble(
      ~variable, ~truth, ~cause,
      "x1", 1, 1,
      "x1", 0.25, 2
    )
  )
}

# D: Various effects but identical in both causes, cause 2 less likely
# Censoring less likely than cause 1, cause 2 less likely than both
sim_d <- function(data, job, n = 1000) {
  xdat <- sim_wrapper_cr(
    formula =
      ~ -2 + 2*dgamma(t, 8, 2) + 1 * x1 + 0.75 * x2 - 0.5 * x3 |
        -4 + 2*dgamma(t, 8, 2) + 1 * x1 + 0.75 * x2 - 0.5 * x3,
    n = n
  )

  list(
    data = xdat,
    effects = tibble::tribble(
      ~variable, ~truth, ~cause,
      "x1", 1, 1,
      "x1", 1, 2,
      "x2", 0.75, 1,
      "x2", 0.75, 2,
      "x3", -0.5, 1,
      "x3", -0.5, 1,
    )
  )
}


# Check simulation settings -----------------------------------------------

if (FALSE) {
  n <- 1000
  simres_a <- replicate(100, table(sim_a(NULL, NULL, n = n)$data$status))
  simres_b <- replicate(100, table(sim_b(NULL, NULL, n = n)$data$status))
  simres_c <- replicate(100, table(sim_c(NULL, NULL, n = n)$data$status))
  simres_d <- replicate(100, table(sim_d(NULL, NULL, n = n)$data$status))


  simres <- list(A = simres_a, B = simres_b, C = simres_c, D = simres_d) |>
    purrr::imap_dfr(~{
      rmeans <- rowMeans(.x)
      rperc <- round(100 * rmeans/1000, 2)

      xres <- as.data.frame(cbind(t(rmeans), t(rperc)))
      names(xres) <- c("Cens", "C1", "C2", "% Cens", "% C1", "% C2")
      cbind(setting = .y, xres)
    })

  simres
}
