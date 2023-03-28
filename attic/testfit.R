sim_surv_binder = expand.grid(n_train = 400, n_test = 200, p = 5000, ce = 0.5, lambda = 0.1, lambda_c = 0.1)


for (i in 1:5) {
  instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 5000, ce = 0.5, lambda = 0.1, lambda_c = 0.1)

  fit <- fwelnet::fwelnet_mt_cox(
    instance[["data"]],
    mt_max_iter = 2,
    z_method = "original",
    alpha = 1,
    t = 100,
    a = 0.5,
    thresh = 1e-5,
    include_mt_beta_history = TRUE
  )

  saveRDS(fit, glue::glue("attic/testfit-{i}.rds"))
}

b1 <- testfit1$beta1[, 3]
b2 <- testfit1$beta2[, 3]
hist(c(b1, b2))
hist(b1)
range(b1)

true_beta1 <- numeric(5000)
true_beta2 <- numeric(5000)
true_beta1[instance$covar_true_effect$block1] <- 0.5
true_beta2[instance$covar_true_effect$block1] <- 0.5

true_beta1[instance$covar_true_effect$block2] <- 0.5
true_beta2[instance$covar_true_effect$block2] <- -0.5

true_beta1[instance$covar_true_effect$block31] <- -0.5
true_beta2[instance$covar_true_effect$block32] <- 0.5

comptbl <- data.table::data.table(
  b1 = testfit1$beta1[, 3],
  b2 = testfit1$beta2[, 3],
  truth1 = true_beta1,
  truth2 = true_beta2
)

library(ggplot2)

comptbl |>
  dplyr::mutate(
    diff1 = truth1 - b1,
    diff2 = truth2 - b2
  ) |>
  filter(b1 > 0, b2 > 0) |>
  tidyr::pivot_longer(cols = starts_with("b"), names_to = "var") |>
  ggplot(aes(x = value)) +
  facet_grid(vars(var)) +
  geom_histogram() +
#  scale_x_log10() +
  geom_vline(xintercept = c(0, 0.5))


comptbl |>
  filter(true_beta1 != 0 | true_beta2 != 0)

comptbl |>
  filter(true_beta1 == 0 | true_beta2 == 0, b1 != 0 | b2 != 0)

comptbl |>
  summarize(
    mse1 = mean((true_beta1 - b1)^2),
    mse2 = mean((true_beta2 - b2)^2)
  )




# test ------------------------------------------------------------------------------------------------------------
source(here::here("R/20-varsel-simulation.R"))
instance <- sim_surv_binder(n_train = 400, n_test = 0, p = 5000, ce = 0.5, lambda = 0.1, lambda_c = 0.1)

fit <- tryCatch(
  fwelnet::fwelnet_mt_cox(
    instance[["data"]],
    mt_max_iter = 3,
    z_method = "original",
    alpha = 1,
    t = 100,
    a = 0.5,
    thresh = 1e-5,
    include_mt_beta_history = TRUE
  ),
  error = function(e) {
    message("Something bork")
    NULL
  }
)

beta_compare <- function(fit, instance, ce = 0.5) {
  true_beta1 <- numeric(5000)
  true_beta2 <- numeric(5000)
  true_beta1[instance$covar_true_effect$block1] <- ce
  true_beta2[instance$covar_true_effect$block1] <- ce

  true_beta1[instance$covar_true_effect$block2] <- ce
  true_beta2[instance$covar_true_effect$block2] <- -ce

  true_beta1[instance$covar_true_effect$block31] <- -ce
  true_beta2[instance$covar_true_effect$block32] <- ce

  comptbl <- data.table::data.table(
    b1 = testfit1$beta1[, fit$mt_max_iter + 1],
    b2 = testfit1$beta2[, fit$mt_max_iter + 1],
    truth1 = true_beta1,
    truth2 = true_beta2
  )
}
