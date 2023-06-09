# try modelling informative censoring with cooper on bladder data?

# dry run varsel prediction with simulated data
source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))
source(here::here("2-variable-selection-sim/20-varsel-simulation.R"))

library(fwelnet)
library(survival)
library(ggplot2)
library(data.table)
library(riskRegression)

set.seed(21537)
instance <- get_bladder_data(type = "both", split = 2/3)

check_status <- function(xs) {
  tab <- as.integer(table(xs))
  mtab <- rbind(events = tab, prop = 100 * round(tab/sum(tab), 2))
  colnames(mtab) <- names(table(xs))
  mtab <- cbind(mtab, sum = rowSums(mtab))
  print(mtab)
}

check_status(instance$train$status)
check_status(instance$test$status)

instance$train$status[instance$train$status == 2] <- 1
instance$train$status <- instance$train$status + 1

instance$test$status[instance$test$status == 2] <- 1
instance$test$status <- instance$test$status + 1

check_status(instance$train$status)
check_status(instance$test$status)


mt_max_iter <- 2
tictoc::tic()
fit <- fwelnet_mt_cox(
  instance$train,
  mt_max_iter = mt_max_iter,
  stratify_by_status = TRUE,
  nfolds = 5,
  alpha = 1,
  t = 100,
  thresh = 1e-7,
  include_mt_beta_history = TRUE
)
tictoc::toc()

fit$fwfit1$glmfit$theta_store
fit$fwfit2$glmfit$theta_store

# last thetas
fit$fwfit1$glmfit$theta
fit$fwfit2$glmfit$theta

# Coefs ----
glmnet_beta1 <- fit$beta1[, 1]
glmnet_beta2 <- fit$beta2[, 1]

# Same for fwelnet estimates
fwel_beta1 <- fit$beta1[, mt_max_iter + 1]
fwel_beta2 <- fit$beta2[, mt_max_iter + 1]

fw_coefs <- list(
  fwelnet = list(
    cause1 = fwel_beta1[fwel_beta1 != 0],
    cause2 = fwel_beta2[fwel_beta2 != 0]
  ),
  glmnet = list(
    cause1 = glmnet_beta1[glmnet_beta1 != 0],
    cause2 = glmnet_beta2[glmnet_beta2 != 0]
  )
)

names(fw_coefs$fwelnet$cause1)
names(fw_coefs$glmnet$cause1)

collapse_coefs <- function(xl) {
  lapply(xl, \(x1) {
    lapply(x1, \(x2) {
      tibble::as_tibble(x2, rownames = "coef")
    }) |>
      data.table::rbindlist(idcol = "cause")
  }) |>
    data.table::rbindlist(idcol = "model")
}

coefdf <- collapse_coefs(fw_coefs)

ggplot(coefdf, aes(y = reorder(coef, value), x = value, color = model, fill = model)) +
  geom_col(position = position_dodge2(preserve = "single"), alpha = 2/3) +
  geom_vline(xintercept = 0, linetype = "solid") +
  scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill"), name = NULL) +
  labs(
    title = "Variables and selected by either method",
    subtitle = "With corresponding coefficient values",
    x = expression(beta[j]), y = "Variable"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title.position = "plot"
  )


# CSC ----
scores <- data.table::rbindlist(list(
  fit_csc_coxph(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1),
  fit_csc_coxph(instance$train, instance$test, model = "glmnet", coefs = fw_coefs$glmnet$cause1, cause = 1),

  fit_csc_coxph(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause2, cause = 2),
  fit_csc_coxph(instance$train, instance$test, model = "glmnet", coefs = fw_coefs$glmnet$cause2, cause = 2)
))


scores |>
  dplyr::filter(metric %in% c("AUC", "Brier")) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  facet_wrap(vars(metric, cause), scales = "free", labeller = labeller(cause = label_both)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = \(x) x * 100) +
  labs(
    title = "Informatice Censoring Test",
    x = "Event Times", y = "Score (%)",
    caption = "censoring = Event 1 // 1,2 -> Event 2"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title.position = "plot")
