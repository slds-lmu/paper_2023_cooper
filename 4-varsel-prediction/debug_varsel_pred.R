# dry run varsel prediction with simulated data
source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))
source(here::here("R/20-varsel-simulation.R"))

library(fwelnet)
library(randomForestSRC)
library(survival)
library(rlang)

# instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 5000, ce = 0.25, lambda1 = 0.1, lambda2 = 0.1, lambda_c = 0.1)
# table(instance$train$status) |> (\(x) {print(x); cat("\n"); x})() |> as.numeric() |> summary()

instance <- get_bladder_data(type = "geno", split = 2/3)
table(instance$train$status) |> (\(x) {print(x); cat("\n"); x})() |> as.numeric() |> summary()

min(table(instance$train$status)) / 20

# xdf <- instance$data
# xdf$status <- factor(xdf$status)
# survfit(Surv(time, status) ~ 1, data = xdf) |> plot()

# from fwel_mt_varselect_pred ----
tictoc::tic()
mt_max_iter <- 3
fit <- fwelnet::fwelnet_mt_cox(
  as.data.frame(instance$train),
  mt_max_iter = mt_max_iter,
  alpha = 1,
  t = 100,
  a = 0.5,
  stratify_by_status = TRUE,
  nfolds = 5,
  thresh = 1e-3,
  include_mt_beta_history = TRUE
)

tictoc::toc()

#saveRDS(object = fit, file = here::here("fwel-testfit.rds"))
saveRDS(object = fit, file = here::here("fwel-testfit-bladder.rds"))


p <- ncol(instance$train)

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

setdiff(names(fw_coefs$fwelnet$cause1), names(fw_coefs$glmnet$cause1))
setequal(names(fw_coefs$fwelnet$cause1), names(fw_coefs$glmnet$cause1))

if (setequal(names(fw_coefs$fwelnet$cause1), names(fw_coefs$glmnet$cause1))) {
  print(names(fw_coefs$fwelnet$cause1))
  stop("fwelnet and glmnet have selected identical variables")
}


# from fit_csc
library(data.table)
library(riskRegression)
train <- as.data.table(instance$train)
test <- as.data.table(instance$test)

# for fwel cause 1 ----
coefs <- fw_coefs$fwelnet$cause1
cause <- 1

# sub-select data
train_cause <- train[, c("time", "status", names(coefs)), with = FALSE]
test_cause <- test[, c("time", "status", names(coefs)), with = FALSE]

csc_model <- CSC(formula = Hist(time, status) ~ ., data = train_cause, cause = cause)


cbind(cox = csc_model$models$`Cause 1`$coefficients, fwelnet = coefs) |>
  round(3 )

# lp
lp <- as.numeric(as.matrix(train_cause[, -c(1, 2)]) %*% coefs)
range(lp)
bfit <- survival::basehaz(csc_model$models$`Cause 1`)
quantile(bfit$time, probs = seq(0.1, 0.5, .1))

bfit |>
  dplyr::filter(time == median(time))

plot(bfit$time, bfit$hazard)

# surv prob
#exp(-0.1075885 * exp(lp)) |> summary()
exp(-0.087395 * exp(lp)) |> summary()

# Get quantiles of time points from full dataset to allow later aggregation per timepoint
eval_times <- quantile(c(train$time, test$time), probs = seq(0.1, 0.75, .05), names = FALSE)
#eval_times <- quantile(test$time, probs = seq(0.1, 0.75, .1), names = FALSE)

# Yikes
# Estimated risk outside the range [0,1].
# Possible cause: incorrect extrapolation, i.e., time and/or covariates used for the prediction differ from those used to fit the Cox models.
model = "fwelnet"
mod_scores_auc <- Score(
  list2(!!model := csc_model), # uses rlang splicing to get list(fwelnet = csc_1)
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = "AUC",
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)


mod_scores_brier <- Score(
  list2(!!model := csc_model),
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = c("Brier"),
  summary = c("ibs", "ipa"),
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)


result_fwelnet <- mod_scores_auc$AUC$score[mod_scores_brier$Brier$score, on = c("model", "times")]
result_fwelnet <- data.table::melt(result_fwelnet, id.vars = c("model", "times"),
                           value.name = "score", variable.name = "metric",
                           meausure.vars = c("AUC", "Brier", "IBS", "IPA"))

result_fwelnet$cause <- cause
result_fwelnet

# for glmnet cause 1 ----
coefs <- fw_coefs$glmnet$cause1

# sub-select data
train_cause <- train[, c("time", "status", names(coefs)), with = FALSE]
test_cause <- test[, c("time", "status", names(coefs)), with = FALSE]

csc_model <- CSC(formula = Hist(time, status) ~ ., data = train_cause, cause = cause)


cbind(cox = csc_model$models$`Cause 1`$coefficients, glmnet = coefs) |>
  round(3 )

# lp
lp <- as.numeric(as.matrix(train_cause[, -c(1, 2)]) %*% coefs)
range(lp)
bfit <- survival::basehaz(csc_model$models$`Cause 1`)
quantile(bfit$time, probs = seq(0.1, 0.5, .1))

bfit |>
  dplyr::filter(time == median(time))

plot(bfit$time, bfit$hazard)

# surv prob
#exp(-0.1075885 * exp(lp)) |> summary()
exp(-0.087395 * exp(lp)) |> summary()

model = "glmnet"
mod_scores_auc <- Score(
  list2(!!model := csc_model), # uses rlang splicing to get list(fwelnet = csc_1)
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = "AUC",
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)


mod_scores_brier <- Score(
  list2(!!model := csc_model),
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = c("Brier"),
  summary = c("ibs", "ipa"),
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)


result_glmnet <- mod_scores_auc$AUC$score[mod_scores_brier$Brier$score, on = c("model", "times")]
result_glmnet <- data.table::melt(result_glmnet, id.vars = c("model", "times"),
                           value.name = "score", variable.name = "metric",
                           meausure.vars = c("AUC", "Brier", "IBS", "IPA"))
result_glmnet$cause <- cause

# Aggregate + plot ----
library(ggplot2)
result <- rbind(result_fwelnet, result_glmnet)

result |>
  #dplyr::filter(metric == "Brier") |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free") +
  geom_line() +
  labs(
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
