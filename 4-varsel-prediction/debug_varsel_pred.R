# dry run varsel prediction with simulated data
source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))
source(here::here("2-variable-selection-sim/20-varsel-algorithms.R"))

library(fwelnet)
library(randomForestSRC)
library(survival)
library(ggplot2)
library(data.table)
library(riskRegression)

# instance <- sim_surv_binder(n_train = 600, n_test = 400, p = 5000, ce = 0.5, lambda1 = 0.1, lambda2 = 0.1, lambda_c = 0.1)
# table(instance$train$status) |> (\(x) {print(x); cat("\n"); x})() |> as.numeric() |> summary()
set.seed(21537)
instance <- get_bladder_data(type = "both", split = 2/3)
table(instance$train$status) |> (\(x) {print(x); cat("\n"); x})() |> as.numeric() |> summary()
table(instance$test$status) |> (\(x) {print(x); cat("\n"); x})() |> as.numeric() |> summary()

min(table(instance$test$status)) / 10

check_status <- function(xs) {
  tab <- as.integer(table(xs))
  mtab <- rbind(events = tab, prop = 100 * round(tab/sum(tab), 2))
  colnames(mtab) <- unique(xs)
  mtab <- cbind(mtab, sum = rowSums(mtab))
  print(mtab)
}

check_status(instance$train$status)
check_status(instance$test$status)

# xdf <- instance$data
# xdf$status <- factor(xdf$status)
# survfit(Surv(time, status) ~ 1, data = xdf) |> plot()

# from fwel_mt_varselect_pred ----
set.seed(21537) # known good
set.seed(2153)

tictoc::tic()
mt_max_iter <- 2
fit <- fwelnet::fwelnet_mt_cox(
  as.data.frame(instance$train),
  mt_max_iter = mt_max_iter,
  alpha = 1,
  t = 100,
  a = 0.5,
  stratify_by_status = TRUE,
  nfolds = 10,
  thresh = 1e-5,
  include_mt_beta_history = TRUE
)

tictoc::toc()

tfolds <- fwelnet::stratified_cv_folds(instance$train, 10)
table(tfolds$fold, tfolds$status)

fit$fwfit1$glmfit$theta_store
fit$fwfit2$glmfit$theta_store

# last thetas
fit$fwfit1$glmfit$theta
fit$fwfit2$glmfit$theta

saveRDS(object = fit, file = here::here("fwel-testfit-binder.rds"))
#saveRDS(object = fit, file = here::here("fwel-testfit-bladder.rds"))


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

names(fw_coefs$fwelnet$cause1)
names(fw_coefs$glmnet$cause1)

setdiff(names(fw_coefs$fwelnet$cause1), names(fw_coefs$glmnet$cause1))

if (setequal(names(fw_coefs$fwelnet$cause1), names(fw_coefs$glmnet$cause1))) {
  print(names(fw_coefs$fwelnet$cause1))
  stop("fwelnet and glmnet have selected identical variables")
}

fit_csc(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1)

fit_csc_coxph(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1)

# from fit_csc

train <- as.data.table(instance$train)
test <- as.data.table(instance$test)

# for fwel cause 1 ----
coefs <- fw_coefs$fwelnet$cause1
cause <- 1

# sub-select data
train_cause <- train[, c("time", "status", names(coefs)), with = FALSE]
test_cause <- test[, c("time", "status", names(coefs)), with = FALSE]

csc_model_fw <- CSC(formula = Hist(time, status) ~ ., data = train_cause, cause = cause)

cbind(cox = csc_model_fw$models$`Cause 1`$coefficients, fwelnet = coefs) |>
  round(3)

# lp
lp <- as.numeric(as.matrix(train_cause[, -c(1, 2)]) %*% coefs)
range(lp)
bfit <- survival::basehaz(csc_model_fw$models$`Cause 1`)
quantile(bfit$time, probs = seq(0.1, 0.5, .1))

bfit |>
  dplyr::filter(time == quantile(time, prob = .5, type = 1))

plot(bfit$time, bfit$hazard)

# surv prob
#exp(-0.1075885 * exp(lp)) |> summary()
exp(-0.66371930 * exp(lp)) |> summary()

# Get quantiles of time points from full dataset to allow later aggregation per timepoint
round(range(c(train$time, test$time)), 2)
eval_times <- quantile(c(train$time, test$time), probs = seq(0.1, 0.6, .1), names = FALSE, type = 1)
#eval_times <- quantile(test$time, probs = seq(0.1, 0.6, .1), names = FALSE, type = 1)

# Yikes
# Estimated risk outside the range [0,1].
# Possible cause: incorrect extrapolation, i.e., time and/or covariates used for the prediction differ from those used to fit the Cox models.
model = "fwelnet"
mod_scores_auc <- Score(
  list(fwelnet = csc_model_fw),
  predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = "AUC",
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)


mod_scores_brier <- Score(
  list(fwelnet = csc_model_fw),
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = c("Brier"),
  summary = c("ibs", "ipa"),
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)

off_preds <- which(mod_scores_brier$off.predictions$negative.values > 0 | mod_scores_brier$off.predictions$values.above.1 > 0)
off_times <- mod_scores_brier[off_preds]$times
summary(test$time)
test[time %in% off_times, -c(1, 2)]
test[time %in% off_times, c(1, 2)]


lpoff <- as.matrix(test[time %in% off_times, -c(1, 2)]) %*% instance$true_model$beta1

exp(-median(test$time) * exp(lpoff)) |> summary()


####
train_cause_1 <- train_cause
train_cause_1$status[train_cause_1$status == 2] <- 0

test_cause <- test_cause
test_cause$status[test_cause$status == 2] <- 0

cph_fw <- coxph(Surv(time, status) ~ ., data = train_cause_1, x = TRUE)

mod_scores_brier_cph <- Score(
  list(fwelnet = cph_fw),
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Surv(time, status) ~ 1,
  data = test_cause,
  metrics = c("Brier"),
  summary = c("ibs", "ipa"),
  se.fit = FALSE,
  times = eval_times
)

mod_scores_brier_cph

all.equal(coef(cph_fw), coef(csc_model_fw$models$`Cause 1`))
####


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

csc_model_glmn <- CSC(formula = Hist(time, status) ~ ., data = train_cause, cause = cause)


cbind(cox = csc_model_glmn$models$`Cause 1`$coefficients, glmnet = coefs) |>
  round(3)

# lp
lp <- as.numeric(as.matrix(train_cause[, -c(1, 2)]) %*% coefs)
range(lp)
bfit <- survival::basehaz(csc_model_glmn$models$`Cause 1`)
quantile(bfit$time, probs = seq(0.1, 0.5, .1))

bfit |>
  dplyr::filter(time == median(time))

plot(bfit$time, bfit$hazard)

# surv prob
#exp(-0.1075885 * exp(lp)) |> summary()
exp(-0.087395 * exp(lp)) |> summary()

mod_scores_auc <- Score(
  list(glmnet = csc_model_glmn),
  #predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = test_cause,
  metrics = "AUC",
  cause = cause,
  se.fit = FALSE,
  times = eval_times
)


mod_scores_brier <- Score(
  list(glmnet = csc_model_glmn),
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

result <- rbind(result_fwelnet, result_glmnet)
result <- result[!(is.na(score) & model == "Null model"), ]
result <- result[!(metric == "IPA" & model == "Null model"), ]

result |>
  dplyr::filter(score < 2 & score > -1) |>
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

# compare different impls

res_csc <- fit_csc(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1)
res_cph <- fit_csc_coxph(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1)
res_cph_glm <- fit_csc_coxph(instance$train, instance$test, model = "glmnet", coefs = fw_coefs$glmnet$cause1, cause = 1)

res_csc |>
  dplyr::filter(score >= 0 & score <= 1) |>
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

res_cph |>
  rbind(res_cph_glm) |>
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
