source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("2-variable-selection-sim/20-varsel-simulation.R"))
source(here::here("5-prediction/5-prediction-wrapper.R"))
source(here::here("R/utils_riskRegression.R"))
# Ensure stuff in R/ is loaded (happens via .Rprofile)

library(fwelnet)
library(randomForestSRC)
library(survival)
library(ggplot2)
library(data.table)
library(riskRegression)

#instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 500, ce = 0.5, lambda1 = 0.1, lambda2 = 0.1, lambda_c = 0.1)
# set.seed(21537)
instance <- get_bladder_data(type = "both", split = 2/3)

check_status(instance$train$status)
check_status(instance$test$status)

# xdf <- instance$data
# xdf$status <- factor(xdf$status)
# survfit(Surv(time, status) ~ 1, data = xdf) |> plot()

# Fit models ----
# set.seed(21537) # known good for bladder resampling censoring issue
set.seed(2153)

tictoc::tic()
mt_max_iter <- 3
fit <- fwelnet::cooper(
  as.data.frame(instance$train),
  mt_max_iter = mt_max_iter,
  alpha = 1,
  t = 100,
  a = 0.5,
  stratify_by_status = TRUE,
  nfolds = 10,
  thresh = 1e-5,
  #
  standardize = TRUE
)
tictoc::toc()
# saveRDS(fit, "testfit-fwel-binder.rds")

# fit <- readRDS("testfit-fwel-binder.rds")

# Fit cs-glmnet with same lambda as cooper on its last iteration
# Proxy for initial glmnet fit, since we can't Score() that easily.
rr_glmnet_c1 <- refit_glmnet(fit, instance$test, event = 1)

# Use eval times based on test data, in 10% - 70% percentiles (arbitrarily), typ = 2 ensures times exist in data
eval_times <- quantile(instance$test$time, probs = seq(0.1, 0.7, .1), type = 2, names = FALSE)

cbfit <- CoxBoost::CoxBoost(
  time = instance$train$time,
  status = instance$train$status,
  x = as.matrix(instance$train[, -c(1, 2)]),
  cmprsk = "csh",
  stepno = 100,
  # Default is 9 * sum(status[subset] == 1), should be approx 9 * 250
  # since in our setting sum(instance$train$status != 0) is approx 233
  penalty = 1000
)

scores <- Score(
  list(cooper = fit, rr_glmnet = rr_glmnet_c1),
  formula = Hist(time, status) ~ 1,
  data = instance$test,
  cause = 1,
  metrics = c("Brier", "AUC"),
  summary = c("ibs", "ipa"),
  se.fit = FALSE,
  times = eval_times
)

scores$Brier$score
scores$AUC$score

scores_tbl <- cleanup_score(scores)


scores_tbl |>
  dplyr::filter(metric %in% c("Brier", "AUC")) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free") +
  geom_line() +
  labs(
    title = "Prediction on simulated data (Binder et al)",
    subtitle = "Cooper vs. cause-specific glmnet fit with\n lambda from cooper fit",
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


# Investigate coefficient estimates -------------------------------------------------------------------------------

coef_tbl <- data.frame(
  # Coefs from Cooper, at lambda.min
  cooper = coef(fit, event = 1),
  # Coefs from initial glmnet solution of cooper, at lambda.min
  glmnet = coef(fit, event = 1, use_initial_fit = TRUE),
  # Refit glmnet using Cooper's final lambda.min as a proxy for the initial fit
  rr_glmnet = as.vector(coef(rr_glmnet_c1$fit)),
  truth = instance$true_model$beta1
)


coef_tbl |>
  # Calc error as e.g. (truth - cooper), in-place
  dplyr::mutate(dplyr::across(cooper:rr_glmnet, \(x) truth - x)) |>
  tidyr::pivot_longer(cols = cooper:rr_glmnet, names_to = "model", values_to = "error") |>
  dplyr::filter(truth != 0) |>
  ggplot(aes(x = error, color = model, fill = model)) +
  geom_density(alpha = 1/3)
  geom_histogram(position = "dodge")


# Cooper's false positives
coef_tbl |>
  dplyr::filter(truth == 0 & cooper != 0) |>
  round(3)

# True effect features
coef_tbl |>
  dplyr::filter(truth != 0) |>
  round(3)

coef_z <- coef_tbl |>
  dplyr::filter(truth == 0)

lapply(coef_z, summary)

coef_tbl |>
  # Calc error as e.g. (truth - cooper), in-place
  tidyr::pivot_longer(cols = cooper:rr_glmnet, names_to = "model", values_to = "estimate") |>
  dplyr::mutate(
    tp = truth != 0 & estimate != 0,
    fp = truth == 0 & estimate != 0,
    tn = truth == 0 & estimate == 0,
    fn = truth != 0 & estimate == 0
  ) |>
  dplyr::group_by(model) |>
  dplyr::summarize(dplyr::across(tp:fn, sum)) |>
  dplyr::mutate(
    tpr = tp/sum(instance$true_model$beta1 != 0),
    fpr = fp/sum(instance$true_model$beta1 == 0),
    ppv = tp/(tp + fp),
    npv = tn/(tn + fn),
    fdr = 1 - ppv,
    f1  = (2 * ppv * tpr) / (ppv + tpr),
    acc = (tp+tn)/length(instance$true_model$beta1)
  )


# Qintuple check that the coefs are where and what I expect them to be.
lambda_min_id <- which(fit$fwelfits[[1]]$lambda.min == fit$fwelfits[[1]]$lambda)
all.equal(
  # Fwelnet CV best beta at given lambda
  fit$fwelfits[[1]]$glmfit$beta[, lambda_min_id],
  # internal glmnet best beta at that same lambda, should be identical (if standardize = FALSE)
  unname(as.matrix(fit$fwelfits[[1]]$glmfit$glmfit$beta))[, lambda_min_id]
)


# Investigate linear predictors? ----------------------------------------------------------------------------------

linpred_tbl <- do.call(cbind, lapply(coef_tbl, \(beta) as.vector(fit$x %*% beta))) |>
  as.data.frame()

linpred_tbl |>
  # Calc error as e.g. (truth - cooper), in-place
  dplyr::mutate(dplyr::across(cooper:rr_glmnet, \(x) truth - x)) |>
  tidyr::pivot_longer(cols = cooper:rr_glmnet, names_to = "model", values_to = "error") |>
  dplyr::filter(truth != 0) |>
  ggplot(aes(x = error, color = model, fill = model)) +
  geom_density(alpha = 1/3) +
  theme_minimal() +
  theme(legend.position = "bottom")


# Same but with standardized coefs? -------------------------------------------------------------------------------
source(here::here("2-variable-selection-sim/20-varsel-simulation.R"))

library(fwelnet)
library(randomForestSRC)
library(survival)
library(ggplot2)
library(data.table)
library(riskRegression)


set.seed(2153)
instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 1000, ce = 0.5, lambda1 = 0.1, lambda2 = 0.1, lambda_c = 0.1)

set.seed(2153)
tictoc::tic()
mt_max_iter <- 3
cooperfit <- fwelnet::fwelnet_mt_cox(
  as.data.frame(instance$train), mt_max_iter = mt_max_iter,
  alpha = 1, t = 100, a = 0.5, stratify_by_status = TRUE, nfolds = 10, thresh = 1e-5,
  standardize = TRUE
)
tictoc::toc()

cooperfit_standardized_bak <- cooperfit

set.seed(2153)
tictoc::tic()
cooperfit_nonst <- fwelnet::fwelnet_mt_cox(
  as.data.frame(instance$train), mt_max_iter = mt_max_iter,
  alpha = 1, t = 100, a = 0.5, stratify_by_status = TRUE, nfolds = 10, thresh = 1e-5,
  standardize = FALSE
)
tictoc::toc()

xn <- dimnames(cooperfit$fwelfits[[1]]$glmfit$glmfit$beta)
xn2 <- dimnames(cooperfit$fwelfits[[2]]$glmfit$glmfit$beta)

# cooperfit$fwelfits[[1]]$glmfit$glmfit$beta <- Matrix::Matrix(cooperfit$fwelfits[[1]]$glmfit$beta, dimnames = xn)
# cooperfit$fwelfits[[2]]$glmfit$glmfit$beta <- Matrix::Matrix(cooperfit$fwelfits[[2]]$glmfit$beta, dimnames = xn2)

identical(
  cooperfit$fwelfits[[1]]$glmfit$glmfit$beta,
  Matrix::Matrix(cooperfit$fwelfits[[1]]$glmfit$beta, dimnames = xn)
)

# Use eval times based on test data, in 10% - 70% percentiles (arbitrarily), typ = 2 ensures times exist in data
eval_times <- quantile(instance$test$time, probs = seq(0.1, 0.7, .1), type = 2, names = FALSE)


rr_glmnet_c1 <- refit_glmnet(cooperfit_nonst, instance$test, event = 1)

scores <- Score(
  list(
    cooper = cooperfit,
    #cooper_nonst = cooperfit_nonst,
    glmnet = rr_glmnet_c1),
  formula = Hist(time, status) ~ 1,
  data = instance$test,
  cause = 1,
  metrics = c("Brier", "AUC"),
  summary = c("ibs", "ipa"),
  se.fit = FALSE,
  times = eval_times
)

cleanup_score(scores) |>
  #dplyr::filter(metric %in% c("Brier", "AUC")) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  #facet_grid(cols = vars(model), rows = vars(metric), scales = "free") +
  facet_grid(rows = vars(metric), scales = "free") +
  geom_line() +
  labs(
    title = "Prediction on simulated data (Binder et al)",
    subtitle = "Cooper vs. cause-specific glmnet fit with\n lambda from cooper fit",
    x = "Event Time (t)", y = "Metric",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

