# Boosting test maybe
renv::install("binderh/CoxBoost")
library(CoxBoost)
library(ggplot2)

instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 400, ce = 0.5, lambda = 0.1, lambda_c = 0.1)

# 9 * sum(instance$data$status != 0)

cbfit <- CoxBoost(
  time = instance$data$time,
  status = instance$data$status,
  x = as.matrix(instance$data[, -c(1, 2)]),
  cmprsk = "sh",
  stepno = 100,
  # Default is 9 * sum(status[subset] == 1), should be approx 9 * 250
  # since in our setting sum(instance$data$status != 0) is approx 233
  penalty = 2000
)

cbfit$xnames
cbfit$model$`1`$linear.predictor
cbfit$model$`1`$coefficients |> str()

cbfit$model$`1`$coefficients[101,]

coef(cbfit) |> str()


test <- cv.CoxBoost(
  time = instance$data$time,
  status = instance$data$status,
  x = as.matrix(instance$data[, -c(1, 2)]),
  cmprsk = "ccsh",
  maxstepno = 1000
)
