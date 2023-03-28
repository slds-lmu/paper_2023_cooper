source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))

instance <- get_bladder_data(split = 2/3, standardize = TRUE)

library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
library(rlang)
library(ggplot2)

tictoc::tic()
fwe <- fwel_mt_varselect_pred(data = NULL, job = NULL, instance = instance)
# Error: from glmnet C++ code (error code 30000); Inititialization numerical error; probably too many censored observations
tictoc::toc()
tictoc::tic()
rf <- rfsrc_varselect_pred(data = NULL, job = NULL, instance = instance)
tictoc::toc()
tictoc::tic()
cb <- coxboost_varselect_pred(data = NULL, job = NULL, instance = instance)
tictoc::toc()

str(fwe)
str(rf)
str(cb)



# ret_glm <- fit_csc(instance$train, instance$test, model = "glmnet", coef_cause1 = fwe$glmnet$cause1, coef_cause2 = fwe$glmnet$cause2)
#
# ret_rf <- fit_csc(instance$train, instance$test, model = "rfsrc", coef_cause1 = rf$cause1, coef_cause2 = rf$cause2)
# ret_cb <- fit_csc(instance$train, instance$test, model = "coxboost", coef_cause1 = cb$cause1, coef_cause2 = cb$cause2)

data.table::rbindlist(list(
  fwe$scores,
  rf$scores,
  cb$scores
), fill = TRUE) |>
  ggplot(aes(x = times, y = score, color = model)) +
  facet_wrap(vars(metric)) +
  geom_path() +
  theme_minimal()

