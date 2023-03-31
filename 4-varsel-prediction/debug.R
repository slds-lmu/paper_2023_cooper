source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))

instance <- get_bladder_data(split = 2/3, standardize = TRUE)

quantile(instance$train$time, probs = seq(0.1, 0.9, .1), names = FALSE)
quantile(instance$test$time, probs = seq(0.1, 0.9, .1), names = FALSE)

quantile(c(instance$train$time, instance$test$time), probs = seq(0.1, 0.9, .1), names = FALSE)

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


data.table::rbindlist(list(
  fwe$scores,
  rf$scores,
  cb$scores
), fill = TRUE) |>
  ggplot(aes(x = times, y = score, color = model)) +
  facet_wrap(vars(metric, cause)) +
  geom_path() +
  theme_minimal()



res <- data.table::rbindlist(lapply(1:1000, function(i) {
  instance <- get_bladder_data(split = 2/3, standardize = TRUE)
  train_props <- as.data.frame(instance$train_event_prop)
  names(train_props) <- c("cause", "train_prop")
  test_props <- as.data.frame(instance$test_event_prop)
  names(test_props) <- c("cause", "test_prop")
  res <- merge(train_props, test_props, by = "cause")
  res$iter <- i
  res
}))

res |>
  mutate(diff = train_prop - test_prop) |>
  distinct(diff, .keep_all = TRUE)
  pull(diff) |>
  unique() |> length()



res <- data.table::rbindlist(lapply(1:1000, function(iter) {
  instance <- get_bladder_data(split = 2/3, standardize = TRUE)
  setequal(setdiff(instance$train$time, instance$test$time), instance$train$time)
  setequal(setdiff(instance$train_row_ids, instance$test_row_ids), instance$train_row_ids)

  instance$train_row_ids
  instance$test_row_ids

  data.frame(
    row_ids = c(sort(instance$train_row_ids), instance$test_row_ids),
    set = c(rep("train", length(instance$train_row_ids)), rep("test", length(instance$test_row_ids))),
    iter = iter
  )
}))

res |>
  count(row_ids, set) |>
  tidyr::pivot_wider(id_cols = "row_ids", names_from = "set", values_from = "n") |>
  mutate(ratio = test/train) |>
  ggplot(aes(x = ratio)) +
  geom_density()
