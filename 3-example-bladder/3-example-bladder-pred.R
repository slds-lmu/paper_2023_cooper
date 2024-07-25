invisible(lapply(list.files("R", pattern = "*.R", full.names = TRUE), source, echo = FALSE))
library(data.table)
library(cooper)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
if (!dir.exists(here::here("results"))) dir.create(here::here("results"))
  set.seed(2023)

bladder_file <- here::here("data/bladder-binder-clinical_geno.rds")
if (!file.exists(bladder_file)) {
  message("Recreating bladder dataset from data-raw/preprocess-bladder-data.R")
  source(here::here("data-raw/preprocess-bladder-data.R"))
}
bladder <- readRDS(bladder_file)

# Create list with $train and $test data
splits <- partition_dt(bladder, train_prop = 0.7)
saveRDS(splits, here::here("data/bladder-split.rds"))

# CooPeR --------------------------------------------------------------------------------------
cli::cli_alert_info("Fitting CooPeR")

cooperfit <- cooper::cooper(
  splits$train, mt_max_iter = 3,
  alpha = 1, t = 100, thresh = 1e-7,
  stratify_by_status = TRUE, nfolds = 5
)

cli::cli_alert_success("Saving CooPeR model")

saveRDS(cooperfit, here::here("results/3-bladder-cooper.rds"))

# rfsrc ---------------------------------------------------------------------------------------
cli::cli_alert_info("Fitting rfsrc")

rf_c1 <- rfsrc_tuned(
  xdat = splits$train,
  splitrule = "logrank",
  importance = "random",
  cause = 1
)

rf_c2 <- rfsrc_tuned(
  xdat = splits$train,
  splitrule = "logrank",
  importance = "random",
  cause = 2
)

cli::cli_alert_success("Saving rfsrc models")

saveRDS(rf_c1, here::here("results/3-bladder-rfsrc-c1.rds"))
saveRDS(rf_c2, here::here("results/3-bladder-rfsrc-c2.rds"))

# Coxboost ----------------------------------------------------------------
cli::cli_alert_info("Fitting CoxBoost")

cbfit <- coxboost_tuned(
  xdat = splits$train,
  cmprsk = "csh",
  # default of 10 iters not sufficient here
  iter.max = 100
)

cli::cli_alert_success("Saving CoxBoost")

saveRDS(cbfit, here::here("results/3-bladder-coxboost.rds"))


# Estimate performance ------------------------------------------------------------------------
cli::cli_alert_info("Doing bladder performance evaluation")

# Read stored models
# cooperfit <- readRDS(here::here("results/3-bladder-cooper.rds"))
# rf_c1 <- readRDS(here::here("results/3-bladder-rfsrc-c1.rds"))
# rf_c2 <- readRDS(here::here("results/3-bladder-rfsrc-c2.rds"))
# cbfit <- readRDS(here::here("results/3-bladder-coxboost.rds"))
splits <- readRDS(here::here("data/bladder-split.rds"))

# Fit CSCs with selected variables and gather results
scores_cmb <- data.table::rbindlist(list(
  fit_csc_coxph(splits, model = "cooper", coefs = selected(cooperfit, "cooper")[["1"]], cause = 1),
  fit_csc_coxph(splits, model = "cooper", coefs = selected(cooperfit, "cooper")[["2"]], cause = 2),
  fit_csc_coxph(splits, model = "coxnet", coefs = selected(cooperfit, "coxnet")[["1"]], cause = 1),
  fit_csc_coxph(splits, model = "coxnet", coefs = selected(cooperfit, "coxnet")[["2"]], cause = 2),
  fit_csc_coxph(splits, model = "rfsrc", coefs = selected(rf_c1)[["1"]], cause = 1),
  fit_csc_coxph(splits, model = "rfsrc", coefs = selected(rf_c2)[["2"]], cause = 2),
  fit_csc_coxph(splits, model = "coxboost", coefs = selected(cbfit)[["1"]], cause = 1),
  fit_csc_coxph(splits, model = "coxboost", coefs = selected(cbfit)[["2"]], cause = 2)
))

cli::cli_alert_success("Saving scores")
saveRDS(scores_cmb, here::here("results/3-bladder-scores.rds"))

