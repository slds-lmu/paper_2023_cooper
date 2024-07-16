source(here::here("R/utils.R"))
library(data.table)
library(cooper)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)

bladder_file <- here::here("data/bladder-binder-clinical_geno.rds")

if (!file.exists(bladder_file)) {
  message("Recreating bladder dataset from data-raw/preprocess-binder.R")
  source(here::here("data-raw/preprocess-binder.R"))
}

if (!dir.exists(here::here("results"))) dir.create(here::here("results"))
set.seed(2023)

bladder <- readRDS(bladder_file)
# Create list with $train and $test data
splits <- partition_dt(bladder, train_prop = 0.7)

saveRDS(splits, here::here("data/bladder-split.rds"))

# Reference data
reference <- readxl::read_excel(
  here::here("data-raw/10780432ccr062940-sup-supplemental_file_2.xls"),
  sheet = "Progression classifier probes"
) |>
  janitor::clean_names()

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