invisible(lapply(list.files("R", pattern = "*.R", full.names = TRUE), source, echo = FALSE))
library(data.table)
library(cooper)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
library(dplyr)
library(ggplot2)
if (!dir.exists(here::here("results"))) dir.create(here::here("results"))
set.seed(2023)

bladder_file <- here::here("data/bladder-binder-clinical_geno.rds")
if (!file.exists(bladder_file)) {
  cli::cli_alert_info("Recreating bladder dataset from data-raw/preprocess-bladder-data.R")
  source(here::here("data-raw/preprocess-bladder-data.R"))
}
bladder <- readRDS(bladder_file)

# Reference data
# Downloaded from supp table 2 (xlsx) https://aacrjournals.org/clincancerres/article/13/12/3545/13137/Gene-Expression-Signatures-Predict-Outcome-in-Non
reference <- readxl::read_excel(
  here::here("data-raw/10780432ccr062940-sup-supplemental_file_2.xls"),
  sheet = "Progression classifier probes"
) |>
  janitor::clean_names()


# CooPeR --------------------------------------------------------------------------------------
cli::cli_alert_info("Fitting CooPeR")

cooperfit <- cooper::cooper(
  bladder, mt_max_iter = 3,
  alpha = 1, t = 100, thresh = 1e-7,
  stratify_by_status = TRUE, nfolds = 5
)

cli::cli_alert_success("Saving CooPeR model")

saveRDS(cooperfit, here::here("results/3-bladder-varsel-cooper.rds"))

# rfsrc ---------------------------------------------------------------------------------------
cli::cli_alert_info("Fitting rfsrc")

rf_c1 <- rfsrc_tuned(
  xdat = bladder,
  splitrule = "logrank",
  importance = "random",
  cause = 1
)

rf_c2 <- rfsrc_tuned(
  xdat = bladder,
  splitrule = "logrank",
  importance = "random",
  cause = 2
)

cli::cli_alert_success("Saving rfsrc models")

saveRDS(rf_c1, here::here("results/3-bladder-varsel-rfsrc-c1.rds"))
saveRDS(rf_c2, here::here("results/3-bladder-varsel-rfsrc-c2.rds"))

# Coxboost ----------------------------------------------------------------
cli::cli_alert_info("Fitting CoxBoost")

cbfit <- coxboost_tuned(
  xdat = bladder,
  cmprsk = "csh",
  # default of 10 iters not sufficient here
  iter.max = 100
)

cli::cli_alert_success("Saving CoxBoost")

saveRDS(cbfit, here::here("results/3-bladder-varsel-coxboost.rds"))

# View selected vars -----------------------------------------------------------------------------

cooperfit <- readRDS(here::here("results/3-bladder-varsel-cooper.rds"))
rf_c1 <- readRDS(here::here("results/3-bladder-varsel-rfsrc-c1.rds"))
rf_c2 <- readRDS(here::here("results/3-bladder-varsel-rfsrc-c2.rds"))
cbfit <- readRDS(here::here("results/3-bladder-varsel-coxboost.rds"))

# Extract coefficients from initial Coxnet
coxnet_beta1 <- coef(cooperfit, event = 1, use_initial_fit = TRUE)
coxnet_beta2 <- coef(cooperfit, event = 2, use_initial_fit = TRUE)
# Extract final coefficients from CooPeR
cooper_beta1 <- coef(cooperfit, event = 1)
cooper_beta2 <- coef(cooperfit, event = 2)

# Covariables shared between causes for
# CooPeR:
(cooper_shared <- intersect(
  names(selected(cooperfit, model = "cooper")[[1]]),
  names(selected(cooperfit, model = "cooper")[[2]])
))
# Coxnet:
(coxnet_shared <- intersect(
  names(selected(cooperfit, model = "coxnet")[[1]]),
  names(selected(cooperfit, model = "coxnet")[[2]])
))

cooper_shared[cooper_shared %in% reference$probe_id]
reference$probe_id[reference$probe_id %in% cooper_shared]

## rfsrc
selected(rf_c1)
selected(rf_c2)

intersect(
  names(selected(rf_c1)),
  names(selected(rf_c2))
)

# CoxBoost
cb_coefs <- coef(cbfit)
coefs_c1 <- selected(cbfit)[[1]]
coefs_c2 <- selected(cbfit)[[2]]

# Count selected per cause
length(coefs_c1)
length(coefs_c2)

# Overlap between causes?
intersect(
  names(selected(cbfit)[[1]]),
  names(selected(cbfit)[[2]])
)
