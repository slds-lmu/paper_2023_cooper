source(here::here("R/utils.R"))
library(data.table)
library(cooper)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
library(dplyr)
library(ggplot2)

if (!dir.exists(here::here("results"))) dir.create(here::here("results"))
set.seed(2023)

bladder <- readRDS(here::here("data/bladder-binder-clinical_geno.rds"))

# Reference data
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

# Extract coefficients from initial Coxnet
coxnet_beta1 <- coef(cooperfit, event = 1, use_initial_fit = TRUE)
coxnet_beta2 <- coef(cooperfit, event = 2, use_initial_fit = TRUE)
# Extract final coefficients from CooPeR
cooper_beta1 <- coef(cooperfit, event = 1)
cooper_beta2 <- coef(cooperfit, event = 2)

selected_vars <- list(
  cooper = list(
    cause1 = cooper_beta1[cooper_beta1 != 0],
    cause2 = cooper_beta2[cooper_beta2 != 0]
  ),
  coxnet = list(
    cause1 = coxnet_beta1[coxnet_beta1 != 0],
    cause2 = coxnet_beta2[coxnet_beta2 != 0]
  )
)

# Covariables shared between causes for
# CooPeR:
(cooper_shared <- intersect(
  names(selected_vars$cooper$cause1),
  names(selected_vars$cooper$cause2)
))
# Coxnet:
(coxnet_shared <- intersect(
  names(selected_vars$coxnet$cause1),
  names(selected_vars$coxnet$cause2)
))

coxnet_beta1[names(coxnet_beta1) == "age"]
coxnet_beta2[names(coxnet_beta2) == "age"]
cooper_beta1[names(cooper_beta1) == "age"]
cooper_beta2[names(cooper_beta2) == "age"]

cooper_shared[cooper_shared %in% reference$probe_id]

reference |>
  dplyr::filter(probe_id %in% cooper_shared)

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

selected(rf_c1)
selected(rf_c2)


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

cb_coefs <- coef(cbfit)
coefs_c1 <- selected(cbfit)[[1]]
coefs_c2 <- selected(cbfit)[[2]]

# Count selected per cause
length(coefs_c1)
length(coefs_c2)

# Overlap between causes?
intersect(names(coefs_c1), names(coefs_c2))


# View selected vars -----------------------------------------------------------------------------

cooperfit <- readRDS(here::here("results/3-bladder-varsel-cooper.rds"))
rf_c1 <- readRDS(here::here("results/3-bladder-varsel-rfsrc-c1.rds"))
rf_c2 <- readRDS(here::here("results/3-bladder-varsel-rfsrc-c2.rds"))
cbfit <- readRDS(here::here("results/3-bladder-varsel-coxboost.rds"))

