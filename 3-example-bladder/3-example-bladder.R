source(here::here("R/utils.R"))
library(data.table)
library(cooper)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
if (!dir.exists(here::here("results"))) dir.create(here::here("results"))
set.seed(2023)

bladder <- readRDS(here::here("data/bladder-binder-clinical_geno.rds"))
# Create list with $train and $test data
splits <- partition_dt(bladder, train_prop = 0.7)

# Reference data
reference <- readxl::read_excel(
  "data-raw/10780432ccr062940-sup-supplemental_file_2.xls",
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

# Extract coefficients from initial Coxnet
coxnet_beta1 <- coef(cooperfit, event = 1, use_initial_fit = TRUE)
coxnet_beta2 <- coef(cooperfit, event = 2, use_initial_fit = TRUE)
# Extract final coefficients from CooPeR
cooper_beta1 <- coef(cooperfit, event = 1)
cooper_beta2 <- coef(cooperfit, event = 2)

selected <- list(
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
(cooper_shared <- intersect(names(selected$cooper$cause1), names(selected$cooper$cause2)))
# Coxnet:
(coxnet_shared <- intersect(names(selected$coxnet$cause1), names(selected$coxnet$cause2)))

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

vimps <- data.table::data.table(
  variable = rownames(rf_c1[["importance"]]),
  vi_c1 = rf_c1[["importance"]][, 1],
  vi_c2 = rf_c2[["importance"]][, 2]
)

vimps[, vita_c1 := ifelse(vi_c1 <= abs(min(vi_c1)), 0, vi_c1)]
vimps[, vita_c2 := ifelse(vi_c2 <= abs(min(vi_c2)), 0, vi_c2)]

saveRDS(vimps, here::here("results/3-bladder-rfsrc-vimps.rds"))

vimps[vita_c1 != 0 | vita_c2 != 0, ]

# Count selected per cause
nrow(vimps[vita_c1 != 0, ])
nrow(vimps[vita_c2 != 0, ])

# Overlap between causes?
vimps[vita_c1 != 0 & vita_c2 != 0, ]


# Coxboost ----------------------------------------------------------------
cli::cli_alert_info("Fitting CoxBoost")

cbfit <- coxboost_tuned(
  xdat = splits$train,
  cmprsk = "csh"
)

cli::cli_alert_success("Saving CoxBoost")

saveRDS(cbfit, here::here("results/3-bladder-coxboost.rds"))

cb_coefs <- coef(cbfit)
coefs_c1 <- cb_coefs[[1]][cb_coefs[[1]] != 0]
coefs_c2 <- cb_coefs[[2]][cb_coefs[[2]] != 0]

coefs_c1
coefs_c2

# Count selected per cause
length(coefs_c1)
length(coefs_c2)

# Overlap between causes?
intersect(names(coefs_c1), names(coefs_c2))


# Performances --------------------------------------------------------------------------------
cli::cli_alert_info("Doing performance evaluation")

rsfs = list(rf_c1, rf_c2)

scores_cmb <- data.table::rbindlist(lapply(1:2, function(cause) {
  eval_times_df <- eval_times_quant(splits$test, cause = cause)

  rr_glmnet <- refit_glmnet(cooperfit, splits$train, event = cause, alpha = 1)

  scores = Score(
    list(
      cooper = cooperfit,
      glmnet = rr_glmnet,
      rsfrc = rsfs[[cause]],
      coxboost = cbfit
    ),
    formula = Hist(time, status) ~ 1,
    data = splits$test,
    metrics = c("Brier", "AUC"),
    summary = c("ibs", "ipa"),
    se.fit = FALSE,
    cause = cause,
    times = eval_times_df$times
  )

  scores_dat <- cleanup_score(scores)
  scores_dat$cause = cause
  merge(scores_dat, eval_times_df, by = "times")
}))

cli::cli_alert_success("Saving scores")

saveRDS(scores_cmb, here::here("results/3-bladder-scores.rds"))
