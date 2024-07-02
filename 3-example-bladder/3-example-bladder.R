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
# Create list with $train and $test data
splits <- partition_dt(bladder, train_prop = 0.7)

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

selected(rf_c1)
selected(rf_c2)


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

cb_coefs <- coef(cbfit)
coefs_c1 <- selected(cbfit)[[1]]
coefs_c2 <- selected(cbfit)[[2]]

# Count selected per cause
length(coefs_c1)
length(coefs_c2)

# Overlap between causes?
intersect(names(coefs_c1), names(coefs_c2))


# Performances --------------------------------------------------------------------------------
cli::cli_alert_info("Doing performance evaluation")

cooperfit <- readRDS(here::here("results/3-bladder-cooper.rds"))
rf_c1 <- readRDS(here::here("results/3-bladder-rfsrc-c1.rds"))
rf_c2 <- readRDS(here::here("results/3-bladder-rfsrc-c2.rds"))
cbfit <- readRDS(here::here("results/3-bladder-coxboost.rds"))


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
scores_cmb = readRDS(here::here("results/3-bladder-scores.rds"))

p = scores_cmb |>
  mutate(
    model = dplyr::case_when(
      model == "cooper" ~ "CooPeR",
      model == "coxboost" ~ "CoxBoost",
      model == "glmnet" ~ "Coxnet",
      model == "rfsrc" ~ "RSF",
      model == "Null model" ~ "Null Model"
    ),
    model = factor(model, levels = rev(c("CooPeR", "Coxnet", "RSF", "CoxBoost", "Null Model")))
  ) |>
  filter(metric %in% c("Brier", "AUC")) |>
  ggplot(aes(x = 100 * time_quant, y = 100 * score, color = model, fill = model)) +
  facet_grid(cols = vars(cause), rows = vars(metric), scales = "free_y", labeller = label_both) +
  geom_line() +
  geom_point() +
  labs(
    title = "Bladder cancer: Performance of CSCs fit with selected variables",
    subtitle = "Performance evaluation based on 70/30 train/test split",
    x = "Time quantile (%)", y = "Score (%)",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot"
  )

ggsave(
  plot = p,
  filename = fs::path(here::here("results"), "3-bladder-performance", ext = ".png"),
  width = 6, height = 6, bg = "white"
)
