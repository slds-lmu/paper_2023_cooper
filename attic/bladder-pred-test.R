# Test prediction on bladder cancer data

# Bladder cancer data from bimj2384-sup-0001-suppmat, survival and genetic data only, standardized.
# FIXME: When doing resampling, standardize as part of the resampling pipeline
surv_geno_std <- readRDS(here::here("data/surv_geno_std.rds"))
dim(surv_geno_std)
table(surv_geno_std$status)

tsk <- mlr3::as_task_classif(surv_geno_std, target = "status", id = "bladder")
tsk$set_col_roles("status", add_to = "stratum")
mlr3::rsmp("holdout")

row_ids <- seq_len(nrow(surv_geno_std))
train_ids <- sample(row_ids, size = floor(2/3 * nrow(surv_geno_std)))
test_ids <- setdiff(row_ids, train_ids)

stopifnot(setequal(union(train_ids, test_ids), row_ids))

surv_geno_train <- surv_geno_std[train_ids, ]
surv_geno_test <- surv_geno_std[test_ids, ]

mt_max_iter <- 3
fit <- fwelnet::fwelnet_mt_cox(
  surv_geno_train,
  mt_max_iter = 3,
  alpha = 1,
  t = 100,
  a = 0.5,
  thresh = 1e-7,
  include_mt_beta_history = TRUE
)

# str(fit)

nzero_c1 <- names(which(fit$beta1[, mt_max_iter + 1] > 0))
nzero_c2 <- names(which(fit$beta2[, mt_max_iter + 1] > 0))

lengths(list(nzero_c1, nzero_c2))
setdiff(nzero_c1, nzero_c2) |> length()
union(nzero_c1, nzero_c2) |> length()

library(riskRegression)
library(survival)

data_fwel_selected <- surv_geno_train[, c("time", "status", nzero_c1)]
csc_1_fwelnet <- CSC(formula = Hist(time, status) ~ ., data = data_fwel_selected, cause = 1)

# glmnet-selected data
nzero_c1_glmnet <- names(which(fit$beta1[, 1] > 0))
nzero_c2_glmnet <- names(which(fit$beta2[, 1] > 0))

lengths(list(nzero_c1_glmnet, nzero_c2_glmnet))
setdiff(nzero_c1_glmnet, nzero_c2_glmnet) |> length()
union(nzero_c1_glmnet, nzero_c2_glmnet) |> length()

data_glmnet_selected <- surv_geno_train[, c("time", "status", nzero_c1_glmnet)]
csc_1_glmnet <- CSC(formula = Hist(time, status) ~ ., data = data_glmnet_selected, cause = 1)

# selected by fwelnet
sort(nzero_c1)
# selected by glmnet
sort(nzero_c1_glmnet)

# selected by fwelnet
sort(nzero_c2)
# selected by glmnet
sort(nzero_c2_glmnet)


mod_scores_auc <- Score(
  list(
    glmnet = csc_1_glmnet,
    fwelnet = csc_1_fwelnet
  ),
  formula = Hist(time, status) ~ 1,
  data = surv_geno_test,
  metrics = c("AUC"),
  cause = 1,
  se.fit = TRUE, # hail mary
  times = quantile(surv_geno_std[["time"]], probs = seq(0.1, 0.9, .1), names = FALSE)
)

mod_scores_brier <- Score(
  list(
    glmnet = csc_1_glmnet,
    fwelnet = csc_1_fwelnet
  ),
  # FIXME: Out of bounds Brier scores, the following suggestion by error msg
  # did not help and caused even worse values. AUC seems plausible though.
  predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
  formula = Hist(time, status) ~ 1,
  data = surv_geno_test,
  metrics = c("Brier"),
  summary = c("ibs", "ipa"),
  cause = 1,
  se.fit = FALSE,
  times = quantile(surv_geno_std[["time"]], probs = seq(0.1, 0.9, .1), names = FALSE)
)

#ggplot2::autoplot(mod_scores_auc)
mod_scores_brier$Brier$score |>
  ggplot(aes(y = Brier, x = times, color = model)) +
    geom_path() +
    labs(title = "bladder cancer survival prediction", subtitle = "Based on fwelnet feature selection and CSC cause-1 prediction") +
    theme_minimal()

mod_scores_auc$AUC$score |>
  ggplot(aes(y = AUC, x = times, color = model)) +
  geom_path() +
  labs(title = "bladder cancer survival prediction", subtitle = "Based on fwelnet feature selection and CSC cause-1 prediction") +
  theme_minimal()



# RFSRC -----------------------------------------------------------------------------------------------------------



library(randomForestSRC)

rf_c1 <- rfsrc(
  Surv(time, status) ~ ., data = surv_geno_std,
  splitrule = "logrank",
  cause = c(1, 0),
  mtry = 500,
  nodesize = 20
)

rf_c2 <- rfsrc(
  Surv(time, status) ~ ., data = surv_geno_std,
  splitrule = "logrank",
  cause = c(0, 1),
  mtry = 500,
  nodesize = 20
)
