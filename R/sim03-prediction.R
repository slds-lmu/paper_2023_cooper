# Prediction ------------------------------------------------------------------------------------------------------

fwel_mt_prediction_wrapper <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 2,
    t = 100,
    a = 0.5,
    thresh = 1e-7
) {
  # suppressMessages(require(survival))
  # suppressMessages(require(riskRegression))

  checkmate::assert_numeric(alpha)
  checkmate::assert_int(mt_max_iter)
  checkmate::assert_int(t)
  checkmate::assert_numeric(a)
  checkmate::assert_numeric(thresh)

  # Only for quick debugging
  # instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 5000)
  message("Fitting fwelnet")
  fit <- fwelnet::fwelnet_mt_cox(
    instance[["data"]],
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh,
    include_mt_beta_history = TRUE
  )

  # Only for debugging
  # saveRDS(fit, "tmp-fweln-mt-testfit.rds")
  # fit <- readRDS(file = "tmp-fweln-mt-testfit.rds")

  # Use riskRegression for dummy model, needs recent GH version
  # Fit sacrifical CSC model to hold our glmnet/fwelnet coefs
  message("Fitting dummy CSC model")
  dummy_csc <- CSC(formula = Hist(time, status) ~ x1 + x2 , data = instance[["data"]], singular.ok = TRUE)

  message("Converting CSC models")
  csc_models <- convert_models_csc(train_data = instance[["data"]], fwel_fit = fit, dummy_csc = dummy_csc)

  checkmate::assert_numeric(csc_models$csc_glmnet$models[[1]]$coefficients, any.missing = FALSE, len = ncol(instance[["data"]]) - 2)
  checkmate::assert_numeric(csc_models$csc_glmnet$models[[2]]$coefficients, any.missing = FALSE, len = ncol(instance[["data"]]) - 2)
  checkmate::assert_numeric(csc_models$csc_fwelnet$models[[1]]$coefficients, any.missing = FALSE, len = ncol(instance[["data"]]) - 2)
  checkmate::assert_numeric(csc_models$csc_fwelnet$models[[2]]$coefficients, any.missing = FALSE, len = ncol(instance[["data"]]) - 2)

  message("Running Score() with AUC")
  mod_scores_auc <- Score(
    list(
      glmnet = csc_models$csc_glmnet,
      fwelnet = csc_models$csc_fwelnet
    ),
    # FIXME: Out of bounds Brier scores, the following suggestion by error msg
    # did not help and caused even worse values. AUC seems plausible though.
    # predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
    formula = Hist(time, status) ~ 1,
    data = instance[["test_data"]],
    metrics = c("AUC"),
    cause = 1,
    se.fit = TRUE, # hail mary
    times = quantile(instance[["test_data"]][["time"]], probs = seq(0.1, 0.9, .1), names = FALSE)
  )
  # in case of runtime issues: se.fit = FALSE

  message("Running Score() with Brier score")
  mod_scores_brier <- Score(
    list(
      glmnet = csc_models$csc_glmnet,
      fwelnet = csc_models$csc_fwelnet
    ),
    # FIXME: Out of bounds Brier scores, the following suggestion by error msg
    # did not help and caused even worse values. AUC seems plausible though.
    # predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
    formula = Hist(time, status) ~ 1,
    data = instance[["test_data"]],
    metrics = c("Brier"),
    summary = c("ibs", "ipa"),
    cause = 1,
    se.fit = FALSE,
    times = quantile(instance[["test_data"]][["time"]], probs = seq(0.1, 0.9, .1), names = FALSE)
  )

  message("Cleaning up results")
  # Extract Brier and AUC score tables and merge them
  scores_Brier <- mod_scores_brier$Brier$score
  auc_names <- names(mod_scores_auc$AUC$score)
  auc_newnames <- c("model", "times", "AUC", "AUC_se", "AUC_lower", "AUC_upper")
  scores_AUC <- data.table::setnames(mod_scores_auc$AUC$score, auc_names, auc_newnames)

  # Return right outer join including Brier scores + null model
  # and AUC with prefixed column names to avoid duplicate names
  res <- scores_AUC[scores_Brier, on = .(model, times)]
  return(data.table::copy(res))
}
