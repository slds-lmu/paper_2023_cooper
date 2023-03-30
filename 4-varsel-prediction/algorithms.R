# Variable selection ----------------------------------------------------------------------------------------------

# fwelnet wrapper for variable selection sim
fwel_mt_varselect_pred <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 2,
    t = 100,
    a = 0.5,
    thresh = 1e-7
) {

  checkmate::assert_numeric(alpha)
  checkmate::assert_int(mt_max_iter)
  checkmate::assert_int(t)
  checkmate::assert_numeric(a)
  checkmate::assert_numeric(thresh)

  fit <- fwelnet::fwelnet_mt_cox(
    instance$train,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh,
    include_mt_beta_history = TRUE
  )

  p <- ncol(instance$train)

  glmnet_beta1 <- fit$beta1[, 1]
  glmnet_beta2 <- fit$beta2[, 1]

  # Same for fwelnet estimates
  fwel_beta1 <- fit$beta1[, mt_max_iter + 1]
  fwel_beta2 <- fit$beta2[, mt_max_iter + 1]

  fw_coefs <- list(
    fwelnet = list(
      cause1 = fwel_beta1[fwel_beta1 != 0],
      cause2 = fwel_beta2[fwel_beta2 != 0]
    ),
    glmnet = list(
      cause1 = glmnet_beta1[glmnet_beta1 != 0],
      cause2 = glmnet_beta2[glmnet_beta2 != 0]
    )
  )

  scores <- data.table::rbindlist(list(
    fit_csc(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause1, cause = 1),
    fit_csc(instance$train, instance$test, model = "fwelnet", coefs = fw_coefs$fwelnet$cause2, cause = 2),
    fit_csc(instance$train, instance$test, model = "glmnet", coefs = fw_coefs$glmnet$cause1, cause = 1),
    fit_csc(instance$train, instance$test, model = "glmnet", coefs = fw_coefs$glmnet$cause2, cause = 2)
  ))

  list(scores = scores, coefs = fw_coefs)

}


rfsrc_varselect_pred <- function(data, job, instance,
                                    mtry = 2000,
                                    nodesize = 30, splitrule = "logrank",
                                    importance = "random", cutoff_method = "vita"
                                    ) {
  # batchtools passes params as factors, need to convert
  importance <- as.character(importance)
  splitrule <- as.character(splitrule)
  cutoff_method <- as.character(cutoff_method)

  # prbly won't use permute, since random should be equivalent
  checkmate::assert_subset(importance, choices = c("random", "permute", "anti"))
  # not sure about other methods yet
  checkmate::assert_subset(cutoff_method, choices = "vita")

  rf_c1 <- rfsrc(
    Surv(time, status) ~ ., data = instance$train,
    splitrule = splitrule,
    cause = c(1, 0),
    importance = importance,
    mtry = mtry,
    nodesize = nodesize
  )

  rf_c2 <- rfsrc(
    Surv(time, status) ~ ., data = instance$train,
    splitrule = splitrule,
    cause = c(0, 1),
    importance = importance,
    mtry = mtry,
    nodesize = nodesize
  )

  vimps <- data.table::data.table(
    variable = rownames(rf_c1[["importance"]]),
    vi_c1 = rf_c1[["importance"]][, 1],
    vi_c2 = rf_c2[["importance"]][, 2]
  )

  # If these fail something is weird, but it would justify using abs(min(x)) rather than -min(x)
  checkmate::assert_true(any(vimps$vi_c1 < 0))
  checkmate::assert_true(any(vimps$vi_c2 < 0))

  if (cutoff_method == "vita") {
    # Apply vita/janitza shortcut
    # Take absolute value of minimal importance value as cut-off for 0-classification
    # Effectively equivalent with vita method in Degenhardt et al. (2019)
    # Selecting
    vimps[, vita_c1 := ifelse(vi_c1 <= abs(min(vi_c1)), 0, vi_c1)]
    vimps[, vita_c2 := ifelse(vi_c2 <= abs(min(vi_c2)), 0, vi_c2)]
  }

  res_c1 <- vimps$vita_c1
  res_c2 <- vimps$vita_c2
  names(res_c1) <- vimps$variable
  names(res_c2) <- vimps$variable

  # Since VIMPs are >=0 we only select the ones greater 0
  res_c1 <- res_c1[res_c1 > 0]
  res_c2 <- res_c2[res_c2 > 0]

  rf_coefs <- list(
    cause1 = res_c1,
    cause2 = res_c2
  )

  scores <- data.table::rbindlist(list(
    fit_csc(instance$train, instance$test, model = "rfsrc", coefs = rf_coefs$cause1, cause = 1),
    fit_csc(instance$train, instance$test, model = "rfsrc", coefs = rf_coefs$cause2, cause = 2)
  ))

  list(scores = scores, coefs = rf_coefs)
}

coxboost_varselect_pred <- function(data, job, instance,
                                       cmprsk = "csh", stepno = 100, penalty = 2000) {
  cbfit <- CoxBoost(
    time = instance$train$time,
    status = instance$train$status,
    x = as.matrix(instance$train[, -c(1, 2)]),
    cmprsk = as.character(cmprsk),
    stepno = stepno,
    # Default is 9 * sum(status[subset] == 1), should be approx 9 * 250
    # since in our setting sum(instance$train$status != 0) is approx 233
    penalty = penalty
  )

  # Extracts coefs at final boosting step, names list per cause (names 1, 2)
  cb_coefs <- coef(cbfit)

  if (cmprsk == "sh") {
    # No cause-specific coefficients with subdistribution approach, so duplicate results (kind of)
    res_c1 <- cb_coefs
    res_c2 <- cb_coefs
  } else {
    res_c1 <- cb_coefs[["1"]]
    res_c2 <- cb_coefs[["2"]]
  }

  # Only get nonzeros
  res_c1 <- res_c1[res_c1 != 0]
  res_c2 <- res_c2[res_c2 != 0]

  cb_coefs <- list(
    cause1 = res_c1,
    cause2 = res_c2
  )

  scores <- data.table::rbindlist(list(
    fit_csc(instance$train, instance$test, model = "coxboost", coefs = cb_coefs$cause1, cause = 1),
    fit_csc(instance$train, instance$test, model = "coxboost", coefs = cb_coefs$cause2, cause = 2)
  ))

  list(scores = scores, coefs = cb_coefs)
}

fit_csc <- function(train, test, model, coefs, cause = 1) {

  checkmate::assert_data_table(train, any.missing = FALSE, min.rows = 10, min.cols = 2)
  checkmate::assert_data_table(test, any.missing = FALSE, min.rows = 10, min.cols = 2)
  checkmate::assert_string(model, min.chars = 2)
  checkmate::assert_numeric(coefs, finite = TRUE, any.missing = FALSE, min.len = 1)
  checkmate::assert_subset(cause, c(1, 2))
  # sub-select data
  train_cause <- train[, c("time", "status", names(coefs)), with = FALSE]
  test_cause <- test[, c("time", "status", names(coefs)), with = FALSE]

  csc_model <- CSC(formula = Hist(time, status) ~ ., data = train_cause, cause = cause)

  # Get quantiles of time points from full dataset to allow later aggregation per timepoint
  eval_times <- quantile(c(train$time, test$time), probs = seq(0.1, 0.9, .1), names = FALSE)

  # Yikes
  # Estimated risk outside the range [0,1].
  # Possible cause: incorrect extrapolation, i.e., time and/or covariates used for the prediction differ from those used to fit the Cox models.
  mod_scores_auc <- Score(
    list2(!!model := csc_model), # uses rlang splicing to get list(fwelnet = csc_1)
    predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
    formula = Hist(time, status) ~ 1,
    data = test_cause,
    metrics = "AUC",
    cause = cause,
    se.fit = FALSE,
    times = eval_times
  )


  mod_scores_brier <- Score(
    list2(!!model := csc_model),
    predictRisk.args = list(CauseSpecificCox = list(product.limit = FALSE)),
    formula = Hist(time, status) ~ 1,
    data = test_cause,
    metrics = c("Brier"),
    summary = c("ibs", "ipa"),
    cause = cause,
    se.fit = FALSE,
    times = eval_times
  )


  result <- mod_scores_auc$AUC$score[mod_scores_brier$Brier$score, on = c("model", "times")]
  result <- data.table::melt(result, id.vars = c("model", "times"),
                             value.name = "score", variable.name = "metric",
                             meausure.vars = c("AUC", "Brier", "IBS", "IPA"))

  result$cause <- cause
  result
}
