source(here::here("R/utils_riskRegression.R"))

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
  if (FALSE) instance <- sim_surv_binder(n_train = 400, n_test = 200, p = 5000)

  message("Fitting fwelnet")
  fit <- fwelnet::fwelnet_mt_cox(
    instance[["train"]],
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh,
    stratify_by_status = TRUE,
    nfolds = 10,
    include_mt_beta_history = TRUE
  )

  # Only for debugging
  # saveRDS(fit, "tmp-fweln-mt-testfit.rds")
  # fit <- readRDS(file = "tmp-fweln-mt-testfit.rds")


  eval_times <- quantile(instance$test$time, probs = seq(0.1, 0.7, .1), type = 2, names = FALSE)

  message("Running Score() with AUC")

  rr_glmnet_c1 <- refit_glmnet(fit, instance$test, event = 1, alpha = alpha)

  scores <- riskRegression::Score(
    list(
      cooper = fit,
      glmnet = rr_glmnet_c1
    ),
    formula = Hist(time, status) ~ 1,
    data = instance$test,
    cause = 1,
    metrics = c("Brier", "AUC"),
    summary = c("ibs", "ipa"),
    se.fit = FALSE,
    times = eval_times
  )
  # in case of runtime issues: se.fit = FALSE


  message("Cleaning up results")
  res <- cleanup_score(scores)
  return(data.table::copy(res))
}
