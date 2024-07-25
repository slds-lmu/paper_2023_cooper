#' Fit a cause-specific Cox model using `survival::coxph()` and evaluate with `riskRegression::Score()`
#'
#' @param instance Simulation data as produced by batchtools simulation function `sim_surv_binder()`.
#' @param model Model type to fit, one of "cooper", "coxboost", "coxnet", "rfsrc".
#' @param coefs Vector of model coefficients (or analogous proeprty) to use in the model fit.
#' @param cause Event of interest, defaults to `1`.
#' @param probs Quantiles to evaluate at, defaults to 10% - 80% in 10% steps.
#' @return A data.table of `riskRegression::Score` results in long format.
#'
fit_csc_coxph <- function(instance, model, coefs, cause = 1, probs = seq(0.1, 0.8, .1)) {
  require(rlang)
  train <- data.table::as.data.table(instance$train)
  test <- data.table::as.data.table(instance$test)

  checkmate::assert_data_table(train, any.missing = FALSE, min.rows = 10, min.cols = 2)
  checkmate::assert_data_table(test, any.missing = FALSE, min.rows = 10, min.cols = 2)
  checkmate::assert_numeric(probs, lower = 1e-5, upper = 1, finite = TRUE, min.len = 1)
  checkmate::assert_string(model, min.chars = 2)
  checkmate::assert_numeric(coefs, finite = TRUE, any.missing = FALSE, min.len = 1, names = "named")
  checkmate::assert_int(cause, lower = 1L, upper = 2L)
  # sub-select data
  train_cause <- train[, c("time", "status", names(coefs)), with = FALSE]
  test_cause <- test[, c("time", "status", names(coefs)), with = FALSE]

  eval_times_df <- eval_times_quant(instance$test, cause = cause, probs = probs)

  if (cause == 1L) {
    train_cause$status[train_cause$status == 2] <- 0
    test_cause$status[test_cause$status == 2] <- 0

  } else if (cause == 2L) {
    train_cause$status[train_cause$status == 1] <- 0
    test_cause$status[test_cause$status == 1] <- 0

    train_cause$status[train_cause$status == 2] <- 1
    test_cause$status[test_cause$status == 2] <- 1
  }

  # Fit coxph via survival rather than riskRegression::CSC b/c it weirdly makes a difference in Score()
  # x = TRUE includes design matrix in output, needed for Score()
  csc_model <- survival::coxph(formula = survival::Surv(time, status) ~ ., data = train_cause, x = TRUE)

  # Get quantiles of time points from full dataset to allow later aggregation per timepoint
  #  eval_times <- quantile(c(train$time, test$time), probs = eval_time_quantiles, names = FALSE)
  eval_times_df <- eval_times_quant(test, cause = cause)

  mod_scores <- riskRegression::Score(
    rlang::list2(!!model := csc_model),
    formula = Surv(time, status) ~ 1,
    data = test_cause,
    metrics = c("Brier", "AUC"),
    summary = c("ibs", "ipa"),
    # cause = cause,
    se.fit = FALSE,
    times = eval_times_df$times
  )

  scores_dat <- cleanup_score(mod_scores)
  scores_dat$cause = cause
  merge(scores_dat, eval_times_df, allow.cartesian = TRUE)
}

#' Cleanup output from riskRegression::Score to make it one long plottable data.frame
#' @param scores As returned by `riskRegression::Score()`, assumed to contain
#'   measures AUC, Brier, IBS and IPA.
#' @return A `data.table` in long format of scores by model.
cleanup_score <- function(scores) {
  result <- scores$AUC$score[scores$Brier$score, on = c("model", "times")]
  result <- data.table::melt(result , id.vars = c("model", "times"),
                             value.name = "score", variable.name = "metric",
                             meausure.vars = c("AUC", "Brier", "IBS", "IPA"))
  # Exclude some superfluous output
  result <- result[!(is.na(score) & model == "Null model"), ]
  #result <- result[!(metric == "IPA" & model == "Null model"), ]
  result
}

#' @param xdf data.frame of test data, assumed to have `time` and `status` cols
#' @param cause Event of interest, used to get event times
#' @param probs Quantiles to evaluate at, defaults to 10% - 80% in 10% steps
#'
#' @return A data.frame with columns `time` and `time_quant` for each quantile
eval_times_quant <- function(xdf, cause = 1, probs = seq(0.1, 0.8, .1)) {
  checkmate::assert_int(cause, lower = 1, upper = max(xdf$status))

  # xdf = instance$test
  event_times = xdf$time[xdf$status == cause]
  eval_times <- quantile(event_times, probs = probs, type = 2, names = FALSE)

  data.table::data.table(times = eval_times, time_quant = probs)
}

#' Partition data into train/test set with stratification
#' by status on the status column
#' @param dt A data.table with `time` and `status` columns.
#' @param train_prop Proportion of data to use for training, defaults to 70%.
#' @return a `list` with `train` and `test` data.tables
partition_dt <- function(dt, train_prop = 0.7) {
  dt <- data.table::as.data.table(dt)
  dt[, rowid := seq_along(time)]
  data.table::setkeyv(dt, "rowid")

  train <- dt[,.SD[sample(.N, ceiling(train_prop * .N))], by = status]
  test <- dt[rowid %in% setdiff(dt$rowid, train$rowid)]

  stopifnot(identical(intersect(train$rowid, test$rowid), integer(0)))

  train[, rowid := NULL]
  test[, rowid := NULL]

  list(train = train, test = test)
}
