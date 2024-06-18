#' Simple status distribution check
#' @param xs A status variable
#' @examples
#'
#' check_status(instance$train$status)
check_status <- function(xs) {
  tab <- as.integer(table(xs))
  mtab <- rbind(events = tab, prop = 100 * round(tab/sum(tab), 2))
  colnames(mtab) <- unique(xs)
  mtab <- cbind(mtab, sum = rowSums(mtab))
  print(mtab)
}


#' Cleanup output from riskRegression::Score to make it one long plottable data.frame
#' @param scores As returned by `riskRegression::Score()`, assumed to contain
#'   measures AUC, Brier, IBS and IPA.
cleanup_score <- function(scores) {
  result <- scores$AUC$score[scores$Brier$score, on = c("model", "times")]
  result <- data.table::melt(result , id.vars = c("model", "times"),
                             value.name = "score", variable.name = "metric",
                             meausure.vars = c("AUC", "Brier", "IBS", "IPA"))
  # Exclude some superfluous output
  result <- result[!(is.na(score) & model == "Null model"), ]
  result <- result[!(metric == "IPA" & model == "Null model"), ]
  result
}

#' Refit a GLMnet with riskRegression using the same predictors and lambda as a cooper fit.
#' This allows evaluating the glmnet internally fit within cooper() with riskRegression::Score(),
#' which is a workaround to use Score() for a cooper model and its internal glmnet simultaneously,
#' using the exact same data and parameters for each.
#' @param cooperfit A model fit of class `"cooper"` as returned by `cooper()`.
refit_glmnet <- function(cooperfit, xnew, event = 1, alpha = 1) {
  checkmate::assert_class(cooperfit, "cooper")
  checkmate::assert_int(event, lower = 1)

  # Recode to censor other event and ensure status in 0,1
  if (event == 1) {
    xnew$status[xnew$status == 2] <- 0
  }
  if (event == 2) {
    xnew$status[xnew$status == 1] <- 0
    xnew$status[xnew$status == 2] <- 1
  }

  riskRegression::GLMnet(
    formula = stats::reformulate(cooperfit$predictors, response = "Surv(time, status)"),
    data = xnew,
    lambda = cooperfit$initial_fits[[event]]$lambda.min,
    cv = FALSE,
    standardize = TRUE, alpha = alpha
  )
}


# S3 method to have riskRegression's predictRisk work with CoxBoost
predictRisk.CoxBoost <- function(object, newdata, times, cause, ...) {
  checkmate::assert_int(cause, lower = 1)
  checkmate::assert_numeric(times, lower = 0, finite = TRUE, min.len = 1)

  # Returns list with predictions named "1" and "2"
  # Entries are vector if times is a scalar, a matrix otherwise
  res = predict(object, newdata = newdata, times = times, type = "CIF")

  # only keep the cause of interest, in matrix form
  res = matrix(res[[as.character(cause)]], ncol = length(times))

  # riskRegression expects colnames to be eval times
  colnames(res) <- times
  res
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

#' Tune rfsrc hyperparameters using `randomForestSRC::tune()`
#' @param xdat data.frame with `time` and `status` columns
#' @param splitrule,importance Passed to `randomForestSRC::tune()` and `randomForestSRC::rfsrc()`
#' @param cause Event of interest, used to set `cause` argument in `rfsrc()`
rfsrc_tuned <- function(xdat, splitrule = "logrank", importance = "random", cause = 1, ...) {

  cause_arg = switch(
    cause,
    `1` = c(1, 0),
    `2` = c(0, 1),
    stop("cause must be 1 or 2")
  )

  tuned = randomForestSRC::tune(
    Surv(time, status) ~ ., data = xdat,
    splitrule = splitrule,
    cause = cause_arg,
    importance = importance,
    # return rf with best params if TRUE, but this does not correctly pass splitrule param!
    # doBest = TRUE,
    ...
  )

  randomForestSRC::rfsrc(
    Surv(time, status) ~ ., data = xdat,
    splitrule = splitrule,
    cause = cause_arg,
    importance = importance,
    mtry = tuned$optimal[["mtry"]],
    nodesize = tuned$optimal[["nodesize"]],
    ...
  )

}


#' CoxBoost wrapper that tunes penalty and step size using `optimCoxBoostPenalty()`
#' @param xdat data.frame with `time` and `status` columns
#' @param cmprsk Passed to `CoxBoost` and `CoxBoost::optimCoxBoostPenalty`
# cbres = coxboost_tuned(instance$train)
coxboost_tuned <- function(xdat, cmprsk = "csh", ...) {
  xmat <- xdat[!(colnames(xdat) %in% c("time", "status"))]
  xmat = as.matrix(xdat)

  optim = CoxBoost::optimCoxBoostPenalty(
    time = xdat$time,
    status = xdat$status,
    x = xmat,
    cmprsk = cmprsk,
    ...
  )

  CoxBoost::CoxBoost(
    time = xdat$time,
    status = xdat$status,
    x = xmat,
    stepno = optim$cv.res$optimal.step,
    penalty = optim$penalty,
    cmprsk = cmprsk,
    ...
  )
}

