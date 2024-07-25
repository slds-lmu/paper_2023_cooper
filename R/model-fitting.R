



#' Refit a GLMnet with riskRegression using the same predictors and lambda as a cooper fit.
#' This allows evaluating the glmnet internally fit within cooper() with `riskRegression::Score()``,
#' which is a workaround to use `Score()` for a cooper model and its internal glmnet simultaneously,
#' using the exact same data and parameters for each.
#' @param cooperfit A model fit of class `"cooper"` as returned by `cooper()`.
#' @return A `riskRegression::GLMnet()` model fit.
refit_glmnet <- function(cooperfit, xnew, event = 1, alpha = 1) {
  checkmate::assert_class(cooperfit, "cooper")
  checkmate::assert_int(event, lower = 1)
  require(survival)

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



#' Tune rfsrc hyperparameters using `randomForestSRC::tune()`
#' @param xdat data.frame with `time` and `status` columns
#' @param splitrule,importance Passed to `randomForestSRC::tune()` and `randomForestSRC::rfsrc()`
#' @param cause Event of interest, used to set `cause` argument in `rfsrc()`
rfsrc_tuned <- function(xdat, splitrule = "logrank", importance = "random", cause = 1, ...) {
  require("survival")

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
#' @param iter.max Passed to `CoxBoost::optimCoxBoostPenalty()`
# cbres = coxboost_tuned(instance$train)
coxboost_tuned <- function(xdat, cmprsk = "csh", iter.max = 10, ...) {
  xdat <- data.table::as.data.table(xdat)
  # remove time/status from predictor matrix
  outcome_vars <-  !(colnames(xdat) %in% c("time", "status"))
  xmat <- xdat[, ..outcome_vars]
  xmat = as.matrix(xmat)
  checkmate::assert_matrix(xmat, mode = "numeric")

  optim = CoxBoost::optimCoxBoostPenalty(
    time = xdat$time,
    status = xdat$status,
    x = xmat,
    cmprsk = cmprsk,
    iter.max = iter.max,
    # Passed to cv.CoxBoost()
    # indicates whether computations in the cross-validation folds should be performed in parallel,
    #  using package parallel. If TRUE, package parallel is employed using the default number of cores.
    #  A value larger than 1 is taken to be the number of cores that should be employed.
    multicore = getOption("mc.cores", default = 1)
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

