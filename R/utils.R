sim_wrapper_cr <- function(
    formula,
    n = 1000,
    time_grid  = seq(0, 10, by = 0.1)) {

  # create data set with covariates
  xdf1 <- tibble::tibble(
    x0 = sample(c(-1,1), n, .3),
    x1 = runif(n, -3, 3),
    x2 = runif(n, -3, 3),
    x3 = runif(n, -3, 3))
  xdf2 <- mvtnorm::rmvnorm(n = nrow(xdf1), mean = rep(0, 10))
  # noise variables
  colnames(xdf2) <- paste0("x", 4:(ncol(xdf2)+3))
  xdf <- cbind(xdf1, xdf2)

  # baseline hazard
  instance <- sim_pexp_cr(
    formula = formula,
    data    = xdf,
    cut     = time_grid
  )
  instance$status <- instance$type
  instance$type <- NULL
  # add censoring
  cens_times <- runif(nrow(instance), 0, 20)
  cens <- instance$time > cens_times
  instance$time <- pmin(instance$time, cens_times)
  instance$status[cens] <- 0
  instance$status[instance$time == 10] <- 0
  instance <- as.data.frame(instance)
  instance$id <- instance$hazard1 <- instance$hazard2 <- NULL
  instance$x0 <- as.factor(instance$x0)

  instance
}


# Copied verbatim from https://raw.githubusercontent.com/adibender/pem.xgb/master/R/sim-cr.R
# to reduce dependencies, and explicitly namespaced functions
#' Simulate competing risks time-to-event data via piece-wise exponential distribution
#'
#' Piece-wise exponential implementation of simulation algorithm described in
#' Beyersmann et al. 2009 (<doi: 10.1002/sim.3516>).
#'
#' @inheritParams pammtools::sim_pexp
#' @importFrom Formula Formula
#' @importFrom rlang is_atomic
#' @importFrom purrr reduce
#' @examples
#' library(pammtools)
#' library(dplyr)
#' library(survival)
#' # Example from Competing Risks book for comparison
#' simul.dat.cp <- function(n, h01, h02, cens.param) {
#'   times <- rexp(n, h01 + h02)
#'   ev <- rbinom(n, size = 1, prob = h01 / (h01 + h02))
#'   ev <- ifelse(ev == 0, 2, 1)
#'   cens.time <- runif(n, cens.param[1], cens.param[2])
#'   obs.time <- pmin(times, cens.time)
#'   obs.cause <- as.numeric(times == obs.time) * ev
#'   data.frame(obs.time, obs.cause)
#' }
#' set.seed(923)
#' n <- 1000
#' dat.exo7 <- simul.dat.cp(n, 0.5, 0.9, c(0.5, 3))
#' # We compute the Nelson-Aalen estimates using survfit()
#' # That’s just to get the number of events and the risk set
#' temp01 <- survfit(Surv(obs.time, obs.cause == 1) ~ 1, dat.exo7)
#' temp02 <- survfit(Surv(obs.time, obs.cause == 2) ~ 1, dat.exo7)
#' na01 <- cumsum(temp01$n.event / temp01$n.risk)
#' na02 <- cumsum(temp02$n.event / temp02$n.risk)
#'
#' plot(temp01$time, na02, type="s", ylab="Cumulative transition hazard", xlab="time")
#' lines(temp01$time, na01, type="s", col=2)
#'
#' # create data set with variables which will affect the hazard rate
#' # (not used here, but usually more complex examples of interest)
#' df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
#'  as_tibble()
#' set.seed(24032018)
#' df <- cbind.data.frame(
#'   x1 = runif (n, -3, 3),
#'   x2 = runif (n, 0, 6))
#' # two component formula specifying cause specific hazards
#' form <- ~ log(0.5)| log(0.9)
#' sim_df <- sim_pexp_cr(form, df, seq(0, 3, by =.25)) %>%
#'  mutate(
#'   cens_time = runif(n(), 0.5, 3),
#'   status = if_else(cens_time < time, 0, 1),
#'   time = pmin(time, cens_time),
#'   type = status * type)
#' temp01_2 <- survfit(Surv(time,type == 1) ~ 1, sim_df)
#' temp02_2 <- survfit(Surv(time, type == 2) ~ 1, sim_df)
#' na01_2 <- cumsum(temp01_2$n.event / temp01_2$n.risk)
#' na02_2 <- cumsum(temp02_2$n.event / temp02_2$n.risk)
#'
#' lines(temp01_2$time, na02_2, type="s", col = 3)
#' lines(temp01_2$time, na01_2, type="s", col = 4)
#' @export
sim_pexp_cr <- function(formula, data, cut) {

  Form <- Formula::Formula(formula)
  F_rhs <- attr(Form, "rhs")
  l_rhs <- length(F_rhs)
  seq_rhs <- seq_len(l_rhs)

  data <- data |>
    dplyr::mutate(
      id     = dplyr::row_number(),
      time   = max(cut),
      status = 1)

  # construct eta for time-constant part
  ped  <- pammtools::split_data(
    formula = survival::Surv(time, status)~.,
    data    = dplyr::select_if (data, rlang::is_atomic),
    cut     = cut,
    id      = "id") |>
    dplyr::rename("t" = "tstart")

  # calculate cause specific hazards
  for(i in seq_rhs) {
    ped[[paste0("hazard", i)]] <-  exp(eval(F_rhs[[i]], ped))
  }
  ped[["rate"]] <- purrr::reduce(ped[paste0("hazard", seq_rhs)], `+`)

  # simulate survival times
  sim_df <- ped |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      time   = pammtools:::rpexp(rate = .data$rate, t = .data$t),
      status = 1L * (.data$time <= max(cut)),
      time   = pmin(.data$time, max(cut))) |>
    dplyr::filter(.data$t < .data$time & .data$time <= .data$tend)

  sim_df$type <- apply(sim_df[paste0("hazard", seq_rhs)], 1,
                       function(probs)
                         sample(seq_rhs, 1, prob = probs))

  sim_df |>
    dplyr::mutate(type = ifelse(.data$status == 1, .data$type, 0)) |>
    dplyr::select(-dplyr::one_of(c("t", "tend", "interval", "offset", "ped_status", "rate")))

}


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
# cbres = coxboost_tuned(instance$train)
coxboost_tuned <- function(xdat, cmprsk = "csh", ...) {
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

#' Partition data into train/test set with stratification
#' by status on the status column
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
