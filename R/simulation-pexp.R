# Simulation functions ------------------------------------------------------------------------

# Copied verbatim from https://raw.githubusercontent.com/adibender/pem.xgb/master/R/sim-cr.R
# to reduce dependencies, and explicitly namespaced functions
#' Simulate competing risks time-to-event data via piece-wise exponential distribution
#'
#' Piece-wise exponential implementation of simulation algorithm described in
#' Beyersmann et al. 2009 (<doi: 10.1002/sim.3516>).
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

#' Wrapper to simulate survival data using piecewise-exponential process
#' Used in proof of concept simulation
#' @param formula Passed to `sim_pexp_cr()`.
#' @param n Number of observations to simulate
#' @param time_grid Grid of time points ot use for survival times.
#'
#' @return A `data.frame` of competing risk survival data with `time`,`status` variables and 14 features.
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
