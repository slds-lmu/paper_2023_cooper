sim_wrapper_cr <- function(
    formula,
    data,
    time_grid  = seq(0, 10, by = 0.1)) {

  # baseline hazard
  instance <- sim_pexp_cr(
    formula = formula,
    data    = data,
    cut     = time_grid)
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
