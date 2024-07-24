#' Wrapper for cooper with coefficient error measurement.
#' @param data,job,instance Parameters required by batchtools.
#' @param alpha,z_method,mt_max_iter,t,thresh Parameters for cooper.
cooper_wrapper <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 5,
    t = 1,
    thresh = 1e-3
) {

  fit <- cooper::cooper(
    instance$data,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    stratify_by_status = TRUE,
    alpha = alpha,
    t = t,
    a = 0.5,
    thresh = thresh,
    include_mt_beta_history = FALSE
  )

  res = data.table::rbindlist(lapply(1:2, \(event) {
    rbind(
      data.table::data.table(
        variable = fit$predictors,
        estimate = coef(fit, event = event, use_initial_fit = TRUE),
        cause = event,
        method = "coxnet"
      ),
      data.table::data.table(
        variable = fit$predictors,
        estimate = coef(fit, event = event),
        cause = event,
        method = "cooper"
      )
    )
  }))

  res = merge(res, instance$effects, by = c("variable", "cause"))
  res[, error := (truth - estimate)]

}
