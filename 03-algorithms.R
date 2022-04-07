# Multi-task fwelnet-----------------------------------------------------------
fwel_mt_wrapper <- function(
    data, job, instance,
    alpha = 1, z_scale = 1, z_method = "original",
    theta = "original",
    mt_max_iter = 2,
    t = 1, a = 0.5
    ) {

  # batchtools stores possible theta values in a factor column, so we need
  # to convert "original" to NULL as fwelnet() expects
  if (theta == "original") {
    theta <- NULL
  } else {
    # .. and have to coerce the factor value to a character first and then to
    # numeric, otherwise we would get the numeric factor level instead of the
    # actual value we want.
    theta <- as.numeric(as.character(theta))
  }

  fwelnet::fwelnet_mt_cox(
    instance$data, mt_max_iter = mt_max_iter,
    z_scale = z_scale, z_method = as.character(z_method),
    alpha = alpha,
    theta = theta,
    t = t,
    a = a
  )
}
