# Multi-task fwelnet-----------------------------------------------------------
fwel_mt_wrapper <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 2,
    t = 1, a = 0.5, thresh = 1e-3
    ) {

  fwelnet::fwelnet_mt_cox(
    instance$data,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh
  )
}
