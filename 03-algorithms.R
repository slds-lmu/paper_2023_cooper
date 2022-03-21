# Multi-task fwelnet-----------------------------------------------------------
fwel_wrapper <- function(
    data, job, instance,
    alpha = 1, z_scale = 1, z_method = "original",
    mt_max_iter = 2) {

  fwelnet::fwelnet_mt_cox(
    instance$data, mt_max_iter = mt_max_iter,
    alpha = alpha,
    z_scale = z_scale, z_method = z_method
  )
}
