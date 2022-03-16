#' Tidy up the output of fwelnet_mt_cox
#'
#'
#' @param mtres As returned by [`fwelnet_mt_cox()`].
#' @param truth A data.frame of true effects, with a structure like
#' |x  |beta  | truth|
#' |:--|:-----|-----:|
#' |x1 |beta1 |   0.5|
#' |x1 |beta2 |   0.5|
#' |x2 |beta2 |   1.0|
#' @return A data.frame in long format suitable for plotting.
tidy_mt_res <- function(mtres, truth) {
  betalist <- mtres[grepl("beta", names(mtres))]

  purrr::map2_df(names(betalist), betalist, ~{
    data.frame(.y) |>
      tibble::rownames_to_column("x") |>
      tidyr::pivot_longer(
        cols = -1,
        names_to = "iter",
        names_transform = ~stringr::str_extract(.x, "\\d+") |> as.integer()
      ) |>
      dplyr::mutate(
        beta = .x
      )
  }) |>
    dplyr::mutate(
      iter = iter - 1,
      is_noise = ifelse(x %in% unique(truth$x), "True effect", "Noise"),
      xcol = ifelse(is_noise == "Noise", "Noise", x) |>
        factor(levels = c(unique(truth$x), "Noise")),
      z_scale = mtres$z_scale,
      z_method = factor(mtres$z_method, levels = c("original", "aligned")),
      converged = mtres$converged
    )
}
