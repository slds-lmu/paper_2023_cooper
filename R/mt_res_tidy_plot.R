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

#' Plot the tidied up results for a single run
#' @param mtres As returned by [`tidy_mt_res()`].
#' @param truth A data.frame of true effects as in [`tidy_mt_res`]
plot_mt_res <- function(mtres, true_effects) {
  ggplot2::ggplot(mtres, ggplot2::aes(x = iter, y = value, color = xcol, alpha = is_noise)) +
    ggplot2::facet_wrap(ggplot2::vars(beta)) +
    ggplot2::geom_hline(data = true_effects,
                        ggplot2::aes(yintercept = truth), lty = "dotted") +
    ggplot2::geom_path() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_alpha_manual(
      values = c("True effect" = 1, "Noise" = 0.2),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "Multi-Task fwelnet",
      subtitle = glue::glue(
        "z scalar: {mtres$z_scale}",
        "z method: {mtres$z_method}.",
        "Conerveged: {tolower(mtres$converged)}", .sep = "\n"
      ),
      x = "# of Multi-Task Iterations",
      y = "Effect estimate",
      color = "Variable"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", plot.title.position = "plot")
}
