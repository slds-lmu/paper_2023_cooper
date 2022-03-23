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
tidy_mt_res <- function(mtres, truth, problem = NULL) {

  if (!is.null(problem)) truth <- truth[truth$problem == problem, ]

  betalist <- mtres[grepl("beta", names(mtres))]

  xdf <- purrr::map2_df(names(betalist), betalist, ~{
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
      converged = mtres$converged,
      theta_c1 =  as.numeric(mtres$fwfit1$glmfit$theta),
      theta_c2 =  as.numeric(mtres$fwfit2$glmfit$theta),
      lambda_c1 = mtres$fwfit1$lambda.min,
      lambda_c2 = mtres$fwfit2$lambda.min
    )

  if (!is.null(problem)) xdf$problem <- problem

  xdf

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

# line plot version of above, but for batchtoolsified output
plot_bt_res <- function(res, truth, problem) {

  truth <- truth[truth$problem == problem, ]
  res <- res[res$problem == problem, ]

  # Sort so Noise comes last in legend
  x_lvl <- c(sort(levels(res$xcol))[-1], "Noise")

  ggplot2::ggplot(res, ggplot2::aes(x = iter, y = value, color = xcol, alpha = is_noise)) +
    ggplot2::facet_grid(
      cols = vars(z_method, z_scale), rows = ggplot2::vars(beta),
      labeller = label_context
    ) +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::geom_hline(data = truth, ggplot2::aes(yintercept = truth), lty = "dashed") +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 1)) +
    ggplot2::scale_color_brewer(palette = "Dark2", breaks = x_lvl) +
    ggplot2::scale_alpha_manual(
      values = c("True effect" = 0.5, "Noise" = 0.05),
      guide = "none"
    ) +
    ggplot2::labs(
      title = glue::glue("Multi-Task fwelnet: Setting {problem}"),
      x = "# of Multi-Task Iterations",
      y = "Effect estimate",
      color = "Variable"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.spacing.y = ggplot2::unit(1, "cm"),
      legend.position = "top",
      plot.title.position = "plot"
    )
}

