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

#' Hacky version of the above to work on the batchtools output
#' I'm so sorry
tidy_mt_res_bt <- function(res, truth) {
  purrr::pmap_dfr(res, ~ {
    job_id <- ..1
    mtres <- ..2
    problem <- ..3

    if (!is.null(problem)) truth <- truth[truth$problem == problem, ]

    betalist <- mtres[grepl("beta", names(mtres))]

    res_tidy <- purrr::map2_df(names(betalist), betalist, ~ {
      data.frame(.y) |>
        tibble::rownames_to_column("x") |>
        tidyr::pivot_longer(
          cols = -1,
          names_to = "iter",
          names_transform = ~ stringr::str_extract(.x, "\\d+") |> as.integer()
        ) |>
        dplyr::mutate(
          beta = .x
        )
    }) |>
      dplyr::mutate(
        job.id = job_id,
        iter = iter - 1,
        is_noise = ifelse(x %in% unique(truth$x), "True effect", "Noise"),
        xcol = ifelse(is_noise == "Noise", "Noise", x) |>
          factor(levels = c(unique(truth$x), "Noise")),
        converged = mtres$converged,
        theta_c1 = as.numeric(mtres$fwfit1$glmfit$theta),
        theta_c2 = as.numeric(mtres$fwfit2$glmfit$theta),
        lambda_c1 = mtres$fwfit1$lambda.min,
        lambda_c2 = mtres$fwfit2$lambda.min
      )

    dplyr::left_join(res_tidy, res[, -2], by = "job.id")
  })
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


#' Box/lineplots of effect estimates for each mt iteration step based on longified batchtools
#' output
#' @param res_long As returned by `tidy_mt_res_bt`.
#' @param true_effects Table of true effects as defined in `02-problems.R`.
#' @param problem Setting to supset by, e.g. `"sim_a"`. Vectorized, multiple
#'   plots will be added together via `patchwork`.
#' @param exclude_noise `[FALSE]` If `TRUE`, only effect variables are shown
lineplot_bt_res <- function(res_long, true_effects, problem, exclude_noise = FALSE, ...) {
  # shoddy vectorization with auto-patchworking
  if (length(problem) > 1) {
    p_list <- purrr::map(
      problem, ~lineplot_bt_res(
        res_long = res_long,
        true_effects = true_effects, problem = .x,
        exclude_noise = exclude_noise
      )
    )

    p <- Reduce(`+`, p_list) +
      patchwork::plot_layout(ncol = 2, guides = "collect") &
      ggplot2::theme(legend.position = "top")

    return(p)
  }

  xdf <- res_long |>
    dplyr::filter(problem == !!problem, z_method == "original", ...)

  if (exclude_noise) {
    xdf <- xdf |> dplyr::filter(is_noise != "Noise")
    x_lvl <- sort(levels(xdf$xcol))
  } else {
    # Sort so Noise comes last in legend
    x_lvl <- c(sort(levels(xdf$xcol))[-1], "Noise")
  }


  xdf |>
    dplyr::mutate(beta = stringr::str_replace_all(beta, "beta", "Cause ")) |>
    ggplot2::ggplot(ggplot2::aes(x = iter, y = value, color = xcol, alpha = is_noise, group = paste0(xcol, job.id))) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(z_scale, theta, t, thresh), rows = ggplot2::vars(beta),
      labeller = label_context, scales = "free"
    ) +
    # stat_summary(aes(group = x), geom = "point", fun = mean, size = 2) +
    ggplot2::geom_boxplot(aes(x = factor(iter), group = NULL)) +
    # ggplot2::geom_path(size = 1.5) +
    # ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(
      data = dplyr::filter(true_effects, problem == !!problem) |>
        dplyr::mutate(beta = stringr::str_replace_all(beta, "beta", "Cause ")),
      ggplot2::aes(yintercept = truth),
      lty = "dashed"
    ) +
    #ggplot2::scale_x_continuous(breaks = seq(0, 100, 1)) +
    ggplot2::scale_y_continuous(breaks = seq(-1, 2, .25)) +
    ggplot2::scale_color_brewer(palette = "Dark2", breaks = x_lvl) +
    ggplot2::scale_alpha_manual(
      values = c("True effect" = 0.75, "Noise" = 0.05),
      guide = "none"
    ) +
    ggplot2::labs(
      title = sim_labels[[problem]],
      x = "# of Multi-Task Iterations",
      y = "Effect estimate",
      color = "Variable"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(1, "cm"),
      legend.position = "top",
      plot.title.position = "plot"
    )
}
