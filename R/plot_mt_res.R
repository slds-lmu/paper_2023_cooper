plot_mt_res <- function(mtres, true_effects) {
  ggplot2::ggplot(mtres, aes(x = iter, y = value, color = xcol, alpha = is_noise)) +
    ggplot2::facet_wrap(vars(beta)) +
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
