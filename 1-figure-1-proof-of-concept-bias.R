library(ggplot2)
library(dplyr)

figure_path <- here::here("results", "figures")
if (!file.exists(figure_path)) dir.create(figure_path)

res_long <- readRDS(here::here("results", "1-results-long-poc.rds"))

# Plot ----------------------------------------------------------------------------------------

p_poc <- res_long |>
  mutate(
    Setting = toupper(sub(x = problem, "sim_", "")),
    Cause = cause
  ) |>
  ggplot(aes(x = variable, y = error, fill = method, color = method)) +
  facet_grid(rows = vars(Cause), cols = vars(Setting), labeller = labeller(Cause = label_both, Setting = label_both), scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  geom_boxplot(alpha = 0.2) +
  scale_color_brewer(
    palette = "Dark2", aesthetics = c("color", "fill"),
    breaks = c("coxnet", "cooper"),
    labels = c(cooper = "CooPeR", coxnet = "Coxnet"),
    name = NULL
  ) +
  labs(
    title = NULL,
    x = NULL, y = "Error (true - estimated coefficient)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.spacing = ggplot2::unit(1, "cm"),
    legend.position = "bottom",
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(plot = p_poc, filename = "1-poc-boxplot-errors.pdf", path = figure_path, width = 10, height = 6)
ggsave(plot = p_poc, filename = "1-poc-boxplot-errors.png", path = figure_path, width = 10, height = 6, bg = "white")

