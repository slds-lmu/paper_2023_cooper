
library(dplyr)
library(ggplot2)

figure_path <- here::here("results", "figures")
if (!file.exists(figure_path)) dir.create(figure_path)

# Read split bladder data
scores_cmb <- readRDS(here::here("results/3-bladder-scores.rds"))

p = scores_cmb |>
  mutate(
    Cause = cause, Metric = metric,
    Metric = ifelse(Metric == "Brier", "Brier Score", "AUC"),
    model = dplyr::case_when(
      model == "cooper" ~ "CooPeR",
      model == "coxboost" ~ "CoxBoost",
      model == "coxnet" ~ "Coxnet",
      model == "rfsrc" ~ "RSF",
      model == "Null model" ~ "Null Model"
    ),
    model = factor(model, levels = rev(c("CooPeR", "Coxnet", "RSF", "CoxBoost", "Null Model")))
  ) |>
  filter(metric %in% c("Brier", "AUC")) |>
  ggplot(aes(x = 100 * time_quant, y = 100 * score, color = model, fill = model)) +
  facet_grid(cols = vars(Cause), rows = vars(Metric), scales = "free_y", labeller = label_both) +
  geom_line() +
  geom_point() +
    scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill")) +
  labs(
    title = "Bladder cancer: Performance of CSC fit with selected variables",
    subtitle = "Evaluation based on 70/30 train/test split",
    x = "Time quantile (%)", y = "Score (%)",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot"
  )

ggsave(
  plot = p,
  filename = fs::path(figure_path, "3-bladder-performance", ext = "png"),
  width = 8, height = 5, bg = "white"
)

ggsave(
  plot = p,
  filename = fs::path(figure_path, "3-bladder-performance", ext = "pdf"),
  width = 8, height = 5
)
