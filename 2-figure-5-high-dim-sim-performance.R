source(here::here("R/utils.R"))
library(data.table)
library(dplyr)
library(ggplot2)

figure_path <- here::here("results", "figures")
if (!file.exists(figure_path)) dir.create(figure_path)

# Postprocess results ---------------------------------------------------------
res_perf <- readRDS(here::here("results", "2-results-varsel-csc-perf.rds"))

aggr_perf = res_perf[, .(
  mean = mean(score),
  median = median(score),
  q25 = quantile(score, probs = .25),
  q75 = quantile(score, probs = .75),
  sd = sd(score)
),
by = .(model, metric, cause, time_quant, lambda1, lambda2)]

p_scores <- aggr_perf |>
  filter(lambda1 == lambda2) |>
  filter(metric %in% c("Brier", "AUC")) |>
  mutate(
    Cause = cause, Metric = metric,
    Metric = ifelse(Metric == "Brier", "Brier Score", "AUC"),
    model = dplyr::case_when(
      model == "cooper" ~ "CooPeR",
      model == "coxboost" ~ "CoxBoost",
      model == "glmnet" ~ "Coxnet",
      model == "rfsrc" ~ "RSF",
      model == "Null model" ~ "Null Model"
    ),
    model = factor(model, levels = rev(c("CooPeR", "Coxnet", "RSF", "CoxBoost", "Null Model")))
  ) |>
  ggplot(aes(x = 100 * time_quant, y = mean, colour = model, fill = model)) +
  facet_grid(cols = vars(Cause), rows = vars(Metric), scales = "free_y", labeller = label_both) +
  geom_line() +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 1/7) +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2", aesthetics = c("color", "fill")) +
  labs(
    title = "Performance scores from variable selection simulation",
    subtitle = paste0(
      "Based on 1000 replications with 400 observations for train and test sets\n",
      "Median with 25th and 75th percentiles shown as ribbons"),
    x = "Time quantile (%)",
    y = "Mean score (%)",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot"
  )

ggsave(
  plot = p_scores,
  filename = fs::path(figure_path, "2-performance-scores", ext = "png"),
  width = 8, height = 6, bg = "white"
)

ggsave(
  plot = p_scores,
  filename = fs::path(figure_path, "2-performance-scores", ext = "pdf"),
  width = 8, height = 6
)
