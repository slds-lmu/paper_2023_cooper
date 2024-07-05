source(here::here("R/utils.R"))
library(batchtools)
library(data.table)
library(ggplot2)

res_varsel <- readRDS(here::here("results", "2-results-varsel-csc-varsel.rds"))
res_perf <- readRDS(here::here("results", "2-results-varsel-csc-perf.rds"))


# vartiable selection -------------------------------------------------------------------------

res_varsel[, thresh := format(thresh, scientific = TRUE)]
res_varsel[, block := gsub("3", "3.", block)]
res_varsel[, block := factor(block, levels = sort(unique(block)))]

res_varsel[, lambda_setting := fcase(
  lambda1 == lambda2, "equal",
  lambda1 < lambda2,  "c1 less"
)]


res_varsel[total_pos == 0, tp := NA]
res_varsel[total_pos == 0, fn := NA]
res_varsel[, tpr := tp/total_pos]
res_varsel[, fpr := fp/total_neg]
res_varsel[, ppv := tp/(tp + fp)]
res_varsel[, npv := tn/(tn + fn)]
res_varsel[, fdr := 1 - ppv]
res_varsel[, f1 := (2 * ppv * tpr) / (ppv + tpr)]
res_varsel[, acc := (tp+tn)/total]


plot_varselect_boxplot(
  res_varsel,
  measure = "PPV",
  #model %in% c("cooper", "glmnet"),
  lambda_setting = "equal",
  blocks = c("block1", "block2", "block3.1", "block3.2")
)

# Plot scores ---------------------------------------------------------------------------------
aggr_perf = res_perf[, .(
  mean = mean(score),
  median = median(score),
  q25 = quantile(score, probs = .25),
  q75 = quantile(score, probs = .75),
  sd = sd(score)
  ),
  by = .(model, metric, cause, time_quant, lambda1, lambda2)]


library(dplyr)

(p_scores <- aggr_perf |>
  filter(lambda1 == lambda2) |>
  filter(metric %in% c("Brier", "AUC")) |>
  mutate(
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
  facet_grid(cols = vars(cause), rows = vars(metric), scales = "free_y", labeller = label_both) +
  geom_line() +
  # geom_errorbar(aes(ymin = q25, ymax = q75)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 1/7) +
  geom_point() +
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
  ))

ggsave(
  plot = p_scores,
  filename = fs::path(here::here("results"), "2-performance-scores", ext = ".png"),
  width = 7, height = 8, bg = "white"
)
