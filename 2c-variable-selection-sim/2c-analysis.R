source(here::here("R/utils.R"))
library(batchtools)
library(data.table)

res <- readRDS(here::here("results", "2-results-varsel.rds"))

# res$result[[1]]$varsel
# res$result[[1]]$scores

res_varsel <- rbindlist(lapply(res$job.id, \(id) {
  varsel_dt = as.data.table(res[job.id == id, result][[1]][["varsel"]])
  varsel_dt[, job.id := id]
  varsel_dt
}))

res_perf <- rbindlist(lapply(res$job.id, \(id) {
  perf_dt = as.data.table(res[job.id == id, result][[1]][["scores"]])
  perf_dt[, job.id := id]
  perf_dt
}))

res_varsel <- ljoin(res_varsel, pars)
res_perf <- ljoin(res_perf, pars)


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
  measure = "FNR",
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

p_scores <- aggr_perf |>
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
  ggplot(aes(x = 100 * time_quant, y = 100 * mean, colour = model, fill = model)) +
  facet_grid(cols = vars(cause), rows = vars(metric), scales = "free_y", labeller = label_both) +
  geom_line() +
  geom_point() +
  labs(
    title = "Performance scores from variable selection simulation",
    subtitle = "Based on 1000 replications with 400 observations for train and test sets",
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
  filename = fs::path(here::here("results"), "2-performance-scores", ext = ".png"),
  width = 7, height = 8, bg = "white"
)
