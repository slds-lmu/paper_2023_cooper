source(here::here("R/utils.R"))
library(batchtools)
library(data.table)
library(ggplot2)

res_varsel <- readRDS(here::here("results", "2-results-varsel-csc-varsel.rds"))
res_perf <- readRDS(here::here("results", "2-results-varsel-csc-perf.rds"))


# Variable selection --------------------------------------------------------------------------

res_varsel[, thresh := format(thresh, scientific = TRUE)]
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



res_varsel_long <- res_varsel |>
  mutate(
    cause = paste0("Cause ", cause),
    block = dplyr::case_when(
      block == "block1" ~ "B1 (Mutual)",
      block == "block2" ~ "B2 (Reversed)",
      block == "block3.1" ~ "B3 (Disjoint 1)",
      block == "block3.2" ~ "B3 (Disjoint 2)",
      block == "block4" ~ "B4 (Cor. Noise)",
      TRUE ~ "Noise"
    ),
    model = case_when(
      model == "cooper" ~ "CooPeR",
      model == "coxboost" ~ "CoxBoost",
      model == "glmnet" ~ "Coxnet",
      model == "rfsrc" ~ "RSF"
    ),
    model = factor(model, levels = rev(c("CooPeR", "Coxnet", "RSF", "CoxBoost")))
  ) |>
  filter(lambda2 == lambda1) |>
  tidyr::pivot_longer(cols = tpr:acc, names_to = "measure", values_to = "value", names_transform = toupper)


p_ppv1 <- res_varsel_long |>
  filter(measure %in% c("PPV")) |>
  filter(!is.na(value)) |>
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  facet_grid(rows = vars(cause), cols = vars(block), scales = "free") +
  coord_flip() +
  geom_boxplot(alpha = 0.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill"), guide = "none") +
  labs(
    x = NULL, y = "PPV [%]", color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14)

ggsave(plot = p_ppv1, filename = "2-ppv-equallambda.pdf", path = here::here("results"), width = 12, height = 4)
ggsave(plot = p_ppv1, filename = "2-ppv-equallambda.png", path = here::here("results"), width = 12, height = 4, bg = "white")

p_fpr1 <- res_varsel_long |>
  filter(measure %in% c("FPR")) |>
  filter(!is.na(value)) |>
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  facet_grid(rows = vars(cause), cols = vars(block), scales = "free") +
  coord_flip() +
  geom_boxplot(alpha = 0.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill"), guide = "none") +
  labs(
    x = NULL, y = "FPR [%]", color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14)

ggsave(plot = p_fpr1, filename = "2-fpr-equallambda.pdf", path = here::here("results"), width = 12, height = 4)
ggsave(plot = p_fpr1, filename = "2-fpr-equallambda.png", path = here::here("results"), width = 12, height = 4, bg = "white")

p_f1 <- res_varsel_long |>
  filter(measure %in% c("F1")) |>
  filter(!is.na(value)) |>
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  facet_grid(rows = vars(cause), cols = vars(block), scales = "free") +
  coord_flip() +
  geom_boxplot(alpha = 0.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill"), guide = "none") +
  labs(
    x = NULL, y = "F1 [%]", color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14)

ggsave(plot = p_fpr1, filename = "2-f1-equallambda.pdf", path = here::here("results"), width = 12, height = 4)
ggsave(plot = p_fpr1, filename = "2-f1-equallambda.png", path = here::here("results"), width = 12, height = 4, bg = "white")


library(kableExtra)

# table_base <- res_varsel_long |>
#   filter(measure == "TPR") |>
#   select(model, value, block, cause) |>
#   group_by(model, block, cause) |>
#   summarize(
#     mean = mean(value, na.rm = TRUE),
#     sd = sd(value, na.rm = TRUE),
#     median = median(value, na.rm = TRUE),
#     q25 = quantile(value, probs = 0.25, na.rm = TRUE),
#     q75 = quantile(value, probs = 0.75, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   mutate(across(where(is.numeric), \(x) {
#     x = round(100 * x, 2)
#     ifelse(is.finite(x), x, "-")
#   })) |>
#   mutate(
#     medianq = glue::glue("{median} [{q25}, {q75}]"),
#     meansd = glue::glue("{mean} ({sd})")
#   )
#
# table_base |>
#   tidyr::pivot_wider(id_cols = c("cause", "model"), names_from = c("block"), values_from = c("medianq")) |>
#   arrange(cause, desc(model)) |>
#   select(Cuase = cause, Model = model, matches("B[123]")) |>
#   kbl(
#     caption = "TPR across 1000 replications with median and 25%, 75% quantile",
#     format = "latex"
#   ) |>
#   kable_styling() |>
#   collapse_rows(1) |>
#   # writeLines()
#   save_kable(file = here::here("results", "tpr-table.pdf"), keep_tex = TRUE)



ppv_tbl_median <- measure_table(res_varsel_long, measure = "PPV", aggr_point = median, aggr_var = IQR, minmax = max)
fpr_tbl_median <- measure_table(res_varsel_long, measure = "FPR", aggr_point = median, aggr_var = IQR, minmax = min)
tpr_tbl_median <- measure_table(res_varsel_long, measure = "TPR", aggr_point = median, aggr_var = IQR, minmax = max)
npv_tbl_median <- measure_table(res_varsel_long, measure = "NPV", aggr_point = median, aggr_var = IQR, minmax = max)


ppv_kbl_median <- ppv_tbl_median |>
  select(-Cause) |>
  kbl(
    caption = "Median (IQR) of PPV scores (\\%) across different covariate blocks for causes 1 and 2.\\label{tab:ppv-median}",
    escape = FALSE, booktabs = TRUE, format = "latex"
  ) |>
  kable_styling(protect_latex = TRUE) |>
  footnote(general = "Method with highest scores within each cause/block", general_title = "Bold:") |>
  pack_rows(index = table(ppv_tbl_median$Cause))

fpr_kbl_median <- fpr_tbl_median |>
  select(-Cause) |>
  kbl(
    caption = "Median (IQR) of FPR scores (\\%) across different covariate blocks for causes 1 and 2.\\label{tab:fpr-median}",
    escape = FALSE, booktabs = TRUE, format = "latex"
  ) |>
  kable_styling(protect_latex = TRUE) |>
  footnote(general = "Method with lowest scores within each cause/block", general_title = "Bold:") |>
  pack_rows(index = table(fpr_tbl_median$Cause))

tpr_kbl_median <- tpr_tbl_median |>
  select(-Cause) |>
  kbl(
    caption = "Median (IQR) of TPR scores (\\%) across different covariate blocks for causes 1 and 2.\\label{tab:tpr-median}",
    escape = FALSE, booktabs = TRUE, format = "latex"
  ) |>
  kable_styling(protect_latex = TRUE) |>
  footnote(general = "Method with highest scores within each cause/block", general_title = "Bold:") |>
  pack_rows(index = table(tpr_tbl_median$Cause))

npv_kbl_median <- npv_tbl_median |>
  select(-Cause) |>
  kbl(
    caption = "Median (IQR) of NPV scores (\\%) across different covariate blocks for causes 1 and 2.\\label{tab:npv-median}",
    escape = FALSE, booktabs = TRUE, format = "latex"
  ) |>
  kable_styling(protect_latex = TRUE) |>
  footnote(general = "Method with highest scores within each cause/block", general_title = "Bold:") |>
  pack_rows(index = table(npv_tbl_median$Cause))


writeLines(ppv_kbl_median, con = here::here("results", "tab-ppv-median.tex"))
writeLines(fpr_kbl_median, con = here::here("results", "tab-fpr-median.tex"))
writeLines(tpr_kbl_median, con = here::here("results", "tab-tpr-median.tex"))
writeLines(npv_kbl_median, con = here::here("results", "tab-npv-median.tex"))


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
  # geom_errorbar(aes(ymin = q25, ymax = q75)) +
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
  ))

ggsave(
  plot = p_scores,
  filename = fs::path(here::here("results"), "2-performance-scores", ext = ".png"),
  width = 8, height = 6, bg = "white"
)
