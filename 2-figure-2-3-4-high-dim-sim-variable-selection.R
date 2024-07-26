library(dplyr)
library(ggplot2)

figure_path <- here::here("results", "figures")
table_path <- here::here("results", "tables")
if (!file.exists(figure_path)) dir.create(figure_path)
if (!file.exists(table_path)) dir.create(table_path)

# Postprocess results ---------------------------------------------------------
res_varsel <- readRDS(here::here("results", "2-results-varsel-csc-varsel.rds"))

res_varsel_long <- res_varsel |>
  mutate(
    cause = paste0("Cause ", cause),
    block = case_when(
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


# Plots ---------------------------------------------------------------------------------------
for (measure in c("PPV", "FPR", "F1")) {

  p <- res_varsel_long |>
    filter(.data$measure %in% .env$measure) |>
    filter(!is.na(value)) |>
    ggplot(aes(x = model, y = value, color = model, fill = model)) +
    facet_grid(rows = vars(cause), cols = vars(block), scales = "free") +
    coord_flip() +
    geom_boxplot(alpha = 0.5) +
    #scale_y_continuous(labels = scales::percent_format()) +
    scale_y_continuous(labels = \(x) 100 * x) +
    scale_color_brewer(palette = "Dark2", aesthetics = c("color", "fill"), guide = "none") +
    labs(
      x = NULL, y = glue::glue("{measure} [%]"), color = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  ggsave(plot = p, filename = glue::glue("2-{tolower(measure)}-equallambda.pdf"), path = figure_path, width = 12, height = 4)
  ggsave(plot = p, filename = glue::glue("2-{tolower(measure)}-equallambda.png"), path = figure_path, width = 12, height = 4, bg = "white")

}
