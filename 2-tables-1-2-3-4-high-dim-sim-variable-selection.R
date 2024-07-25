source(here::here("R/tables.R"))
library(data.table)
library(dplyr)
library(kableExtra)

table_path <- here::here("results", "tables")
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

# Tables ----------------------------------------------------------------------
# Produces all tables in the manuscript, all located in the appendix

for (measure in c("PPV", "FPR", "TPR", "NPV")) {
  tab <- measure_table(res_varsel_long, measure = measure, aggr_point = median, aggr_var = IQR, minmax = max)

  tab_kbl <- tab |>
    select(-Cause) |>
    kbl(
      caption = sprintf(
        "Median (IQR) of %s scores (\\%%) across different covariate blocks for causes 1 and 2.\\label{tab:%s-median}",
        measure, tolower(measure)),
      escape = FALSE, booktabs = TRUE, format = "latex"
    ) |>
    kable_styling(protect_latex = TRUE) |>
    footnote(general = "Method with best scores within each cause/block", general_title = "Bold:") |>
    pack_rows(index = table(tab$Cause))

  writeLines(tab_kbl, con = fs::path(table_path, sprintf("tab-%s-median", tolower(measure)), ext = "tex"))
}
