# Get results ---------------------------------------------------------------------------------
# source(here::here("1-proof-of-concept/1-problems.R"))

library(batchtools)
library(ggplot2)
library(dplyr)
# res <- readRDS(here::here("results", "1-results-poc.rds"))
res_long <- readRDS(here::here("results", "1-results-long-poc.rds"))


# sim_labels <- c(
#   "sim_a" = "A: x1 has equal effect on both causes / same prevalence",
#   "sim_b" = "B: x1 has effect on cause 1, x2 has equal effect on cause 2 / same prevalence",
#   "sim_c" = "C: x1 has effect on cause 1 and smaller effect on cause2 / cause 2 less prevalent",
#   "sim_d" = "D: x{1,2,3} have equal effects in both causes / cause 2 less prevalent"
# )

# Plot ----------------------------------------------------------------------------------------

(p_poc <- res_long |>
  mutate(
    Setting = toupper(sub(x = problem, "sim_", "")),
    Cause = cause
  ) |>
  ggplot(aes(x = variable, y = error, fill = method, color = method)) +
  facet_grid(rows = vars(Cause), cols = vars(Setting), labeller = labeller(Cause = label_both, Setting = label_both), scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 1) +
  #geom_jitter(alpha = 0.1, position = position_dodge2(width = .5, preserve = "single")) +
  geom_boxplot(alpha = 0.2) +
  scale_color_brewer(
    palette = "Dark2", aesthetics = c("color", "fill"),
    breaks = c("coxnet", "cooper"),
    labels = c(cooper = "CooPeR", coxnet = "Coxnet"),
    name = NULL
  ) +
  labs(
    title = NULL,
    #subtitle = sim_labels[[problem]],
    x = NULL, y = "Error (true - estimated coefficient)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.spacing = ggplot2::unit(1, "cm"),
    legend.position = "bottom",
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ))

ggsave(plot = p_poc, filename = "1-poc-boxplot-errors.pdf", path = here::here("results"), width = 10, height = 6)
ggsave(plot = p_poc, filename = "1-poc-boxplot-errors.png", path = here::here("results"), width = 10, height = 6)

