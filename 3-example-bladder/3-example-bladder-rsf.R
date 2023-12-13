library(data.table)
library(randomForestSRC)

set.seed(2023)
dat <- data.table::as.data.table(readRDS(here::here("data/bladder-binder-clinical_geno.rds")))

# RSF ---------------------------------------------------------------------
rf_c1 <- rfsrc(
  Surv(time, status) ~ ., data = dat,
  splitrule = "logrank",
  cause = c(1, 0),
  importance = "random",
  mtry = 2000,
  nodesize = 30
)

rf_c2 <- rfsrc(
  Surv(time, status) ~ ., data = dat,
  splitrule = "logrank",
  cause = c(0, 1),
  importance = "random",
  mtry = 2000,
  nodesize = 30
)

vimps <- data.table::data.table(
  variable = rownames(rf_c1[["importance"]]),
  vi_c1 = rf_c1[["importance"]][, 1],
  vi_c2 = rf_c2[["importance"]][, 2]
)

vimps[, vita_c1 := ifelse(vi_c1 <= abs(min(vi_c1)), 0, vi_c1)]
vimps[, vita_c2 := ifelse(vi_c2 <= abs(min(vi_c2)), 0, vi_c2)]

vimps[vita_c1 != 0 | vita_c2 != 0, ]

# Count selected per cause
nrow(vimps[vita_c1 != 0, ])
nrow(vimps[vita_c2 != 0, ])

# Overlap between causes?
vimps[vita_c1 != 0 & vita_c2 != 0, ]

