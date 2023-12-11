
library(data.table)
library(CoxBoost)

set.seed(2023)

bladder_file <- here::here("data/bladder-binder-clinical_geno.rds")
dat <- data.table::as.data.table(readRDS(bladder_file))


# Coxboost ----------------------------------------------------------------
cbfit <- CoxBoost(
  time = dat$time,
  status = dat$status,
  x = as.matrix(dat[, -c(1, 2)]),
  cmprsk = "csh",
  stepno = 100,
  penalty = 2000
)
cb_coefs <- coef(cbfit)
coefs_c1 <- cb_coefs[[1]][cb_coefs[[1]] != 0]
coefs_c2 <- cb_coefs[[2]][cb_coefs[[2]] != 0]

coefs_c1
coefs_c2

# Count selected per cause
length(coefs_c1)
length(coefs_c2)

# Overlap between causes?
intersect(names(coefs_c1), names(coefs_c2))

