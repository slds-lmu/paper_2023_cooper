
library(data.table)

set.seed(2023)

bladder_file <- here::here("data/bladder-binder-clinical_geno.rds")
dat <- data.table::as.data.table(readRDS(bladder_file))

mt_max_iter <- 3
fit <- fwelnet::fwelnet_mt_cox(
  dat,
  mt_max_iter = mt_max_iter,
  z_method = "original",
  alpha = 1,
  t = 100,
  a = 0.5,
  thresh = 1e-7,
  stratify_by_status = TRUE,
  nfolds = 5,
  include_mt_beta_history = TRUE
)

# Coefs
glmnet_beta1 <- fit$beta1[, 1]
glmnet_beta2 <- fit$beta2[, 1]
fwel_beta1 <- fit$beta1[, mt_max_iter + 1]
fwel_beta2 <- fit$beta2[, mt_max_iter + 1]

fw_coefs <- list(
  fwelnet = list(
    cause1 = fwel_beta1[fwel_beta1 != 0],
    cause2 = fwel_beta2[fwel_beta2 != 0]
  ),
  glmnet = list(
    cause1 = glmnet_beta1[glmnet_beta1 != 0],
    cause2 = glmnet_beta2[glmnet_beta2 != 0]
  )
)

# Same in causes?
intersect(names(fw_coefs$fwelnet$cause1), names(fw_coefs$fwelnet$cause2))
intersect(names(fw_coefs$glmnet$cause1), names(fw_coefs$glmnet$cause2))


