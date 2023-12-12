set.seed(2023)
bladder <- readRDS(here::here("data/bladder-binder-clinical_geno.rds"))

# Fit CooPeR
fit <- cooper::cooper(
  bladder, mt_max_iter = 3,
  alpha = 1, t = 100, thresh = 1e-7,
  stratify_by_status = TRUE, nfolds = 5
)

# Extract coefficients from initial Coxnet
coxnet_beta1 <- coef(fit, event = 1, use_initial_fit = TRUE)
coxnet_beta2 <- coef(fit, event = 2, use_initial_fit = TRUE)
# Extract final coefficients from CooPeR
cooper_beta1 <- coef(fit, event = 1)
cooper_beta2 <- coef(fit, event = 2)

selected <- list(
  cooper = list(
    cause1 = cooper_beta1[cooper_beta1 != 0],
    cause2 = cooper_beta2[cooper_beta2 != 0]
  ),
  coxnet = list(
    cause1 = coxnet_beta1[coxnet_beta1 != 0],
    cause2 = coxnet_beta2[coxnet_beta2 != 0]
  )
)

# Covariables shared between causes for
# CooPeR:
intersect(names(selected$cooper$cause1), names(selected$cooper$cause2)) |> sort()
# Coxnet:
intersect(names(selected$coxnet$cause1), names(selected$coxnet$cause2))


