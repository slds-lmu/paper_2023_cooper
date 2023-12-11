set.seed(2023)

bladder <- readRDS(here::here("data/bladder-binder-clinical_geno.rds"))

table(bladder$status)

fit <- cooper::cooper(
  bladder,
  mt_max_iter = 3,
  alpha = 1,
  t = 100,
  thresh = 1e-7,
  stratify_by_status = TRUE,
  nfolds = 5
)

# Extract
glmnet_beta1 <- coef(fit, event = 1, use_initial_fit = TRUE)
glmnet_beta2 <- coef(fit, event = 2, use_initial_fit = TRUE)
cooper_beta1 <- coef(fit, event = 1)
cooper_beta2 <- coef(fit, event = 2)

selected <- list(
  cooper = list(
    cause1 = cooper_beta1[cooper_beta1 != 0],
    cause2 = cooper_beta2[cooper_beta2 != 0]
  ),
  glmnet = list(
    cause1 = glmnet_beta1[glmnet_beta1 != 0],
    cause2 = glmnet_beta2[glmnet_beta2 != 0]
  )
)

# Covariables shared between causes for
# cooper:
intersect(names(selected$cooper$cause1), names(selected$cooper$cause2))
# glmnet:
intersect(names(selected$glmnet$cause1), names(selected$glmnet$cause2))


