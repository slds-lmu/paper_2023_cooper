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
(cooper_shared <- intersect(names(selected$cooper$cause1), names(selected$cooper$cause2)))
# Coxnet:
(coxnet_shared <- intersect(names(selected$coxnet$cause1), names(selected$coxnet$cause2)))

coxnet_beta1[names(coxnet_beta1) == "age"]
coxnet_beta2[names(coxnet_beta2) == "age"]
cooper_beta1[names(cooper_beta1) == "age"]
cooper_beta2[names(cooper_beta2) == "age"]

reference <- readxl::read_excel("data-raw/10780432ccr062940-sup-supplemental_file_2.xls", sheet = "Progression classifier probes") |>
  janitor::clean_names()

cooper_shared[cooper_shared %in% reference$probe_id]

reference |>
  dplyr::filter(probe_id %in% cooper_shared)
