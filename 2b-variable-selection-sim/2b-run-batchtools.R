source(here::here("2b-variable-selection-sim/2b-simulation.R"))
source(here::here("2b-variable-selection-sim/2b-algorithms.R"))
source(here::here("R/utils.R"))

library(batchtools)
library(randomForestSRC)
library(CoxBoost)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_cache = FALSE,
  repls = 1000
)

set.seed(config$global_seed)
# Set to FALSE to remove + recreate registry
continue_bt <- FALSE

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "varsel-sim-pred"
#reg_name <- "DEBUG"
reg_dir <- here::here("registries", reg_name)

if (continue_bt) {
  loadRegistry(reg_dir, writeable = TRUE)
} else {
  unlink(reg_dir, recursive = TRUE)
  makeExperimentRegistry(
    file.dir = reg_dir,
    packages = c("cooper", "randomForestSRC", "CoxBoost", "survival", "riskRegression"),
    seed = config$global.seed,
    source = c(here::here("2b-variable-selection-sim/2b-algorithms.R"),
               here::here("2b-variable-selection-sim/2b-simulation.R"),
               here::here("R/utils.R"))
  )
}

# Problems -----------------------------------------------------------
addProblem(name = "sim_surv_binder", fun = sim_surv_binder, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "cooper", fun = cooper_varsel_wrapper)
addAlgorithm(name = "rfsrc", fun = rfsrc_varselect_wrapper)
addAlgorithm(name = "coxboost", fun = coxboost_varselect_wrapper)


# Experiments -----------------------------------------------------------
prob_design <- list(
  sim_surv_binder = expand.grid(
    n_train = 400,
    n_test = 400,
    p = 5000,
    ce = c(0.5),
    # cause 1 either slightly higher prevalence or noticeably lower
    lambda1 = c(0.1, 0.01),
    lambda2 = c(0.1),
    lambda_c = 0.1
  )
)

algo_design <- list(
  cooper = expand.grid(
    mt_max_iter = 3,
    alpha = c(1),
    t = c(100),
    thresh = c(1e-7)
  ),
  rfsrc = expand.grid(
    importance = "random",
    cutoff_method = "vita",
    splitrule = "logrank"
  ),
  coxboost = expand.grid(
    cmprsk = "csh" # c("sh", "csh", "ccsh"),
  )
)


addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
jobtbl <- unwrap(getJobPars(), c("algo.pars", "prob.pars"))
jobtbl[, chunk := chunk(job.id, chunk.size = 400, shuffle = TRUE)]

# Test jobs -----------------------------------------------------------
if (FALSE) testJob(id = 273)  # random cooper
if (FALSE) testJob(id = 785)  # random rfsrc
if (FALSE) testJob(id = 1171) # random coxboost
# Submit -----------------------------------------------------------

if (Sys.info()[["nodename"]] %in% c("blog1", "blog2")) {
  resources = list(
    partition = "teton-knl",
    memory = 4096,
    comment = "cooper-varsel",
    walltime = 3600 * 24 * 2
  )
} else {
  resources = list(comment = "cooper-varsel")
}

submitJobs(
  jobtbl,
  resources = resources
)
waitForJobs()

res <-  ijoin(reduceResultsDataTable(), flatten(getJobPars()))
jobs_total <- nrow(getJobTable())
jobs_completed <- 100 * nrow(res)/jobs_total
jobs_errored <- 100 * nrow(findErrors())/jobs_total

c(completed = jobs_completed, errored = jobs_errored)

saveRDS(res, here::here("results", "results-varsel.rds"))
