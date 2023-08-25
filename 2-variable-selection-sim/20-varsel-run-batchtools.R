library(batchtools)
library(randomForestSRC)
library(CoxBoost)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_cache = FALSE,
  repls = 100
)

set.seed(config$global_seed)
# Set to FALSE to remove + recreate registry
continue_bt <- TRUE

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "fwel_sim_varsel_2309"
reg_dir <- here::here("registries", reg_name)

if (continue_bt) {
  loadRegistry(reg_dir, writeable = TRUE)
} else {
  unlink(reg_dir, recursive = TRUE)
  makeExperimentRegistry(file.dir = reg_dir, packages = c("randomForestSRC", "CoxBoost"))

}

# Problems -----------------------------------------------------------
addProblem(name = "sim_surv_binder", fun = sim_surv_binder, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "fwel_mt", fun = fwel_mt_varselect_wrapper)
addAlgorithm(name = "rfsrc", fun = rfsrc_varselect_wrapper)
addAlgorithm(name = "coxboost", fun = coxboost_varselect_wrapper)


# Experiments -----------------------------------------------------------
prob_design <- list(
  sim_surv_binder = expand.grid(
    n_train = 400,
    p = 5000,
    ce = 0.5,
    lambda1 = c(0.1, 0.01),
    lambda2 = c(0.1),
    lambda_c = 0.1
  )
)

algo_design <- list(
  fwel_mt = expand.grid(
    mt_max_iter = 3,
    alpha = c(1),
    t = c(100),
    thresh = c(1e-7)
  ),
  rfsrc = expand.grid(
    importance = "random",
    cutoff_method = "vita",
    mtry = c(500, 2000, 3000),
    nodesize = 30,
    splitrule = "logrank"
  ),
  coxboost = expand.grid(
    cmprsk = c("sh", "csh", "ccsh"),
    stepno = c(100, 300),
    penalty = c(1000, 2000, 3000)
  )
)


addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
jobtbl <- unwrap(getJobPars(), c("algo.pars", "prob.pars"))

# Test jobs -----------------------------------------------------------
if (interactive()) testJob(id = 200)

# Submit -----------------------------------------------------------

submitJobs(jobtbl[algorithm == "coxboost"])
submitJobs(jobtbl[algorithm == "fwel_mt"])
submitJobs(jobtbl[algorithm == "rfsrc"])

# ids <- findNotStarted()
# submitJobs(ids)



# Monitor jobs ------------------------------------------------------------
if (interactive()) {
  getStatus()

  ijoin(
    unwrap(getJobPars(findErrors()), c("algo.pars", "prob.pars")),
    getErrorMessages(findErrors())
  )
}

# Get results -------------------------------------------------------------
# res <-  ijoin(reduceResultsDataTable(), flatten(getJobPars()))
# res

