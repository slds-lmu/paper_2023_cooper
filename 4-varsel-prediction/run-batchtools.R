source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))

library(batchtools)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
library(rlang)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_cache = FALSE,
  repls = 200
)

set.seed(config$global_seed)
# Set to FALSE to remove + recreate registry
continue_bt <- FALSE

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "fwel_sim_varsel_predict"
reg_dir <- here::here("registries", reg_name)

if (continue_bt) {
  loadRegistry(reg_dir, writeable = TRUE)
} else {
  unlink(reg_dir, recursive = TRUE)
  makeExperimentRegistry(file.dir = reg_dir, packages = c("randomForestSRC", "CoxBoost", "rlang", "data.table", "riskRegression"), seed = config$global.seed)
}

# Problems -----------------------------------------------------------
addProblem(name = "bladder_geno", fun = get_bladder_data, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "fwel_mt", fun = fwel_mt_varselect_pred)
addAlgorithm(name = "rfsrc", fun = rfsrc_varselect_pred)
addAlgorithm(name = "coxboost", fun = coxboost_varselect_pred)


# Experiments -----------------------------------------------------------
prob_design <- list(
  bladder_geno = expand.grid(
    split = 2/3, standardize = TRUE
  )
)

algo_design <- list(
  fwel_mt = expand.grid(
    mt_max_iter = 3,
    alpha = 1,
    t = 100,
    thresh = 1e-7
  ),
  rfsrc = expand.grid(
    importance = "random",
    cutoff_method = "vita",
    mtry = c(500, 2000, 3000),
    nodesize = 30,
    splitrule = "logrank"
  ),
  coxboost = expand.grid(
    cmprsk = c("sh", "csh"),
    stepno = c(100, 300),
    penalty = c(1000, 3000)
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

submitJobs(c(1, 456, 23, 156))


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

